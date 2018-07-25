/* Python module for fast sampling from a finite population, with probability
 * proportional to size (PPS), with and without replacement.
 *
 * For sampling *with* replacement, we use
 * Marsaglia et al.'s compact 5-table method, as described in their paper 
 * 'Fast generation of discrete random variables' in the Journal of Statistical
 * Software, July 2004, vol 11, issue 3.  The implementation here is based
 * on Ed Schofield's IntSampler for fast sampling with replacement, and is
 * in the file sampler5tbl.c (copyright Ed Schofield, 2005-2006, as is
 * much of the interface code in this file).  This code is based
 * upon the C implementation that accompanies the JSS paper, but is simpler, 
 * and uses a different random number generator, the Mersenne Twister in
 * Jean-Sebastien Roy's RandomKit.
 *
 * For sampling *without* replacement, we use Sampford's rejection algorithm,
 * using the 5-table method to accelerate the rejection steps.
 *
 *
 * Modified from IntSampler, June 2006 by Tom Loredo for the Inference package:
 *   Removed use of NumPy's copy of RandomKit (not cross-platform); now
 *       use our own copy, NumPy's RandomKit state is explicitly passed in 
 *       and out as needed.
 *   Added sampling *without* replacement capability using Sampford's
 *       algorithm, using Ed's implementation of the Marsaglia et al.
 *       5-table algorithm to accelerate the rejection part of
 *       Sampford's algorithm.
 *   Added module-level functions; some maintain and allow testing of
 *       the copy of RandomKit, others enable unweighted sampling
 *       without replacement.
 *
 * TODO:
 *   Convert equalprob/i to use RandomKit.
 *   Make the 5-table version of SampfordPPS.
 *
 * License: BSD-style (see LICENSE.txt at root of scipy tree)
 *
 * Modifications:
 *   2007-07-07:  LONG -> ULONG in rng state get/set, since current numpy
 *                refuses to cast ulong to long
 */


/*******************************************************************************
* Headers, prototypes, globals, and macros.
*******************************************************************************/

#include "Python.h"
#include "numpy/arrayobject.h"
#include "sampford.h"
#include <string.h>

/* Error handler macro for use in any extension; adapted from Dubois & Yang. 
   This requires that the function using it have the FAIL label defined. */
#define Py_TRY(BOOLEAN) {if (!(BOOLEAN)) goto FAIL;}

/* Macros to cast numpy array data to a 1-d arrays; borrowed
   from lapack_litemodule.c in numpy. */
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

// Function prototypes:
static PyArrayObject* PyArray_Intsample(Sampler* mysampler, unsigned long size);
static PyObject* sample(PyObject *self, PyObject *args, PyObject *keywords);

/*------------------------------------------------------------------------------
* Local storage of the state of the PRNG, to be copied from numpy
* to ensure reproducibility when the user uses random.set_state to
* repeat runs.
------------------------------------------------------------------------------*/
static rk_state rng_state;


/*******************************************************************************
* Support for the new type.
*******************************************************************************/

// PPSSampler type
typedef struct {
    PyObject_HEAD
    Sampler* pSampler;
    Sampford* pSampford;
    Sampler* pRatioSampler;  /* For 5-table acceleration of Sampford */
} PPSSampler;

//    "destroy" should be called automatically upon deletion
static void PPSSampler_destroy(PPSSampler* self)
{
    // printf("[Started destroying sampler]\n");
    if (self->pSampler != NULL)
        destroy_sampler5tbl(self->pSampler);
    if (self->pSampford != NULL)
        destroy_sampford(self->pSampford);
    self->ob_type->tp_free((PyObject*)self);
    // printf("[Finished destroying sampler]\n");
}

static PyObject *
PPSSampler_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PPSSampler *self;

    // printf("[Started new() method]\n");
    self = (PPSSampler *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->pSampler = NULL;
        self->pSampford = NULL;
        self->pRatioSampler = NULL;
    }
    // printf("[Finished new() method]\n");
    return (PyObject *) self;
}

static int
PPSSampler_init(PPSSampler *self, PyObject *args, PyObject *kwds)
{
    PyArrayObject *pmf_table;
    int k;                              /* size of the table */

    // printf("[Started initializing sampler]\n");

    /* parse the arguments  */
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &pmf_table)) {
      return -1;  /* Error indicator */
    }

    // printf("[Parsed arguments]\n");

    /* check that we have a one-dimensional array */
    if (pmf_table->nd != 1) {
      /* throw a Python exception (ValueError) */
      PyErr_SetString(PyExc_ValueError,
                      "the weight table must be a one-dimensional array");
      return -1;
    }

    // printf("[Array is 1D]\n");

    /* check that the data type is float64, (C double) */
    if (PyArray_DESCR(pmf_table)->type_num != NPY_DOUBLE) {
      PyErr_SetString(PyExc_ValueError,
                      "the pmf table must be of type float64");
      return -1;
    }

    // printf("[Data type is float64]\n");

    k = pmf_table->dimensions[0];    /* length of the array */

    self->pSampler = init_sampler5tbl((double*) pmf_table->data, k);
    if (self->pSampler == NULL)
    {
        Py_DECREF(self);
        return -1;
    }
    
    // printf("[Finished initializing sampler]\n");

    // Eventually, we need to do this:
    // destroy(self->pSampler);

    return 0;
}


static char sample__doc__[] = \
  "sample(size): return an array with a random discrete sample\n"\
  "of the given size from the probability mass function specified when\n"\
  "initializing the sampler.\n";
                                                      
static PyObject*
sample(PyObject *self, PyObject *args, PyObject *keywords)
{
    PyArrayObject *samplearray;

    int size;
    
    /* parse the arguments  */
    if (!PyArg_ParseTuple(args, "i", &size))
        return NULL;

    // printf("[Parsed arguments]\n");

    // Check that size > 0
    if (size <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "sample size must be positive");
        return NULL;
    }

    // printf("[size is > 0]\n");

    samplearray = PyArray_Intsample(((PPSSampler*)self)->pSampler, size);

    return (PyObject*) samplearray;
}



static PyArrayObject*
PyArray_Intsample(Sampler* mysampler, unsigned long size)
{
    PyArrayObject* samplearray;
    unsigned long* ptr;
    
    int ndim = 1;
    npy_intp dims[1] = {size};
    samplearray = (PyArrayObject*) PyArray_SimpleNew(ndim, dims, NPY_LONG);
    if (samplearray == NULL) return NULL;

    ptr = (unsigned long*) PyArray_DATA(samplearray);
    Dran_array(mysampler, ptr, size, &rng_state);

    return samplearray;
}


/*------------------------------------------------------------------------------
* Prepare weights for PPS sampling (sort, accumulate).
------------------------------------------------------------------------------*/
static char prepwts_doc[] =
"prepwts(wts):  Prepare weights for PPS sampling (sort, accumulate)." ;

static PyObject *
prepwts (PyObject *self, PyObject *args, PyObject *keywords) {
    int npopn;
    PyObject *win;
    PyArrayObject *wts, *sorted, *cumwt, *indices;
    PyArrayObject *no_out=NULL;
    npy_intp dims[1];
    int *idata;
    int i, j, itemp;
    double *wdata, *sdata, tot;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    wts = NULL;
    sorted = NULL;
    cumwt = NULL;
    indices = NULL;

/* Parse the input weights. */
    Py_TRY(PyArg_ParseTuple(args, "O", &win));
    wts = (PyArrayObject *)
            PyArray_ContiguousFromAny(win, NPY_DOUBLE, 1, 1);
    if (wts == NULL) goto FAIL;
    npopn = wts->dimensions[0];
    dims[0] = npopn;
    sorted = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (sorted == NULL) goto FAIL;

/* Sort and reverse the weights. */
	indices = (PyArrayObject *) PyArray_ArgSort(wts, 0, PyArray_QUICKSORT);
	idata = IDATA(indices);

/*
	if (1==1) return Py_BuildValue("NN", PyArray_Return(indices), PyArray_Return(cumwt));
*/

	for (i=0; i<npopn/2; i++) {
		j = npopn-i-1;
		itemp = idata[i];
		idata[i] = idata[j];
		idata[j] = itemp;
	}
	wdata = DDATA(wts);
	sdata = DDATA(sorted);
	for (i=0; i<npopn; i++) {
		sdata[i] = wdata[idata[i]];
	}

/* Accumulate the reverse-sorted weights; save total weight. */
	cumwt = (PyArrayObject *) PyArray_CumSum(sorted, 0, NPY_DOUBLE, no_out);
	wdata = DDATA(cumwt);
	tot = wdata[npopn-1];

/* Copy the info to this instance. */
    ((PPSSampler*)self)->pSampford = init_sampford(npopn, DDATA(cumwt), idata);
    Py_XDECREF(cumwt);

/* Return the sorted indices, weights and total weight. */
    return Py_BuildValue("NNf", PyArray_Return(indices), PyArray_Return(sorted), tot);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(wts);
    Py_XDECREF(sorted);
    Py_XDECREF(cumwt);
    Py_XDECREF(indices);
    return NULL;
}


/*------------------------------------------------------------------------------
* Calculate weight ratios for PPS sampling.
------------------------------------------------------------------------------*/
static char prepratios_doc[] =
"prepratios(nsamp, wts, tot):  Prepare weight ratios for PPS sampling." ;

static PyObject *
prepratios (PyObject *self, PyObject *args, PyObject *keywords) {
    int nsamp, npopn;
    PyObject *win;
    PyArrayObject *wts, *ratios, *cumrat;
    PyArrayObject *no_out=NULL;
    npy_intp dims[1];
    int i;
    double tot, *wdata, *rdata;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    wts = NULL;
    ratios = NULL;
    cumrat = NULL;

/* Parse the input weights. */
    Py_TRY(PyArg_ParseTuple(args, "iOd", &nsamp, &win, &tot));
    wts = (PyArrayObject *)
            PyArray_ContiguousFromAny(win, NPY_DOUBLE, 1, 1);
    if (wts == NULL) goto FAIL;
    npopn = wts->dimensions[0];
    /* Note that the nsamp=npopn case can fail; handle this in Python. */
    if (nsamp >= npopn) {
    	PyErr_SetString(PyExc_ValueError, "Requested sample size >= population!");
    	goto FAIL;
    }
    dims[0] = npopn;
    ratios = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (ratios == NULL) goto FAIL;

/* Calculate the weight ratios used in Sampford's PPS algorithm. */
	wdata = DDATA(wts);
	rdata = DDATA(ratios);
	for (i=0; i<npopn; i++) {
		rdata[i] = (wdata[i]/tot) / (1-nsamp*(wdata[i]/tot));
	}

/* Accumulate the ratios. */
	cumrat = (PyArrayObject *) PyArray_CumSum(ratios, 0, NPY_DOUBLE, no_out);

/* Copy results into this instance. */
    set_sample(((PPSSampler*)self)->pSampford, nsamp, rdata, DDATA(cumrat));

/* Return the cumulative ratios array. */
	Py_XDECREF(ratios);
    return PyArray_Return(cumrat);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(wts);
    Py_XDECREF(ratios);
    Py_XDECREF(cumrat);
    return NULL;
}


/*------------------------------------------------------------------------------
* Get a sample without replacement.
------------------------------------------------------------------------------*/
static char samplenr_doc[] = \
  "samplenr(): return an array of sampled indices (PPS without replacement).";
                                                      
static PyObject*
samplenr(PyObject *self, PyObject *args, PyObject *keywords)
{
    PyArrayObject *sample;
    int ntry;
    npy_intp dims[1];
    
    /** Initialize objects to NULL that we may have to XDECREF if a constructor
     ** fails so that XDECREF has an argument to work with.  */
    sample = NULL;

    /* Parse the arguments. */
    Py_TRY(PyArg_ParseTuple(args, ""));

    /* Make an array to hold the sample. */
    dims[0] = ((PPSSampler*)self)->pSampford->nsamp;
    sample = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT);
    if (sample == NULL) goto FAIL;

    ntry = sampford_pps(((PPSSampler*)self)->pSampford, IDATA(sample), &rng_state);

    /* Return the sample. */
    return Py_BuildValue("iN", ntry, PyArray_Return(sample));

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(sample);
    return NULL;
}


/*------------------------------------------------------------------------------
* Calculate weight ratio tables for accelerated PPS sampling.
------------------------------------------------------------------------------*/
static char prepratiotables_doc[] =
"prepratiotables():  Prepare weight ratio tables for PPS sampling.\n"\
"    prepratios(...) should be called first.";

static PyObject *
prepratiotables (PyObject *self, PyObject *args, PyObject *keywords) {
    int npopn;
    double *rdata;

/* No inputs... */
    Py_TRY(PyArg_ParseTuple(args, ""));

/* Build a 5-table sampler for the weight ratios.
   Note that the weight sampler is in *unsorted* order, but the
   Sampford ratios are stored in *sorted* order. */
    npopn = ((PPSSampler*)self)->pSampford->npopn;
    rdata = ((PPSSampler*)self)->pSampford->ratios;
    ((PPSSampler*)self)->pRatioSampler = init_sampler5tbl(rdata, npopn);
    if ( ((PPSSampler*)self)->pRatioSampler == NULL ) goto FAIL;

/* No outputs... */
    Py_INCREF(Py_None);
    return Py_None;

/** Misbehavior ends up here! */
FAIL:
    return NULL;
}


/*------------------------------------------------------------------------------
* Get a sample without replacement, using 5-table samplers.
------------------------------------------------------------------------------*/
static char samplenr5_doc[] = \
"samplenr5(): return an array of sampled indices (PPS without replacement).\n"\
"    5-table lookup is used to accelerate the sampling for large populations.";
                                                      
static PyObject*
samplenr5(PyObject *self, PyObject *args, PyObject *keywords)
{
    PyArrayObject *sample;
    int ntry;
    npy_intp dims[1];
    
    /** Initialize objects to NULL that we may have to XDECREF if a constructor
     ** fails so that XDECREF has an argument to work with.  */
    sample = NULL;

    /* Parse the arguments. */
    Py_TRY(PyArg_ParseTuple(args, ""));

    /* Make an array to hold the sample. */
    dims[0] = ((PPSSampler*)self)->pSampford->nsamp;
    sample = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT);
    if (sample == NULL) goto FAIL;

    ntry = sampford_pps_tables(((PPSSampler*)self)->pSampford, 
               ((PPSSampler*)self)->pSampler, 
               ((PPSSampler*)self)->pRatioSampler, 
               IDATA(sample), &rng_state);

    /* Return the sample. */
    return Py_BuildValue("iN", ntry, PyArray_Return(sample));

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(sample);
    return NULL;
}


/*******************************************************************************
* The new type class methods tables.
*******************************************************************************/

static PyMethodDef PPSSampler_methods[] = {
        {"sample", (PyCFunction)sample, METH_VARARGS | METH_KEYWORDS,
         sample__doc__},
        {"prepwts", (PyCFunction)prepwts, METH_VARARGS | METH_KEYWORDS,
         prepwts_doc},
        {"prepratios", (PyCFunction)prepratios, METH_VARARGS | METH_KEYWORDS,
         prepratios_doc},
        {"samplenr", (PyCFunction)samplenr, METH_VARARGS | METH_KEYWORDS,
         samplenr_doc},
        {"prepratiotables", (PyCFunction)prepratiotables, METH_VARARGS | METH_KEYWORDS,
         prepratiotables_doc},
        {"samplenr5", (PyCFunction)samplenr5, METH_VARARGS | METH_KEYWORDS,
         samplenr5_doc},
        {NULL,          NULL}           /* sentinel */
};


static PyTypeObject PPSSamplerType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_ppssampler._ppssampler",             /*tp_name*/
    sizeof(PPSSampler), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PPSSampler_destroy,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "PPSSampler objects",      /*tp_doc*/
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    PPSSampler_methods,        /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PPSSampler_init,      /* tp_init */
    0,                         /* tp_alloc */
    PPSSampler_new,            /* tp_new */
};


/*******************************************************************************
* Methods for this module.
*******************************************************************************/

/*------------------------------------------------------------------------------
* Set the PRNG state.
------------------------------------------------------------------------------*/
static char set_rng_state_doc[] =
"set_rng_state(key, pos):  Set the RNG state; assumes randomkit's Mersenne twister." ;

static PyObject *
set_rng_state (PyObject *self, PyObject *args) {
    PyObject *keyin;
    PyArrayObject *key;
    int pos;
    char *id;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    key = NULL;

/* Parse the arguments. */
    Py_TRY(PyArg_ParseTuple(args, "sOi", &id, &keyin, &pos));
    key = (PyArrayObject *)
            PyArray_ContiguousFromAny(keyin, PyArray_ULONG, 1, 1);
    if (key == NULL) goto FAIL;
    /** We only support Mersenne twister. */
    if (strcmp(id, "MT19937") != 0) {
        PyErr_SetString(PyExc_ValueError, "Only supports RNG id MT19937!");
            goto FAIL;
    }
    /** Check that the length is right. */
    if (key->dimensions[0] != RK_STATE_LEN) {
        PyErr_SetString(PyExc_ValueError, "Input key length incorrect!");
            goto FAIL;
    }

/* Copy the state data and return. */
    memcpy(rng_state.key, (void*)(key->data), RK_STATE_LEN*sizeof(unsigned long));
    rng_state.pos = pos;
    Py_INCREF(Py_None);
    return Py_None;

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(key);
    return NULL;
}


/*------------------------------------------------------------------------------
* Get the PRNG state.
------------------------------------------------------------------------------*/
static char get_rng_state_doc[] =
"get_rng_state():  Get the RNG state; assumes randomkit's Mersenne twister." ;

static PyObject *
get_rng_state (PyObject *self, PyObject *args) {
    PyArrayObject *key;
    npy_intp dims[1];

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    key = NULL;

/* Parse the arguments. */
    Py_TRY(PyArg_ParseTuple(args, ""));

/* Create the key array. */
    dims[0] = RK_STATE_LEN;
    key = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_ULONG);
    if (key == NULL) goto FAIL;

/* Copy the state data and return. */
    memcpy((void*)(key->data), rng_state.key, RK_STATE_LEN*sizeof(unsigned long));
    return Py_BuildValue("sNi", "MT19937", PyArray_Return(key), rng_state.pos);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(key);
    return NULL;
}


/*------------------------------------------------------------------------------
* Get an integer pseudo random integer in [0,hi].
------------------------------------------------------------------------------*/
static char local_irand_doc[] =
"local_irand(hi):  Get an integer pseudo random number in [0,hi]." ;

static PyObject *
local_irand (PyObject *self, PyObject *args) {
	int hi;

/* Parse the arguments. */
    Py_TRY(PyArg_ParseTuple(args, "i", &hi));

/* Return the value. */
    return Py_BuildValue("l", (long) rk_interval((unsigned long) hi, &rng_state));

/** Misbehavior ends up here! */
FAIL:
    PyErr_SetString(PyExc_RuntimeError, "local_irand parse error!");
    return NULL;
}


/*------------------------------------------------------------------------------
* Get a uniform pseudo random number over [0.,1.).
------------------------------------------------------------------------------*/
static char local_rand_doc[] =
"local_rand():  Get a uniform pseudo random number over [0.,1.)." ;

static PyObject *
local_rand (PyObject *self, PyObject *args) {

/* Parse the arguments. */
    Py_TRY(PyArg_ParseTuple(args, ""));

/* Return the value. */
    return Py_BuildValue("d", rk_double(&rng_state));

/** Misbehavior ends up here! */
FAIL:
    PyErr_SetString(PyExc_RuntimeError, "local_rand parse error!");
    return NULL;
}


/*------------------------------------------------------------------------------
* Sample nsamp indices from integers up to (npopn-1), without replacement.
------------------------------------------------------------------------------*/
static char equalprobi_doc[] =
"equalprobi(ns, np, pool):  Sample ns indices from [0,np-1], without replacement." ;

static PyObject *
equalprobi (PyObject *self, PyObject *args) {
    int nsamp, npopn;
    PyObject *poolin;
    PyArrayObject *sample, *pool;
    int *sdata, *pdata;
    npy_intp dims[1];
    int i, pick;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    sample = NULL;
    pool = NULL;

/* Parse the arguments. */
    Py_TRY(PyArg_ParseTuple(args, "iiO", &nsamp, &npopn, &poolin));
    if (nsamp <= 0) {
    	PyErr_SetString(PyExc_ValueError, "Requested sample size must be > 0!");
    	goto FAIL;
    }
    if (nsamp > npopn) {
    	PyErr_SetString(PyExc_ValueError, "Requested sample larger than population!");
    	goto FAIL;
    }
    /* pool is an int array used as workspace. */
    pool = (PyArrayObject *)
            PyArray_ContiguousFromAny(poolin, PyArray_INT, 1, 1);
    if (pool == NULL) goto FAIL;
    pdata = IDATA(pool);
    /** Check that the dimensions are consistent. */
    if (pool->dimensions[0] < npopn) {
        PyErr_SetString(PyExc_ValueError, "Input pool space too small!");
            goto FAIL;
    }

    dims[0] = nsamp;
    sample = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT);
    if (sample == NULL) goto FAIL;
    sdata = IDATA(sample);

/* Initialize the pool. */
    for (i=0; i<npopn; i++) pdata[i] = i;

/* Pick from the pool at random; shuffle the last item down to replace the 
   one picked. */
    for (i=0; i<nsamp; i++) {
   		pick = npopn * (1.*rand())/(RAND_MAX+1.);  /* Always < npopn. */
/*    		PySys_WriteStdout("pick, npopn: %i  %i  %i\n", pick, pdata[pick], npopn); */
   		sdata[i] = pdata[pick];
   		npopn--;  /* Now points to the last item (at npopn-1). */
   		pdata[pick] = pdata[npopn]; /* Note if pick=npopn, it just rolls off the end. */
    }

/* Return the sample. */
    return PyArray_Return(sample);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(pool);
    Py_XDECREF(sample);
    return NULL;
}


/*------------------------------------------------------------------------------
* Sample nsamp integers from an input sequence, without replacement.
------------------------------------------------------------------------------*/
static char equalprob_doc[] =
"equalprob(ns, pool):  Sample ns integers from the input pool, without replacement." ;

static PyObject *
equalprob (PyObject *self, PyObject *args) {
    int nsamp, npopn;
    PyObject *poolin;
    PyArrayObject *sample, *pool;
    int *sdata, *pdata;
    npy_intp dims[1];
    int i, pick;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    sample = NULL;
    pool = NULL;

/* Parse the arguments. */
    Py_TRY(PyArg_ParseTuple(args, "iO", &nsamp, &poolin));
    if (nsamp <= 0) {
    	PyErr_SetString(PyExc_ValueError, "Requested sample size must be > 0!");
    	goto FAIL;
    }
    /* pool contains the collection of integers.  Use a copy since we'll
       shuffle it. */
    pool = (PyArrayObject *)
            PyArray_CopyFromObject(poolin, PyArray_INT, 1, 1);
    if (pool == NULL) goto FAIL;
    pdata = IDATA(pool);
    npopn = pool->dimensions[0];
    /* Check population size. */
    if (nsamp > npopn) {
    	PyErr_SetString(PyExc_ValueError, "Requested sample larger than population!");
    	goto FAIL;
    }

    dims[0] = nsamp;
    sample = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT);
    if (sample == NULL) goto FAIL;
    sdata = IDATA(sample);

/* Pick from the pool at random; shuffle the last item down to replace the 
   one picked. */
    for (i=0; i<nsamp; i++) {
   		pick = npopn * (1.*rand())/(RAND_MAX+1.);  /* Always < npopn. */
/*    		PySys_WriteStdout("pick, npopn: %i  %i  %i\n", pick, pdata[pick], npopn); */
   		sdata[i] = pdata[pick];
   		npopn--;  /* Now points to the last item (at npopn-1). */
   		pdata[pick] = pdata[npopn]; /* Note if pick=npopn, it just rolls off the end. */
    }

/* Return the sample. */
    return PyArray_Return(sample);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(pool);
    Py_XDECREF(sample);
    return NULL;
}


/*******************************************************************************
* The methods table defining how Python accesses this module's methods.  This
* must have an element for each new method and end with the "sentinel" element.
* It should have 4 entries for each method (except the sentinel): 
*         {call name, function to call, argument type, doc string},
* with each entry ending with a comma.
*******************************************************************************/

/* Module functions */
static PyMethodDef module_methods[] = {
    {"set_rng_state", set_rng_state, METH_VARARGS, set_rng_state_doc},
    {"get_rng_state", get_rng_state, METH_VARARGS, get_rng_state_doc},
    {"local_irand", local_irand, METH_VARARGS, local_irand_doc},
    {"local_rand", local_rand, METH_VARARGS, local_rand_doc},
    {"equalprobi", equalprobi, METH_VARARGS, equalprobi_doc},
    {"equalprob", equalprob, METH_VARARGS, equalprob_doc},
    {NULL}  /* Sentinel */
};


/*******************************************************************************
* The module doc string and initialization function.
*
* The init function should be the only non-static global object in this file.
*******************************************************************************/

/* Doc strings: */
static char ppssampler__doc__[] = \
  "A module allowing fast sampling from a finite population with probabilities\n"\
  "proportional to weights.\n"\
  "\n"\
  "Use the syntax:\n"\
  ">>> s = _ppssampler(table)\n"\
  "to create an object for sampling from a distribution with probability\n"\
  "mass function given by 'table', where:\n"\
  "\n"\
  "x       0           1           ...     k-1\n"\
  "p(x)    table[0]    table[1]            table[k-1]\n"\
  "\n"\
  "The values of table[i] need not be normalized to sum to 1, but must be\n"\
  "non-negative.\n";


#ifndef PyMODINIT_FUNC  /* declarations for shared library import/export */
#define PyMODINIT_FUNC DL_EXPORT(void)
#endif
PyMODINIT_FUNC
init_ppssampler(void)
{
    PyObject *m, *d, *s;
    
    /* Initialize scipy */
    import_array();   
    /* Import the ufunc objects */
    // import_ufunc();

    // PPSSamplerType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PPSSamplerType) < 0)
        return;

    /* Create the module and add the functions */
    m = Py_InitModule3("_ppssampler", module_methods, ppssampler__doc__); 
    
    if (m == NULL)
        return;
    
    // /* Add our class */
    // PyStructSequence_InitType(&SamplerType, &sampler_type_desc);
    Py_INCREF(&PPSSamplerType);
    PyModule_AddObject(m, "_ppssampler", (PyObject*) &PPSSamplerType);


    /* Add some symbolic constants to the module */
    d = PyModule_GetDict(m);

    s = PyString_FromString("0.1-alpha");
    PyDict_SetItemString(d, "__version__", s);
    Py_DECREF(s);

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module _ppssampler");
}






