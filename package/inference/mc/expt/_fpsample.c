/*******************************************************************************
* _fpsample
*
* Finite population samplers for Python.
*
* Created:  18 May 2006 by Tom Loredo
* Modified: 19 May 2006
*******************************************************************************/

/*******************************************************************************
* Headers, prototypes, and macros.
*******************************************************************************/
#include "Python.h"

/* Other includes needed for this extension. */
#include "numpy/arrayobject.h"  /* Header needed to use scipy array types. */
#include <stdlib.h>
#include "randomkit.h"


/* Error handler macro for use in any extension; adapted from Dubois & Yang. 
   This requires that the function using it have the FAIL label defined. */
#define Py_TRY(BOOLEAN) {if (!(BOOLEAN)) goto FAIL;}

/* Macros to cast numpy array data to a 1-d arrays; borrowed
   from lapack_litemodule.c in numpy. */
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

/*------------------------------------------------------------------------------
*  A handy function from lapack_litemodule that verifies that an array
*  argument passed in from Python has the right type and is contiguous. 
*
*  SciPy_CheckObject(object, type, name, typename, funcname) 
*		object	-	A Python object
*		type	-	A numpy type as defined in array.h (e.g., PyArray_DOUBLE)
*		name 	-	Diagnostic string for printout; the name of the object
*		typename-	Diagnostic string for printout; the name of the type expected
*		funcname-	Diagnostic string for printout; the name of the function
*
*  The function returns 0 if there is a problem; it will have stored a
*  Python string object in the global variable ErrorObject describing the error.
*  This string object can be accessed if necessary via SciPy_CheckObject_Error();
*  use this rather than directly using ErrorObject to keep the proper refcount.
------------------------------------------------------------------------------*/
static PyObject *ErrorObject;

/* This isn't used here; see comment above for use case. */ 
static PyObject *SciPy_CheckObject_Error(void)
{  if (! ErrorObject)
      ErrorObject = PyString_FromString("SciPy_CheckObject Error");
   Py_INCREF(ErrorObject);
   return ErrorObject;
}

static PyObject *ErrorReturn(char *mes)
{  if (!ErrorObject)
      ErrorObject = PyString_FromString("SciPy_CheckObject error");
   PyErr_SetString(ErrorObject,mes);
   return NULL;
}

static int SciPy_CheckObject(PyObject *ob, int t, char *obname,
        char *tname, char *funname) {

        char buf[255];

        if (! PyArray_Check(ob))
        {       sprintf(buf,"Expected an array for parameter %s in %s",obname, funname);
        ErrorReturn(buf);
        return 0;
    }
    if (!(((PyArrayObject *)ob)->flags & CONTIGUOUS)) 
    {    sprintf(buf,"Parameter %s is not contiguous in %s",obname, funname);
        ErrorReturn(buf);
        return 0;
    }
    if (!(((PyArrayObject *)ob)->descr->type_num == t))
    {    sprintf(buf,"Parameter %s is not of type %s in %s",obname,tname, funname);
        ErrorReturn(buf);
        return 0;
    }
    return 1;
}

/*******************************************************************************
* Methods for this module.
*
* The functions that will actually be called by Python to implement the methods
* are defined here, as well as any global variables they may use to communicate
* with each other and any other "helper" functions not defined elsewhere.  The
* functions that will actually be called by Python must return PyObject 
* pointers.
*
* Be sure to define a documentation string for each method.
*
* The functions and global variables here should all be static (i.e., accessible
* only from within this module) to prevent side effects from name collisions
* outside this extension.  They can have any names (the names Python will
* use to access them are defined separately, below); they here follow a 
* convention that should be self-explanatory.  
*******************************************************************************/

/*------------------------------------------------------------------------------
* Local storage of the state of the PRNG, to be copied from numpy
* to ensure reproducibility when the user uses random.set_state to
* repeat runs.
------------------------------------------------------------------------------*/
static rk_state rng_state;

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

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    key = NULL;

/* Parse the arguments. */
    Py_TRY(PyArg_ParseTuple(args, "Oi", &keyin, &pos));
    key = (PyArrayObject *)
            PyArray_ContiguousFromObject(keyin, PyArray_LONG, 1, 1);
    if (key == NULL) goto FAIL;
    /** Check that the length is right. */
    if (key->dimensions[0] != RK_STATE_LEN) {
        PyErr_SetString(PyExc_ValueError, "Input key length incorrect!");
            goto FAIL;
    }

/* Copy the state data and return. */
    memcpy(rng_state.key, (void*)(key->data), RK_STATE_LEN*sizeof(long));
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
    int dims[1];

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    key = NULL;

/* Parse the arguments. */
    Py_TRY(PyArg_ParseTuple(args, ""));

/* Create the key array. */
    dims[0] = RK_STATE_LEN;
    key = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_LONG);
    if (key == NULL) goto FAIL;

/* Copy the state data and return. */
    memcpy((void*)(key->data), rng_state.key, RK_STATE_LEN*sizeof(long));
    return Py_BuildValue("Ni", PyArray_Return(key), rng_state.pos);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(key);
    return NULL;
}


/*------------------------------------------------------------------------------
* Get an integer pseudo random number in [0,hi].
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
* Return a 5-table sampler state as a tuple.
------------------------------------------------------------------------------*/
static char get_5tbl_state_doc[] =
"" ;

static PyObject *
get_5tbl_state (Sampler* sampler) {
    PyArrayObject *AA, *BB, *CC, *DD, *EE;
    int dims[0];

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    AA = NULL;
    BB = NULL;
    CC = NULL;
    DD = NULL;
    EE = NULL;

	dims[0] = sampler->sizeAA;
    AA = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_LONG);
    if (AA == NULL) goto FAIL;
	dims[0] = sampler->sizeBB;
    BB = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_LONG);
    if (BB == NULL) goto FAIL;
	dims[0] = sampler->sizeCC;
    CC = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_LONG);
    if (CC == NULL) goto FAIL;
	dims[0] = sampler->sizeCC;
    DD = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_LONG);
    if (DD == NULL) goto FAIL;
	dims[0] = sampler->sizeEE;
    EE = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_LONG);
    if (EE == NULL) goto FAIL;

    return Py_BuildValue("NNNNNl", 
               Py_ArrayReturn(AA), Py_ArrayReturn(BB), Py_ArrayReturn(CC), 
               Py_ArrayReturn(DD), Py_ArrayReturn(EE), 
               sampler->prob1event);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(AA);
    Py_XDECREF(BB);
    Py_XDECREF(CC);
    Py_XDECREF(DD);
    Py_XDECREF(EE);
    return NULL;

}


/*------------------------------------------------------------------------------
* Set a 5-table sampler state from a tuple.
------------------------------------------------------------------------------*/
static char set_5tbl_state_doc[] =
"" ;

void
set_5tbl_state (PyObject* args, Sampler* sampler) {
	PyObject *Ain, *Bin, *Cin, *Din, *Ein;
    PyArrayObject *AA, *BB, *CC, *DD, *EE;
    int dims[0];
    long prob1event;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    AA = NULL;
    BB = NULL;
    CC = NULL;
    DD = NULL;
    EE = NULL;

/* Parse the arguments. */
    Py_TRY(PyArg_ParseTuple(args, "OOOOOl", &Ain, &Bin, &Cin, &Din, &Ein, &prob1event));

    AA = (PyArrayObject *)PyArray_ContiguousFromObject(Ain, PyArray_LONG, 1, 1);
    if (AA == NULL) goto FAIL;
    sampler->sizeAA = (int32_t) AA->dimensions[0];

    BB = (PyArrayObject *)PyArray_ContiguousFromObject(Bin, PyArray_LONG, 1, 1);
    if (BB == NULL) goto FAIL;
    sampler->sizeBB = (int32_t) BB->dimensions[0];

    CC = (PyArrayObject *)PyArray_ContiguousFromObject(Cin, PyArray_LONG, 1, 1);
    if (CC == NULL) goto FAIL;
    sampler->sizeCC = (int32_t) CC->dimensions[0];

    DD = (PyArrayObject *)PyArray_ContiguousFromObject(Din, PyArray_LONG, 1, 1);
    if (DD == NULL) goto FAIL;
    sampler->sizeDD = (int32_t) DD->dimensions[0];

    EE = (PyArrayObject *)PyArray_ContiguousFromObject(Ein, PyArray_LONG, 1, 1);
    if (EE == NULL) goto FAIL;
    sampler->sizeEE = (int32_t) EE->dimensions[0];

	sampler->prob1event = prob1event;

    sampler->t1 = sampler->sizeAA<<24;
    sampler->t2 = sampler->t1+(sampler->sizeBB<<18);
    sampler->t3 = sampler->t2+(sampler->sizeCC<<12);
    sampler->t4 = sampler->t3+(sampler->sizeDD<<6);

    return;

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(AA);
    Py_XDECREF(BB);
    Py_XDECREF(CC);
    Py_XDECREF(DD);
    Py_XDECREF(EE);
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
    int dims[1];
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
            PyArray_ContiguousFromObject(poolin, PyArray_INT, 1, 1);
    if (pool == NULL) goto FAIL;
    pdata = IDATA(pool);
    /** Check that the dimensions are consistent. */
    if (pool->dimensions[0] < npopn) {
        PyErr_SetString(PyExc_ValueError, "Input pool space too small!");
            goto FAIL;
    }

    dims[0] = nsamp;
    sample = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_INT);
    if (sample == NULL) goto FAIL;
    sdata = IDATA(sample);

/* Initialize the pool. */
    for (i=0; i<npopn; i++) pdata[i] = i;

/* Pick from the pool at random; shuffle the last item down to replace the 
   one picked. */
    for (i=0; i<nsamp; i++) {
   		pick = npopn * (1.*rand())/(RAND_MAX+1.);  /* Always < npopn. */
   		PySys_WriteStdout("pick, npopn: %i  %i  %i\n", pick, pdata[pick], npopn);
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
    int dims[1];
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
    sample = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_INT);
    if (sample == NULL) goto FAIL;
    sdata = IDATA(sample);

/* Pick from the pool at random; shuffle the last item down to replace the 
   one picked. */
    for (i=0; i<nsamp; i++) {
   		pick = npopn * (1.*rand())/(RAND_MAX+1.);  /* Always < npopn. */
   		PySys_WriteStdout("pick, npopn: %i  %i  %i\n", pick, pdata[pick], npopn);
   		sdata[i] = pdata[pick];
   		npopn--;  /* Now points to the last item (at npopn-1). */
   		pdata[pick] = pdata[npopn]; /* Note if pick=npopn, 
   									   it just rolls off the end. */
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
* Prepare weights for PPS sampling (sort, accumulate).
------------------------------------------------------------------------------*/
static char prepwts_doc[] =
"prepwts(wts):  Prepare weights for PPS sampling (sort, accumulate)." ;

static PyObject *
prepwts (PyObject *self, PyObject *args) {
    int npopn;
    PyObject *win;
    PyArrayObject *wts, *sorted, *cumwt, *indices;
    int dims[1], *idata;
    int i, j, itemp;
    double *wdata, *sdata;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    wts = NULL;
    sorted = NULL;
    cumwt = NULL;
    indices = NULL;

/* Parse the input weights. */
    Py_TRY(PyArg_ParseTuple(args, "O", &win));
    wts = (PyArrayObject *)
            PyArray_ContiguousFromObject(win, PyArray_DOUBLE, 1, 1);
    if (wts == NULL) goto FAIL;
    npopn = wts->dimensions[0];
    dims[0] = npopn;
    sorted = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
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

/* Accumulate the reverse-sorted weights. */
	cumwt = (PyArrayObject *) PyArray_CumSum(sorted, 0, PyArray_DOUBLE);

/* Return the arrays. */
    return Py_BuildValue("NNN", PyArray_Return(indices), PyArray_Return(sorted), 
    	PyArray_Return(cumwt));

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
static char wtratios_doc[] =
"wtratios(nsamp, wts, tot):  Prepare weight ratios for PPS sampling." ;

static PyObject *
wtratios (PyObject *self, PyObject *args) {
    int nsamp, npopn;
    PyObject *win;
    PyArrayObject *wts, *ratios, *cumrat;
    int dims[1];
    int i;
    double tot, *wdata, *rdata;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    wts = NULL;
    ratios = NULL;
    cumrat = NULL;

/* Parse the input weights. */
    Py_TRY(PyArg_ParseTuple(args, "iOf", &nsamp, &win, &tot));
    wts = (PyArrayObject *)
            PyArray_ContiguousFromObject(win, PyArray_DOUBLE, 1, 1);
    if (wts == NULL) goto FAIL;
    npopn = wts->dimensions[0];
    if (nsamp > npopn) {
    	PyErr_SetString(PyExc_ValueError, "Requested sample larger than population!");
    	goto FAIL;
    }
    dims[0] = npopn;
    ratios = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
    if (ratios == NULL) goto FAIL;

/* Calculate the weight ratios used in Sampford's PPS algorithm. */
	wdata = DDATA(wts);
	rdata = DDATA(ratios);
	for (i=0; i<npopn; i++) {
		rdata[i] = (wdata[i]/tot) / (1-nsamp*(wdata[i]/tot));
	}

/* Accumulate the ratios. */
	cumrat = (PyArrayObject *) PyArray_CumSum(ratios, 0, PyArray_DOUBLE);

/* Return the arrays. */
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
* Generate a PPS sample via Sampford's rejection algorithm.
------------------------------------------------------------------------------*/
static char SampfordPPS_doc[] =
"SampfordPPS(nsamp, cumwt, cumrat):  Generate a PPS sample via Sampford's \n"\
"rejection algorithm." ;

static PyObject *
SampfordPPS (PyObject *self, PyObject *args) {
    int nsamp, npopn;
    PyObject *win, *rin;
    PyArrayObject *cumwt, *cumrat, *sample;
    int dims[1];
    int i, j, k, duplicate, ntry;
    int *sdata;
    double tot, val, *wdata, *rdata;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    cumwt = NULL;
    cumrat = NULL;
    sample = NULL;

/* Parse the inputs. */
    Py_TRY(PyArg_ParseTuple(args, "iOO", &nsamp, &win, &rin));
    cumwt = (PyArrayObject *)
            PyArray_ContiguousFromObject(win, PyArray_DOUBLE, 1, 1);
    if (cumwt == NULL) goto FAIL;
    npopn = cumwt->dimensions[0];
    cumrat = (PyArrayObject *)
            PyArray_ContiguousFromObject(rin, PyArray_DOUBLE, 1, 1);
    if (cumrat == NULL) goto FAIL;
    if (npopn != cumrat->dimensions[0]) {
    	PyErr_SetString(PyExc_ValueError, "Mismatched weight and ratio arrays!");
    	goto FAIL;
    }
    if (nsamp < 1) {
    	PyErr_SetString(PyExc_ValueError, "Sample size must be >= 1!");
    	goto FAIL;
    }
    if (nsamp > npopn) {
    	PyErr_SetString(PyExc_ValueError, "Requested sample larger than population!");
    	goto FAIL;
    }
    dims[0] = nsamp;
    sample = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_INT);
    if (sample == NULL) goto FAIL;
    wdata = DDATA(cumwt);
    rdata = DDATA(cumrat);
    sdata = IDATA(sample);

/* Cycle the rejection loop until nsamp distinct members get selected. */
	ntry = 0;
	do {
		ntry++;

/* Pick the first value based on the raw weights. */
		val = wdata[npopn-1] * (1.*rand())/(RAND_MAX+1.);
		for (i=0; i<npopn; i++) {
			if (val <= wdata[i]) break;
		}
		sdata[0] = i;
		/* PySys_WriteStdout("First sample unit:  %i\n", i); */
	
/* Pick subsequent samples based on the ratios, rejecting sequences with duplicates. */
		duplicate = 0;
		tot = rdata[npopn-1];
		for (j=1; j<nsamp; j++) {
			val = tot * (1.*rand())/(RAND_MAX+1.);
			for (i=0; i<npopn; i++) {
				if (val <= rdata[i]) break;
			}
			/* If i was already selected, note a duplicate and end the ratio loop. */
			for (k=0; k<j; k++) {
				if (sdata[k] == i) {
					duplicate = 1;
					break;
				}
			}
			/* PySys_WriteStdout("j, choice, dup:  %i  %i  %i\n", j, i, duplicate); */
			if (duplicate) break;
			sdata[j] = i;
		}

/* End the rejection loop if there are no duplicates; else start over. */
	} while (duplicate);

/* Return the sample. */
    return Py_BuildValue("iN", ntry, PyArray_Return(sample));

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(cumwt);
    Py_XDECREF(cumrat);
    Py_XDECREF(sample);
    return NULL;
}


/*------------------------------------------------------------------------------
* Generate a PPS sample via Sampford's rejection algorithm, using 5-table
* lookup for the rejection loop part of the algorithm.
------------------------------------------------------------------------------*/
static char SampfordPPS5_doc[] =
"SampfordPPS(nsamp, cumwt, cumrat):  Generate a PPS sample via Sampford's \n"\
"rejection algorithm, using 5-table lookup to accelerate the rejection loop." ;

static PyObject *
SampfordPPS (PyObject *self, PyObject *args) {
    int nsamp, npopn;
    PyObject *win, *rin;
    PyArrayObject *cumwt, *cumrat, *sample;
    int dims[1];
    int i, j, k, duplicate, ntry;
    int *sdata;
    double tot, val, *wdata, *rdata;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  */
    cumwt = NULL;
    cumrat = NULL;
    sample = NULL;

/* Parse the inputs. */
    Py_TRY(PyArg_ParseTuple(args, "iOO", &nsamp, &win, &rin));
    cumwt = (PyArrayObject *)
            PyArray_ContiguousFromObject(win, PyArray_DOUBLE, 1, 1);
    if (cumwt == NULL) goto FAIL;
    npopn = cumwt->dimensions[0];
    cumrat = (PyArrayObject *)
            PyArray_ContiguousFromObject(rin, PyArray_DOUBLE, 1, 1);
    if (cumrat == NULL) goto FAIL;
    if (npopn != cumrat->dimensions[0]) {
    	PyErr_SetString(PyExc_ValueError, "Mismatched weight and ratio arrays!");
    	goto FAIL;
    }
    if (nsamp < 1) {
    	PyErr_SetString(PyExc_ValueError, "Sample size must be >= 1!");
    	goto FAIL;
    }
    if (nsamp > npopn) {
    	PyErr_SetString(PyExc_ValueError, "Requested sample larger than population!");
    	goto FAIL;
    }
    dims[0] = nsamp;
    sample = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_INT);
    if (sample == NULL) goto FAIL;
    wdata = DDATA(cumwt);
    rdata = DDATA(cumrat);
    sdata = IDATA(sample);

/* Cycle the rejection loop until nsamp distinct members get selected. */
	ntry = 0;
	do {
		ntry++;

/* TODO:  For large npopn, use Marsaglia et al.'s 5-table algorithm for index slxn. */

/* Pick the first value based on the raw weights. */
		val = wdata[npopn-1] * (1.*rand())/(RAND_MAX+1.);
		for (i=0; i<npopn; i++) {
			if (val <= wdata[i]) break;
		}
		sdata[0] = i;
		/* PySys_WriteStdout("First sample unit:  %i\n", i); */
	
/* Pick subsequent samples based on the ratios, rejecting sequences with duplicates. */
		duplicate = 0;
		tot = rdata[npopn-1];
		for (j=1; j<nsamp; j++) {
			val = tot * (1.*rand())/(RAND_MAX+1.);
			for (i=0; i<npopn; i++) {
				if (val <= rdata[i]) break;
			}
			/* If i was already selected, note a duplicate and end the ratio loop. */
			for (k=0; k<j; k++) {
				if (sdata[k] == i) {
					duplicate = 1;
					break;
				}
			}
			/* PySys_WriteStdout("j, choice, dup:  %i  %i  %i\n", j, i, duplicate); */
			if (duplicate) break;
			sdata[j] = i;
		}

/* End the rejection loop if there are no duplicates; else start over. */
	} while (duplicate);

/* Return the sample. */
    return Py_BuildValue("iN", ntry, PyArray_Return(sample));

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(cumwt);
    Py_XDECREF(cumrat);
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
static PyMethodDef methods[] = {
    {"equalprobi", equalprobi, METH_VARARGS, equalprobi_doc},
    {"equalprob", equalprob, METH_VARARGS, equalprob_doc},
    {"prepwts", prepwts, METH_VARARGS, prepwts_doc},
    {"wtratios", wtratios, METH_VARARGS, wtratios_doc},
    {"SampfordPPS", SampfordPPS, METH_VARARGS, SampfordPPS_doc},
    {"SampfordPPS5", SampfordPPS, METH_VARARGS, SampfordPPS_doc},
    {"set_rng_state", set_rng_state, METH_VARARGS, set_rng_state_doc},
    {"get_rng_state", get_rng_state, METH_VARARGS, get_rng_state_doc},
    {"local_irand", local_irand, METH_VARARGS, local_irand_doc},
    {"local_rand", local_rand, METH_VARARGS, local_rand_doc},
    {NULL,		NULL, 0}		/* sentinel marking end of methods list */
};


/*******************************************************************************
* The initialization function---it must be named for the module ("init" followed
* by the module name).  Change its name, and the two occurrences of "example"
* in the function definition, to match the module name.
*
* Note the import_array() call!!!!  This is only needed if you use scipy
* arrays in an extension; non-array extensions should omit this.  But if you use 
* scipy arrays and forget this, you'll likely crash when you try to use your
* extension.
*
* This should be the only non-static global object in this file.
*******************************************************************************/
void init_fpsample(void) {
  PyObject *m;
  
  /* Create the module and add the functions */
  m = Py_InitModule("_fpsample", methods);
  import_array();

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module _fpsample");
}

