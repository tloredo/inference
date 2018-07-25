
/*******************************************************************************
* vba: "Vector Bretthorst Algorithm"
*
* A C extension module for Python.
*
* The core of the Bretthorst algorithm: calculating the metric and
* finding its Jacobian and the sufficient statistic.
*
* Created:  12 Apr 2000  Tom Loredo
* Modified: 11 May 2000  TL - Fixed memory leaks
*           12 May 2000  Speedup
*******************************************************************************/

/*******************************************************************************
* Headers, prototypes, and macros.
*******************************************************************************/

/* Include Python headers; needed for *all* extensions. */
#include "Python.h"

/* Other includes needed for this extension. */
#include "numpy/arrayobject.h"
#include <math.h>

/* Prototypes for this extension that aren't in an include file. */
/* None are needed for this example. */

/* Error handler macro for use in any extension; adapted from Dubois & Yang. */
#define Py_TRY(BOOLEAN) {if (!(BOOLEAN)) goto FAIL;}

/* Macro to cast NumPy array data to a 1-d double array. */
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))

/*------------------------------------------------------------------------------
* Some handy error handling functions, from lapack_litemodule.c in old Numeric.
------------------------------------------------------------------------------*/
static PyObject *ErrorObject;

/* Not currently used.  See ExampleScipyModule for use case. */
static PyObject *vbaError(void)
{  if (! ErrorObject)
      ErrorObject = PyString_FromString("vbaError");
   Py_INCREF(ErrorObject);
   return ErrorObject;
}

static PyObject *ErrorReturn(char *mes)
{  if (!ErrorObject)
      ErrorObject = PyString_FromString("mkarrayError");
   PyErr_SetString(ErrorObject,mes);
   return NULL;
}

#define TRY(E) if( ! (E)) return NULL

/* Not currently used.  See ExampleScipyModule for use case. */
static int SciPy_CheckObject(PyObject *ob, int t, char *obname,
        char *tname, char *funname)
{       char buf[255];
        if (! PyArray_Check(ob))
        {       sprintf(buf,"Expected an array for parameter %s in lapack_dge.%s",obname, funname);
        ErrorReturn(buf);
        return 0;
    }
    if (!(((PyArrayObject *)ob)->flags & NPY_CONTIGUOUS)) 
    {    sprintf(buf,"Parameter %s is not contiguous in lapack_dge.%s",obname, funname);
        ErrorReturn(buf);
        return 0;
    }
    if (!(((PyArrayObject *)ob)->descr->type_num == t))
    {    sprintf(buf,"Parameter %s is not of type %s in lapack_lite.%s",obname,tname, funname);
        ErrorReturn(buf);
        return 0;
    }
    return 1;
}



/*******************************************************************************
* Methods for this module.

* The functions that will actually be called by python to implement the methods
* are defined here.
* Be sure to define a documentation string for each method.
* These static objects can have any names; they here follow a convention
* that should be self-explanatory. 
*******************************************************************************/

/*------------------------------------------------------------------------------
* Calculate the model metric & related quantities needed for the marginal
* distribution for the nonlinear parameters.  Return them in a tuple:
*
* (metric, L, rdet, P, A, S)
*
*    metric = model metric (MxM matrix)
*    L      = Cholesky factorization of metric (lower triangular MxM)
*    rdet   = Sqrt of determinant of the metric
*    P      = Projections of data on models (M vector)
*    A      = Amplitude estimates (M vector)
*    S      = Sufficient statistic
*
* The Cholesky algorithms are adapted from *Numerical Recipes*.
------------------------------------------------------------------------------*/
static char vba_metricAnalysis_doc[] =
"metricAnalysis(dmat,data):  Calculate the metric from the design matrix; get various derived quantities from it & the data via Cholesky.";

static PyObject *
vba_metricAnalysis (PyObject *self, PyObject *args) {
    int mdims[2], pdims[1];
    int nmod, nsamp, a, b, i, m, an, bn;
    PyObject *input1, *input2;
    PyArrayObject *dmat, *smpls;
    PyArrayObject *metric, *cfmetric, *proj, *amps;
    double *dmdata, *mdata, *cfdata, *diag, *pdata, *sdata, *adata;
    double sum, rdet;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  Strictly, this
 ** isn't needed for every object. */
    dmat = NULL;
    smpls = NULL;
    metric = NULL;
    proj = NULL;
    cfmetric = NULL;
    amps = NULL;
    diag = NULL;  /* This is for malloc/free, not Python. */

/** Parse the design matrix & data; make them contiguous. */
    Py_TRY(PyArg_ParseTuple(args, "OO", &input1, &input2));
    dmat = (PyArrayObject *)
            PyArray_ContiguousFromAny(input1, PyArray_DOUBLE, 2, 2);
    if (dmat == NULL) goto FAIL;
    smpls = (PyArrayObject *)
            PyArray_ContiguousFromAny(input2, PyArray_DOUBLE, 1, 1);
    if (smpls == NULL) goto FAIL;

/** Check that their dimensions are consistent. */
    if (dmat->dimensions[1] != smpls->dimensions[0]) {
        PyErr_SetString(PyExc_ValueError, "inconsistent data dimensions!");
            goto FAIL;
    }

/** Store the # of models & samples; make vector references to dmat & data. */
    nmod = dmat->dimensions[0];
    nsamp = dmat->dimensions[1];
    dmdata = DDATA(dmat);
    sdata = DDATA(smpls);

/** Calculate the metric. */
    mdims[0] = nmod;
    mdims[1] = nmod;
    metric = (PyArrayObject *)PyArray_FromDims(2, mdims, PyArray_DOUBLE);
    if (metric == NULL) goto FAIL;
    mdata = DDATA(metric);
    for (a=0; a<nmod; a++) {
        an = a*nsamp;
        for (b=0; b<nmod; b++) {
            bn = b*nsamp;
            sum = 0.;
            for (i=0; i<nsamp; i++) {
                sum += *(dmdata+an+i) * *(dmdata+bn+i);
            }
            *(mdata+a*nmod+b) = sum;
        }
    }

/** Calculate the projections of the data on the models. */
    pdims[0] = nmod;
    proj = (PyArrayObject *)PyArray_FromDims(1, pdims, PyArray_DOUBLE);
    if (proj == NULL) goto FAIL;
    pdata = DDATA(proj);
    for (a=0; a<nmod; a++) {
        an = a*nsamp;
        sum = 0.;
        for (i=0; i<nsamp; i++) {
            sum += *(dmdata+an+i) * *(sdata+i);
        }
        *(pdata+a) = sum;
    }

/** Perform a Cholesky factorization of the metric. Calculate
*** the sqrt(determinant) along the way.  */
    cfmetric = (PyArrayObject *)PyArray_Copy(metric);
    if (cfmetric == NULL) goto FAIL;
    cfdata = DDATA(cfmetric);
    diag = (double *) PyMem_Malloc(nmod*sizeof(double));
    rdet = 1.;
    for (a=0; a<nmod; a++) {
        an = a*nmod;
        for (b=a; b<nmod; b++) {
            bn = b*nmod;
            sum = *(cfdata+a*nmod+b);
            for (m=a-1; m>=0; m--) {
                sum -= *(cfdata+an+m) * *(cfdata+bn+m);
            }
            if (a == b) {
                if (sum <= 0.) {
                    PyErr_SetString(PyExc_ValueError,
                        "metric is not positive definite!");
                    goto FAIL;
                }
                *(diag+a) = sqrt(sum);
                rdet = rdet * *(diag+a);
            } else *(cfdata+bn+a) = sum/ *(diag+a);
        }
    }

/** We have the factorization in the lower triangle and diag;
*** replace the diagonal and zero the upper triangle. */
    for (a=0; a<nmod; a++) {
        an = a*nmod;
        for (b=a; b<nmod; b++) {
            if (a == b) {
                *(cfdata+an+a) = *(diag+a);
            } else *(cfdata+an+b) = 0.;
        }
    }

/** Calculate the amplitude estimates by backsubstitution. */
    pdims[0] = nmod;
    amps = (PyArrayObject *)PyArray_FromDims(1, pdims, PyArray_DOUBLE);
    if (amps == NULL) goto FAIL;
    adata = DDATA(amps);
    for (a=0; a<nmod; a++) {
        an = a*nmod;
        sum = *(pdata+a);
        for (b=a-1; b>=0; b--) {
            sum -= *(cfdata+an+b) * *(adata+b);
        }
        *(adata+a) = sum / *(cfdata+an+a);
    }
    for (a=nmod-1; a>=0; a--) {
        sum = *(adata+a);
        for (b=a+1; b<nmod; b++) {
            sum -= *(cfdata+b*nmod+a) * *(adata+b);
        }
        *(adata+a) = sum / *(cfdata+a*nmod+a);
    }

/** Calculate the sufficient statistic. */
    sum = 0.;
    for (a=0; a<nmod; a++) {
        sum += *(adata+a) * *(pdata+a);
    }

/** Clean up and return everything in a tuple. */
    PyMem_Free(diag);
    Py_DECREF(dmat);
    Py_DECREF(smpls);
    return Py_BuildValue("NNfNNf",PyArray_Return(metric),
        PyArray_Return(cfmetric),rdet,
        PyArray_Return(proj),
        PyArray_Return(amps),sum);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(dmat);
    Py_XDECREF(smpls);
    Py_XDECREF(metric);
    Py_XDECREF(proj);
    Py_XDECREF(cfmetric);
    Py_XDECREF(amps);
    if (diag != NULL) PyMem_Free(diag);
    return NULL;
}

/*------------------------------------------------------------------------------
* Calculate covariance matrix for the amplitudes from the Cholesky
* factorization of the metric.
------------------------------------------------------------------------------*/
static char vba_covar_doc[] =
"covar(L):  Calculate the covariance matrix for the amplitudes from the Cholesky factorization of the metric.";

static PyObject *
vba_covar (PyObject *self, PyObject *args) {
    int nmod, a, b, m, an, bn;
    PyObject *input;
    PyArrayObject *cfmetric, *covar;
    double *cvdata;
    double sum;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  Strictly, this
 ** isn't needed for every object. */
    cfmetric = NULL;
    covar = NULL;

/** Parse the Cholesky factor & make it contiguous. */
    Py_TRY(PyArg_ParseTuple(args, "O", &input));
    cfmetric = (PyArrayObject *)
            PyArray_ContiguousFromAny(input, PyArray_DOUBLE, 2, 2);
    if (cfmetric == NULL) goto FAIL;

/** Copy the matrix and gather needed info. */
    covar = (PyArrayObject *)PyArray_Copy(cfmetric);
    if (covar == NULL) goto FAIL;
    nmod = cfmetric->dimensions[0];
    cvdata = DDATA(covar);

/** Find the inverse of the metric to use as covariances.  First find
*** the inverse of L; then multiply in place for the full upper inverse;
*** then fill in the lower triangle. */
    for (a=0; a<nmod; a++) {
        an = a*nmod;
        *(cvdata+an+a) = 1. / *(cvdata+an+a);
        for (b=a+1; b<nmod; b++) {
            bn = b*nmod;
            sum = 0.;
            for (m=a; m<b; m++) sum -= *(cvdata+bn+m) * *(cvdata+m*nmod+a);
            *(cvdata+bn+a) = sum / *(cvdata+bn+b);
        }
    }
    for (a=0; a<nmod; a++) {
        an = a*nmod;
        for (b=a; b<nmod; b++) {
            sum = 0.;
            for (m=b; m<nmod; m++) sum += *(cvdata+m*nmod+a) * *(cvdata+m*nmod+b);
            *(cvdata+an+b) = sum;
        }
    }
    for (a=0; a<nmod; a++) {
        an = a*nmod;
        for (b=a+1; b<nmod; b++) {
            *(cvdata+b*nmod+a) = *(cvdata+an+b);
        }
    }

/** Clean up and return the covariance matrix. */
    Py_DECREF(cfmetric);
    return PyArray_Return(covar);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(cfmetric);
    Py_XDECREF(covar);
    return NULL;
}


/*******************************************************************************
* Methods table; must have an element for each new method
* and end with the "sentinel" element.
* It should have 4 entries for each method: 
*         {call name, function to call, argument type, doc string},
* with each entry ending with a comma.
*******************************************************************************/

static PyMethodDef methods[] = {
  {"metricAnalysis", vba_metricAnalysis, METH_VARARGS, vba_metricAnalysis_doc},
  {"covar", vba_covar, METH_VARARGS, vba_covar_doc},
  {NULL,		NULL, 0}		/* sentinel */
};


/*******************************************************************************
* The initialization function---it must be named for the module.
* This should be the only non-static global object in this file.
*******************************************************************************/

void initvba(void) {
  PyObject *m;
  
  /* Create the module and add the functions */
  m = Py_InitModule("vba", methods);
  import_array();

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module vba");
}

