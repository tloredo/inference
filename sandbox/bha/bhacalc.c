
/*******************************************************************************
* bhacalc: "Bayesian Harmonic Analysis CALCulations"
*
* A C extension module for Python.
*
* The core of the Bretthorst algorithm for harmonic sinusoids: calculating the 
* metric and finding its Jacobian and the sufficient statistic.
*
* Created:  12 Apr 2000  Tom Loredo
* Modified: 11 May 2000  TL - Fixed memory leaks
*           12 May 2000  Speedup
*           21 Dec 2005  Convert to NumPy
*
* Copyright 2000, 2005 by Tom Loredo; all rights reserved.
* Only noncommercial use permitted.
*******************************************************************************/

/*******************************************************************************
* Headers, prototypes, and macros.
*******************************************************************************/

/* Include Python headers; needed for *all* extensions. */
#include "Python.h"

/* Other includes needed for this extension. */
#include "numpy/arrayobject.h"
#include <math.h>

/* Prototypes for functions defined and used in this extension file that aren't 
   taken care of in an include file. */
static void sameSums (int n, double gamma, double *c, double *s, double *x);
static void diffSums (int n, double hgp, double hgm, double *c, double *s, double *x12, double *x21);


/* Error handler macro for use in any extension; adapted from Dubois & Yang. */
#define Py_TRY(BOOLEAN) {if (!(BOOLEAN)) goto FAIL;}

/* Macro to cast SciPy array data to a 1-d double array. */
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))

/*------------------------------------------------------------------------------
* Some handy error handling functions, from lapack_litemodule.c in Numeric.
------------------------------------------------------------------------------*/
static PyObject *ErrorObject;

static PyObject *bhaError(void)
{  if (! ErrorObject)
      ErrorObject = PyString_FromString("bhaError");
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

static int SciPy_CheckObject(PyObject *ob, int t, char *obname,
        char *tname, char *funname)
{       char buf[255];
        if (! PyArray_Check(ob))
        {       sprintf(buf,"Expected an array for parameter %s in lapack_dge.%s",obname, funname);
        ErrorReturn(buf);
        return 0;
    }
    if (!(((PyArrayObject *)ob)->flags & CONTIGUOUS)) 
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

/*..............................................................................
. These "helper" functions calculate sums analytically; they aren't directly
. accessed by Python.
..............................................................................*/
static void sameSums (int n, double gamma, double *c, double *s, double *x) {
    double d, cd, hn;

	d = 0.5*sin(n*gamma)/sin(gamma);
	cd = cos((n-1.)*gamma) * d;
	hn = 0.5*n;
	*c = hn + cd;
	*s = hn - cd;
	*x = sin((n-1.)*gamma) * d;
}

static void diffSums (int n, double hgp, double hgm, double *c, double *s, double *x12, double *x21) {
    double dp, cdp, sdp, dm, cdm, sdm;

	dp = 0.5*sin(n*hgp)/sin(hgp);
	cdp = cos((n-1.)*hgp) * dp;
	sdp = sin((n-1.)*hgp) * dp;
	dm = 0.5*sin(n*hgm)/sin(hgm);
	cdm = cos((n-1.)*hgm) * dm;
	sdm = sin((n-1.)*hgm) * dm;
	*c = cdp + cdm;
	*s = -cdp + cdm;
	*x12 = sdp + sdm;
	*x21 = sdp - sdm;
}

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
static char bha_harmonicAnalysis_doc[] =
"Calculate the metric; get various derived quantities from it & the data projections.";

static PyObject *
bha_harmonicAnalysis (PyObject *self, PyObject *args) {
    int mdims[2], pdims[1];
    int nharm, ndat, nmod, a, b, i, j, m, an, bn;
    int i11, i12, i21, i22;
    PyObject *input2;
    PyArrayObject *proj;
    PyArrayObject *metric, *cfmetric, *amps;
    double *mdata, *cfdata, *diag, *pdata, *adata;
    double gamma, g, c, s, x, hgp, hgm, x12, x21, sum, rdet;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  Strictly, this
 ** isn't needed for every object. */
    metric = NULL;
    proj = NULL;
    cfmetric = NULL;
    amps = NULL;
    diag = NULL;  /* This is for malloc/free, not Python. */

/** Parse the number of harmonics (incl. fundamental) & prjxns; make prjxns contiguous. */
    Py_TRY(PyArg_ParseTuple(args, "iidO", &nharm, &ndat, &gamma, &input2));
    proj = (PyArrayObject *)
            PyArray_ContiguousFromObject(input2, PyArray_DOUBLE, 1, 1);
    if (proj == NULL) goto FAIL;

/** Check that dimensions are consistent. */
    if (2*nharm != proj->dimensions[0]) {
        PyErr_SetString(PyExc_ValueError, "Inconsistent projxn dimensions!");
            goto FAIL;
    }

/** Store the # of models; make vector references to proj. */
    nmod = 2*nharm;
    pdata = DDATA(proj);

/** Calculate the metric. */
    mdims[0] = nmod;
    mdims[1] = nmod;
    metric = (PyArrayObject *)PyArray_FromDims(2, mdims, PyArray_DOUBLE);
    if (metric == NULL) goto FAIL;
    mdata = DDATA(metric);
	for (i=1; i<=nharm; i++) {
		for (j=i; j<=nharm; j++) {
			if (i == j) {
				g = i * gamma;
				sameSums(ndat, g, &c, &s, &x);
				i11 = (2*i-2)*nmod + 2*j-2;
				i12 = (2*i-2)*nmod + 2*j-1;
				i21 = (2*i-1)*nmod + 2*j-2;
				i22 = (2*i-1)*nmod + 2*j-1;
				*(mdata+i11) = c;
				*(mdata+i12) = x;
				*(mdata+i21) = *(mdata+i12);
				*(mdata+i22) = s;
			} else {
				hgp = 0.5 * (i+j) * gamma;
				hgm = 0.5 * (j-i) * gamma;
				diffSums(ndat, hgp, hgm, &c, &s, &x12, &x21);
				i11 = (2*i-2)*nmod + 2*j-2;
				i12 = (2*i-2)*nmod + 2*j-1;
				i21 = (2*i-1)*nmod + 2*j-2;
				i22 = (2*i-1)*nmod + 2*j-1;
				*(mdata+i11) = c;
				*(mdata+i12) = x12;
				*(mdata+i21) = x21;
				*(mdata+i22) = s;
				/* Copy this below the diagonal */
				for (a=2*i-2; a<=2*i-1; a++) {
					for (b=2*j-2; b<=2*j-1; b++) {
						*(mdata+b*nmod+a) = *(mdata+a*nmod+b);
					}
				}
			}
		}
	}

/** Perform a Cholesky factorization of the metric. Calculate
*** the sqrt(determinant) along the way.  */
    cfmetric = (PyArrayObject *)PyArray_Copy(metric);
    if (cfmetric == NULL) goto FAIL;
    cfdata = DDATA(cfmetric);
    diag = (double *) malloc(nmod*sizeof(double));
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
    free(diag);
    Py_DECREF(proj);
    return Py_BuildValue("NNfNf",PyArray_Return(metric),
        PyArray_Return(cfmetric),rdet,
        PyArray_Return(amps),sum);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(metric);
    Py_XDECREF(proj);
    Py_XDECREF(cfmetric);
    Py_XDECREF(amps);
    if (diag != NULL) free(diag);
    return NULL;
}

/*------------------------------------------------------------------------------*/
static char bha_multipletAnalysis_doc[] =
"Calculate the metric; get various derived quantities from it & the data projections.";

static PyObject *
bha_multipletAnalysis (PyObject *self, PyObject *args) {
    int mdims[2], pdims[1];
    int nsin, ndat, nmod, a, b, i, j, m, an, bn;
    int i11, i12, i21, i22;
    PyObject *input1, *input2;
    PyArrayObject *gammas, *proj;
    PyArrayObject *metric, *cfmetric, *amps;
    double *gdata, *pdata, *mdata, *cfdata, *diag, *adata;
    double g, c, s, x, hgp, hgm, x12, x21, sum, rdet;
    /* const double pi = 3.1415926535897931; */

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  Strictly, this
 ** isn't needed for every object. */
    metric = NULL;
    proj = NULL;
    cfmetric = NULL;
    amps = NULL;
    diag = NULL;  /* This is for malloc/free, not Python. */

/** Parse the number of harmonics (incl. fundamental) & prjxns; make prjxns contiguous. */
    Py_TRY(PyArg_ParseTuple(args, "iOO", &ndat, &input1, &input2));
    gammas = (PyArrayObject *)
            PyArray_ContiguousFromObject(input1, PyArray_DOUBLE, 1, 1);
    if (gammas == NULL) goto FAIL;
    proj = (PyArrayObject *)
            PyArray_ContiguousFromObject(input2, PyArray_DOUBLE, 1, 1);
    if (proj == NULL) goto FAIL;

/** Check that dimensions are consistent. */
	nsin = gammas->dimensions[0];
	nmod = proj->dimensions[0];
    if (2*nsin != nmod) {
        PyErr_SetString(PyExc_ValueError, "Inconsistent gamma/prjxn dimensions!");
            goto FAIL;
    }

/** Make vector references to gammas & proj. */
    gdata = DDATA(gammas);
    pdata = DDATA(proj);

/** Calculate the metric. */
    mdims[0] = nmod;
    mdims[1] = nmod;
    metric = (PyArrayObject *)PyArray_FromDims(2, mdims, PyArray_DOUBLE);
    if (metric == NULL) goto FAIL;
    mdata = DDATA(metric);
	for (i=1; i<=nsin; i++) {
		for (j=i; j<=nsin; j++) {
			if (i == j) {
				g = *(gdata+i-1);
				sameSums(ndat, g, &c, &s, &x);
				i11 = (2*i-2)*nmod + 2*j-2;
				i12 = (2*i-2)*nmod + 2*j-1;
				i21 = (2*i-1)*nmod + 2*j-2;
				i22 = (2*i-1)*nmod + 2*j-1;
				*(mdata+i11) = c;
				*(mdata+i12) = x;
				*(mdata+i21) = *(mdata+i12);
				*(mdata+i22) = s;
			} else {
				hgp = 0.5 * (*(gdata+i-1) + *(gdata+j-1));
				hgm = 0.5 * (*(gdata+j-1) - *(gdata+i-1));
				diffSums(ndat, hgp, hgm, &c, &s, &x12, &x21);
				i11 = (2*i-2)*nmod + 2*j-2;
				i12 = (2*i-2)*nmod + 2*j-1;
				i21 = (2*i-1)*nmod + 2*j-2;
				i22 = (2*i-1)*nmod + 2*j-1;
				*(mdata+i11) = c;
				*(mdata+i12) = x12;
				*(mdata+i21) = x21;
				*(mdata+i22) = s;
				/* Copy this below the diagonal */
				for (a=2*i-2; a<=2*i-1; a++) {
					for (b=2*j-2; b<=2*j-1; b++) {
						*(mdata+b*nmod+a) = *(mdata+a*nmod+b);
					}
				}
			}
		}
	}

/** Perform a Cholesky factorization of the metric. Calculate
*** the sqrt(determinant) along the way.  */
    cfmetric = (PyArrayObject *)PyArray_Copy(metric);
    if (cfmetric == NULL) goto FAIL;
    cfdata = DDATA(cfmetric);
    diag = (double *) malloc(nmod*sizeof(double));
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
    free(diag);
    Py_DECREF(proj);
    return Py_BuildValue("NNfNf",PyArray_Return(metric),
        PyArray_Return(cfmetric),rdet,
        PyArray_Return(amps),sum);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(metric);
    Py_XDECREF(proj);
    Py_XDECREF(cfmetric);
    Py_XDECREF(amps);
    if (diag != NULL) free(diag);
    return NULL;
}

/*------------------------------------------------------------------------------
* Calculate covariance matrix for the amplitudes from the Cholesky
* factorization of the metric.
------------------------------------------------------------------------------*/
static char bha_covar_doc[] =
"covar(L):  Calculate the covariance matrix for the amplitudes from the Cholesky factorization of the metric.";

static PyObject *
bha_covar (PyObject *self, PyObject *args) {
    int nmod, a, b, m, an, bn;
    PyObject *input;
    PyArrayObject *cfmetric, *covar;
    double *cvdata;
    double sum;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails, so that XDECREF has an argument to work with.  Strictly, this
 ** isn't needed for every object. */
    cfmetric = NULL;
    covar = NULL;

/** Parse the Cholesky factor & make it contiguous. */
    Py_TRY(PyArg_ParseTuple(args, "O", &input));
    cfmetric = (PyArrayObject *)
            PyArray_ContiguousFromObject(input, PyArray_DOUBLE, 2, 2);
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
  {"harmonicAnalysis", bha_harmonicAnalysis, METH_VARARGS, bha_harmonicAnalysis_doc},
  {"multipletAnalysis", bha_multipletAnalysis, METH_VARARGS, bha_multipletAnalysis_doc},
  {"covar", bha_covar, METH_VARARGS, bha_covar_doc},
  {NULL,		NULL, 0}		/* sentinel */
};


/*******************************************************************************
* The initialization function---it must be named for the module.
* This should be the only non-static global object in this file.
*******************************************************************************/

void initbhacalc(void) {
  PyObject *m;
  
  /* Create the module and add the functions */
  m = Py_InitModule("bhacalc", methods);
  import_array();

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module bhacalc");
}

