
/*******************************************************************************
* raystat: "RAYleigh STATistic and related functions"
*
* A C extension module for Python.
*
* Created:  30 Aug 2000  Tom Loredo
*           (For SGR project with Adam Kruger & Ira Wasserman)
* Modified: 25 May 2006  TJL
*           Change from Numeric to numpy, other small changes for Inference
*           package.
*     
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
static double lbessi0(double x);

/* Error handler macro for use in any extension; adapted from Dubois & Yang. */
#define Py_TRY(BOOLEAN) {if (!(BOOLEAN)) goto FAIL;}

/* Macro to cast NumPy array data to a 1-d double array. */
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))

/*------------------------------------------------------------------------------
* Some handy error handling functions, from lapack_litemodule.c in NumPy.
------------------------------------------------------------------------------*/
static PyObject *ErrorObject;

static PyObject *ErrorReturn(char *mes)
{  if (!ErrorObject)
      ErrorObject = PyString_FromString("mkarrayError");
   PyErr_SetString(ErrorObject,mes);
   return NULL;
}

#define TRY(E) if( ! (E)) return NULL

static int NumPy_CheckObject(PyObject *ob, int t, char *obname,
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
*******************************************************************************/

/*------------------------------------------------------------------------------
* Calculate the Rayleigh statistic and related quantities at a single frequency.
* Return them in a tuple:
*
* (S, C, P)
*
*    S	= Sum of sin(w*t)
*    C	= Sum of cos(w*t)
*    P	= Rayleigh power
------------------------------------------------------------------------------*/
static char rayleigh_doc[] =
"rayleigh(w, events):  Calculate the Rayleigh statistic and related quantities.";

static PyObject *
rayleigh (PyObject *self, PyObject *args) {
    PyObject *input;
    PyArrayObject *data;
    int ndat, i;
    double *times;
    double w;
    double wt, S, C, P;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  Strictly, this
 ** isn't needed for every object. */
    data = NULL;

/** Parse the args; make the data array contiguous. */
    Py_TRY(PyArg_ParseTuple(args, "dO", &w, &input));
    data = (PyArrayObject *)
            PyArray_ContiguousFromObject(input, PyArray_DOUBLE, 1, 1);
    if (data == NULL) goto FAIL;

/** Store the # of data; make vector references to event times. */
    ndat = data->dimensions[0];
    times = DDATA(data);

/** Calculate the statistics. */
	S = 0.;
	C = 0.;
    for (i=0; i<ndat; i++) {
        wt = w * *(times+i);
        S = S + sin(wt);
        C = C + cos(wt);
    }
    P = (S*S + C*C)/ndat;

/** Clean up and return everything in a tuple. */
    Py_DECREF(data);
    return Py_BuildValue("fff", S, C, P);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(data);
    return NULL;
}

/*------------------------------------------------------------------------------
* Calculate the squared Rayleigh statistic on an evenly spaced grid of 
* angular frequencies using trigonometric recurrences to speed things up.
*
------------------------------------------------------------------------------*/
static char rsgrid_doc[] =
"wvals, powers = rsgrid(n,w_l,w_u,events):  Calculate the squared Rayleigh\n"\
"statistic on a uniformly spaced frequency grid of n points from w_l to w_u.";

static PyObject *
rsgrid (PyObject *self, PyObject *args) {
    PyObject *input;
    PyArrayObject *data, *P, *wvals;
    int ndat, nw, i, j, Pdims[1];
    double *times, *Pdata, *wdata;
    double *sd, *cd1, *swt, *cwt;
    double w_l, w_u, dw, dwt, wt, w, temp;
    double S, C;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  Strictly, this
 ** isn't needed for every object. */
    data = NULL;
    P = NULL;
    wvals = NULL;
    sd = NULL;
    cd1 = NULL;
    swt = NULL;
    cwt = NULL;

/** Parse the args; make the data array contiguous. */
    Py_TRY(PyArg_ParseTuple(args, "iddO", &nw, &w_l, &w_u, &input));
    data = (PyArrayObject *)
            PyArray_ContiguousFromObject(input, PyArray_DOUBLE, 1, 1);
    if (data == NULL) goto FAIL;

/** Store the # of data; make vector references to event times. */
    ndat = data->dimensions[0];
    times = DDATA(data);

/** Make space for the trig. */
	sd = (double *) malloc(ndat*sizeof(double));
	cd1 = (double *) malloc(ndat*sizeof(double));
	swt = (double *) malloc(ndat*sizeof(double));
	cwt = (double *) malloc(ndat*sizeof(double));

/** Setup for the trig recurrences we'll use in the frequency loop. */
	dw = (w_u - w_l) / (nw-1);
    for (j=0; j<ndat; j++) {
        dwt = dw * *(times+j);
	    sd[j] = sin(dwt);
	    cd1[j] = sin(0.5*dwt);
	    cd1[j] = -2.*cd1[j]*cd1[j];
	    wt = w_l * *(times+j);
	    swt[j] = sin(wt);
	    cwt[j] = cos(wt);
	}

/** Loop over frequencies. */
	Pdims[0] = nw;
    P = (PyArrayObject *)PyArray_FromDims(1, Pdims, PyArray_DOUBLE);
    if (P == NULL) goto FAIL;
    Pdata = DDATA(P);
    wvals = (PyArrayObject *)PyArray_FromDims(1, Pdims, PyArray_DOUBLE);
    if (wvals == NULL) goto FAIL;
    wdata = DDATA(wvals);
    w = w_l;
	for (i=0; i<nw; i++) {
	    wdata[i] = w;
	    w += dw;
		S = C = 0.;
		
		/* Loop over data to calculate the Rayleigh power;
		   also, use some trig to prepare for the next frequency. */
        for (j=0; j<ndat; j++) {
            S = S + swt[j];
            C = C + cwt[j];

		    temp = cwt[j];
		    cwt[j] = cwt[j] + (cwt[j]*cd1[j] - swt[j]*sd[j]);
		    swt[j] = swt[j] + (swt[j]*cd1[j] + temp*sd[j]);
        }
        Pdata[i] = (S*S + C*C)/ndat;
	}

/** Clean up and return everything in a tuple. */
    Py_DECREF(data);
    free(sd);
    free(cd1);
    free(swt);
    free(cwt);
    return Py_BuildValue("NN", PyArray_Return(wvals), PyArray_Return(P));

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(data);
    Py_XDECREF(P);
    Py_XDECREF(wvals);
    if (sd != NULL) free(sd);
    if (cd1 != NULL) free(cd1);
    if (swt != NULL) free(swt);
    if (cwt != NULL) free(cwt);
    return NULL;
}


/*------------------------------------------------------------------------------
* Calculate the log marginal likelihood for w & kappa given S and ndata.
------------------------------------------------------------------------------*/
static char lml_wk_doc[] =
"lml_wk(k,n,S):  Calculate the log marginal likelihood for k and w.";

static PyObject *
lml_wk (PyObject *self, PyObject *args) {
    double k, S, lml;
    int n;

/* Parse the double argument. */
    Py_TRY(PyArg_ParseTuple(args, "did", &k, &n, &S));

/* Do the calculation. */
	lml = lbessi0(k*S) - n*lbessi0(k);

/* Return a PyObject. */
    return Py_BuildValue("d", lml);

/** Misbehavior ends up here! */
FAIL:
    return NULL;
}


/*-----------------------------------------------------------------------------*
* Logarithm of the I_0 modified Bessel function, using polynomial
* approximations from Abram. & Stegun.
------------------------------------------------------------------------------*/
static double lbessi0(double x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
		ans=log(ans);
	} else {
		y=3.75/ax;
		ans=(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2)))))))) / sqrt(ax);
		ans=ax + log(ans);
	}
	return ans;
}

/*******************************************************************************
* Methods table; must have an element for each new method
* and end with the "sentinel" element.
* It should have 4 entries for each method: 
*         {call name, function to call, argument type, doc string},
* with each entry ending with a comma.
*******************************************************************************/

static PyMethodDef methods[] = {
  {"rayleigh", rayleigh, METH_VARARGS, rayleigh_doc},
  {"rsgrid", rsgrid, METH_VARARGS, rsgrid_doc},
  {"lml_wk", lml_wk, METH_VARARGS, lml_wk_doc},
  {NULL,		NULL, 0}		/* sentinel */
};


/*******************************************************************************
* The initialization function---it must be named for the module.
* This should be the only non-static global object in this file.
*******************************************************************************/

void init_rayleigh() {
  PyObject *m;
  
  /* Create the module and add the functions */
  m = Py_InitModule("_rayleigh", methods);
  import_array();

  /* 
  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module _rayleigh");
}

