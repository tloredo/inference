
/*******************************************************************************
* fkeplermodule: "Fast KEPLER"
*
* Solves the transcendental Kepler equation, using a combination of Newton-
* Raphson iteration and bisection.  "Fast" in comparison with its original
* pure Python implementation.  Methods are provided for finding a single
* solution, and for finding a vector of solutions (returned as a NumPy array).
*
* Created:  07 Apr 2000  Tom Loredo
* Last mod: 09 Apr 2000  TL
*           12 May 2000  TL - Cleaned up error handling
*******************************************************************************/

/*******************************************************************************
* Headers, prototypes, and macros.
*******************************************************************************/

/* Include Python headers; needed for *all* extensions. */
#include "Python.h"

/* Other includes needed for this extension. */
#include "numpy/arrayobject.h"  /* Header needed to use Numeric array types. */

/* We need the math standard C library; but according
   to the NumPy folks, for the Mac we have to wrap it up with a bit of
   protection; this is done in "mymath.h" which is in the NumPy include
   directory.  This may only be needed for CFM68k builds. */
#ifdef macintosh
#include "mymath.h"
#else
#include <math.h>
#endif

/* Prototypes for functions defined and used in this extension file that aren't 
   taken care of in an include file. */
static void t2TA (double t, double *c, double *s);
static void t2EA (double t, double *c, double *s);

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

static int NumPy_CheckObject(PyObject *ob, int t, char *obname,
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

/* Hidden globals for this module, modified by the next two functions. */

static double tau = 1., e = 0., T0 = 0.;  /* These define the orbit. */
static double err = 1.e-8;

/*------------------------------------------------------------------------------
* Set the globals defining the orbit.
------------------------------------------------------------------------------*/
static char fkepler_setup_doc[] =
"setup(tau, e, T0):  Define the Keplerian orbit.";

static PyObject *
fkepler_setup (PyObject *self, PyObject *args) {

/*--- Parse the 3 double arguments. */
    Py_TRY(PyArg_ParseTuple(args, "ddd", &tau, &e, &T0));

    /* Return None, counting the new reference. */
    Py_INCREF(Py_None);
    return Py_None;

/** Misbehavior ends up here! */
FAIL:
    return NULL;
}

/*------------------------------------------------------------------------------
* Set the global defining the Newton-Raphson tolerance.
------------------------------------------------------------------------------*/
static char fkepler_settol_doc[] =
"settol(err):  Set the Newton-Raphson tolerance.";

static PyObject *
fkepler_settol (PyObject *self, PyObject *args) {

/*--- Parse the double argument. */
    Py_TRY(PyArg_ParseTuple(args, "d", &err));

    /* Return None, counting the new reference. */
    Py_INCREF(Py_None);
    return Py_None;

/** Misbehavior ends up here! */
FAIL:
    return NULL;
}

/*------------------------------------------------------------------------------
* Convert time to eccentric anomaly.
------------------------------------------------------------------------------*/
static char fkepler_t2EA_doc[] =
"t2EA(t):  Convert time to eccentric anomaly.";

static PyObject *
fkepler_t2EA (PyObject *self, PyObject *args) {
    double t;
    double f, hi, lo, EA, W, EA_p, dE;
    int i;
    double imodf;

    double pi = 3.14159265359;
    const int maxit = 20;

/*--- Parse the double argument. */
    Py_TRY(PyArg_ParseTuple(args, "d", &t));

/*--- Convert time to phase mod 2*pi. */
    f = modf((t-T0)/tau, &imodf);
    /* printf("phase, f, imodf: %g %g %g\n", (t-T0)/tau, f, imodf); */
    if (f > 0.5) f = f - 1.;
    else if (f < -0.5) f = f + 1.;
    f = 2*pi*f;

/*--- Initial guess and bracket for root. */
    EA = f;
    hi = pi;
    lo = - pi;

/*--- Use Newton-Raphson; but if we jump too far, bisect. */
    for (i=1; i<=maxit; i++) {
        W = EA - e*sin(EA) - f;
        if (W < 0.) lo = EA;
        else hi = EA;
        dE = - W/(1. - e*cos(EA));
        EA_p = EA;
        EA = EA_p + dE;
        if (EA < lo || EA > hi) {
            EA = lo + 0.5*(hi-lo);
            dE = EA - EA_p;
        }
        /* printf("%2i %12.8g %12.8g %5.3e\n", i, EA, EA_p, dE); */
        if (fabs(dE) < err || EA == EA_p)
            return Py_BuildValue("d", EA);
    }
    /* printf("t2EA did not converge!"); */
    return Py_BuildValue("d", -99.);

/** Misbehavior ends up here! */
FAIL:
    return NULL;
}

/*------------------------------------------------------------------------------
* Convert time to true anomaly; return cosine & sine.
------------------------------------------------------------------------------*/
static char fkepler_t2TA_doc[] =
"t2TA(t):  Convert time to true anomaly; return cosine & sine.";

static PyObject *
fkepler_t2TA (PyObject *self, PyObject *args) {
    double t;
    double f, hi, lo, EA, W, EA_p, dE;
    int i;
    double imodf;
    double c, s, th, th2;

    double pi = 3.141592653;
    const int maxit = 20;

/*** Parse the double argument. */
    Py_TRY(PyArg_ParseTuple(args, "d", &t));

/*** First, convert time to eccentric anomaly.  We just copy the t2EA
**** code here to avoid overhead. */

/*--- Convert time to phase mod 2*pi. */
    f = modf((t-T0)/tau, &imodf);
    if (f > 0.5) f = f - 1.;
    else if (f < -0.5) f = f + 1.;
    f = 2*pi*f;

/*--- Initial guess and bracket for root. */
    EA = f;
    hi = pi;
    lo = - pi;

/*--- Use Newton-Raphson; but if we jump too far, bisect. */
    for (i=1; i<=maxit; i++) {
        W = EA - e*sin(EA) - f;
        if (W < 0.) lo = EA;
        else hi = EA;
        dE = - W/(1. - e*cos(EA));
        EA_p = EA;
        EA = EA_p + dE;
        if (EA < lo || EA > hi) {
            EA = lo + 0.5*(hi-lo);
            dE = EA - EA_p;
        }
        if (fabs(dE) < err || EA == EA_p) goto TACalc;
    }
    PySys_WriteStdout("t2EA did not converge in t2TA!");
    return NULL;

/*** Now get cos & sin of TA directly from tan(phi/2). */
TACalc:
    if (EA == 0.) {
        c = 1.;
        s = 0.;
    }
    else if (abs(EA) == pi) {
        c = -1.;
        s = 0.;
    }
    else {
        th = sqrt((1.+e)/(1.-e)) * tan(0.5*EA);
        th2 = th*th;
        c = (1.-th2)/(1.+th2);
        s = (1.-c)/th;
    }
    return Py_BuildValue("dd", c, s);

/** Misbehavior ends up here! */
FAIL:
    return NULL;
}


/*------------------------------------------------------------------------------
* Convert a vector of times to eccentric anomaly; return cosines & sines.
------------------------------------------------------------------------------*/
static char fkepler_vt2EA_doc[] =
"vt2EA(t):  Convert times to eccentric anomaly; return cosines & sines.";

/*..............................................................................
. Here is the method itself.
..............................................................................*/
static PyObject *
fkepler_vt2EA (PyObject *self, PyObject *args) {
    PyObject *input;
    PyArrayObject *t_array;
    PyArrayObject *EA_array;
    int n, i;
    int dims[2];
    double t, c, s;
    double *t_vec, *EA_vec;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  Strictly, this
 ** isn't needed for every object. */
     EA_array = NULL;

/** Parse the argument; make it contiguous. */
    Py_TRY(PyArg_ParseTuple(args, "O", &input));
    t_array = (PyArrayObject *)
            PyArray_ContiguousFromObject(input, PyArray_DOUBLE, 1, 1);
    if (t_array == NULL) goto FAIL;

/** Create a 2 x n array for the result. */
    n = t_array->dimensions[0];
    dims[0] = n;
    dims[1] = 2;
    EA_array = (PyArrayObject *)PyArray_FromDims(2, dims, PyArray_DOUBLE);
    if (EA_array == NULL) goto FAIL;

/** Step through the times. */
    t_vec = DDATA(t_array);
    EA_vec = DDATA(EA_array);
    for (i=0; i<n; i++) {
        t2EA(*(t_vec+i), &c, &s);
        *(EA_vec + 2*i) = c;
        *(EA_vec + 2*i + 1) = s;
    }

/** Clean up and return a NumPy array. */
    Py_DECREF(t_array);
    return PyArray_Return(EA_array);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(t_array);
    Py_XDECREF(EA_array);
    return NULL;
}

/*..............................................................................
. This "helper" function does the actual computations; it isn't directly
. accessed by Python.
..............................................................................*/
static void t2EA (double t, double *c, double *s) {
    double f, hi, lo, EA, W, EA_p, dE;
    int i;
    double imodf;
    double th, th2;

    double pi = 3.141592653;
    const int maxit = 20;

/*** First, convert time to eccentric anomaly.  We just copy the t2EA
**** code here to avoid overhead. */

/*--- Convert time to phase mod 2*pi. */
    f = modf((t-T0)/tau, &imodf);
    if (f > 0.5) f = f - 1.;
    else if (f < -0.5) f = f + 1.;
    f = 2*pi*f;

/*--- Initial guess and bracket for root. */
    EA = f;
    hi = pi;
    lo = - pi;

/*--- Use Newton-Raphson; but if we jump too far, bisect. */
    for (i=1; i<=maxit; i++) {
        W = EA - e*sin(EA) - f;
        if (W < 0.) lo = EA;
        else hi = EA;
        dE = - W/(1. - e*cos(EA));
        EA_p = EA;
        EA = EA_p + dE;
        if (EA < lo || EA > hi) {
            EA = lo + 0.5*(hi-lo);
            dE = EA - EA_p;
        }
        if (fabs(dE) < err || EA == EA_p) goto EACalc;
    }
    PySys_WriteStdout("t2EA did not converge in t2EA!");
    *c = -1.;
    *s = -1.;
    return;

/*** Now get cos & sin of TA directly from tan(phi/2). */
EACalc:
	*c = cos(EA);
	*s = sin(EA);
}


/*------------------------------------------------------------------------------
* Convert a vector of times to true anomaly; return cosines & sines.
------------------------------------------------------------------------------*/
static char fkepler_vt2TA_doc[] =
"vt2TA(t):  Convert times to true anomaly; return cosines & sines.";

/*..............................................................................
. Here is the method itself.
..............................................................................*/
static PyObject *
fkepler_vt2TA (PyObject *self, PyObject *args) {
    PyObject *input;
    PyArrayObject *t_array;
    PyArrayObject *TA_array;
    int n, i;
    int dims[2];
    double t, c, s;
    double *t_vec, *TA_vec;

/** Initialize objects to NULL that we may have to XDECREF if a constructor
 ** fails so that XDECREF has an argument to work with.  Strictly, this
 ** isn't needed for every object. */
     TA_array = NULL;

/** Parse the argument; make it contiguous. */
    Py_TRY(PyArg_ParseTuple(args, "O", &input));
    t_array = (PyArrayObject *)
            PyArray_ContiguousFromObject(input, PyArray_DOUBLE, 1, 1);
    if (t_array == NULL) goto FAIL;

/** Create a 2 x n array for the result. */
    n = t_array->dimensions[0];
    dims[0] = n;
    dims[1] = 2;
    TA_array = (PyArrayObject *)PyArray_FromDims(2, dims, PyArray_DOUBLE);
    if (TA_array == NULL) goto FAIL;

/** Step through the times. */
    t_vec = DDATA(t_array);
    TA_vec = DDATA(TA_array);
    for (i=0; i<n; i++) {
        t2TA(*(t_vec+i), &c, &s);
        *(TA_vec + 2*i) = c;
        *(TA_vec + 2*i + 1) = s;
    }

/** Clean up and return a NumPy array. */
    Py_DECREF(t_array);
    return PyArray_Return(TA_array);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(t_array);
    Py_XDECREF(TA_array);
    return NULL;
}

/*..............................................................................
. This "helper" function does the actual computations; it isn't directly
. accessed by Python.
..............................................................................*/
static void t2TA (double t, double *c, double *s) {
    double f, hi, lo, EA, W, EA_p, dE;
    int i;
    double imodf;
    double th, th2;

    double pi = 3.141592653;
    const int maxit = 20;

/*** First, convert time to eccentric anomaly.  We just copy the t2EA
**** code here to avoid overhead. */

/*--- Convert time to phase mod 2*pi. */
    f = modf((t-T0)/tau, &imodf);
    if (f > 0.5) f = f - 1.;
    else if (f < -0.5) f = f + 1.;
    f = 2*pi*f;

/*--- Initial guess and bracket for root. */
    EA = f;
    hi = pi;
    lo = - pi;

/*--- Use Newton-Raphson; but if we jump too far, bisect. */
    for (i=1; i<=maxit; i++) {
        W = EA - e*sin(EA) - f;
        if (W < 0.) lo = EA;
        else hi = EA;
        dE = - W/(1. - e*cos(EA));
        EA_p = EA;
        EA = EA_p + dE;
        if (EA < lo || EA > hi) {
            EA = lo + 0.5*(hi-lo);
            dE = EA - EA_p;
        }
        if (fabs(dE) < err || EA == EA_p) goto TACalc;
    }
    printf("t2EA did not converge in t2TA!");
    *c = -1.;
    *s = -1.;
    return;

/*** Now get cos & sin of TA directly from tan(phi/2). */
TACalc:
    if (EA == 0.) {
        *c = 1.;
        *s = 0.;
    }
    else if (abs(EA) == pi) {
        *c = -1.;
        *s = 0.;
    }
    else {
        th = sqrt((1.+e)/(1.-e)) * tan(0.5*EA);
        th2 = th*th;
        *c = (1.-th2)/(1.+th2);
        *s = (1.-*c)/th;
    }
}


/*******************************************************************************
* The methods table defining how Python accesses this module's methods.  This
* must have an element for each new method and end with the "sentinel" element.
* It should have 4 entries for each method (except the sentinel): 
*         {call name, function to call, argument type, doc string},
* with each entry ending with a comma.
*******************************************************************************/

static PyMethodDef methods[] = {
  {"setup", fkepler_setup, METH_VARARGS, fkepler_setup_doc},
  {"settol", fkepler_settol, METH_VARARGS, fkepler_settol_doc},
  {"t2EA", fkepler_t2EA, METH_VARARGS, fkepler_t2EA_doc},
  {"t2TA", fkepler_t2TA, METH_VARARGS, fkepler_t2TA_doc},
  {"vt2TA", fkepler_vt2TA, METH_VARARGS, fkepler_vt2TA_doc},
  {"vt2EA", fkepler_vt2EA, METH_VARARGS, fkepler_vt2EA_doc},
  {NULL,		NULL, 0}		/* sentinel */
};


/*******************************************************************************
* The initialization function---it must be named for the module ("init" followed
* by the module name).
*
* This should be the only non-static global object in this file.
*******************************************************************************/

void initfkepler() {
  PyObject *m;
  
  /* Create the module and add the functions */
  m = Py_InitModule("fkepler", methods);
  import_array();

  /* 
  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module fkepler");
}

