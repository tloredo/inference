
/*******************************************************************************
* fkeplermodule: "Fast KEPLER"
*
* Solves the transcendental Kepler equation, using a combination of Newton-
* Raphson iteration and bisection.  "Fast" in comparison with its original
* pure Python implementation.  Methods are provided for finding a single
* solution, and for finding a vector of solutions (returned as a SciPy array).
*
* Created:  07 Apr 2000  Tom Loredo
* Last mod: 09 Apr 2000  TL
*           12 May 2000  TL - Cleaned up error handling
*           21 Dec 2005  TL - Converted to numpy
*******************************************************************************/

/*******************************************************************************
* Headers, prototypes, and macros.
*******************************************************************************/

/* Include Python headers; needed for *all* extensions. */
#include "Python.h"

/* Other includes needed for this extension. */
#include "numpy/arrayobject.h"  /* Header needed to use SciPy array types. */
#include <math.h>

/* Prototypes for functions defined and used in this extension file that aren't 
   taken care of in an include file. */
static void t2TA (double t, double *c, double *s);
static void t2EA (double t, double *c, double *s);
static void t2anoms (double t, double *cea, double *sea, 
                     double *cta, double *sta);
void initfkepler (void);

/* Error handler macro for use in any extension; adapted from Dubois & Yang. */
#define Py_TRY(BOOLEAN) {if (!(BOOLEAN)) goto FAIL;}

/* Macro to cast SciPy array data to a 1-d double array. */
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))


/*------------------------------------------------------------------------------
* Some handy error handling functions, from lapack_litemodule.c in old Numeric.
------------------------------------------------------------------------------*/
static PyObject *ErrorObject;

static PyObject *ErrorReturn(char *mes)
{  if (!ErrorObject)
      ErrorObject = PyString_FromString("fkeplerError");
   PyErr_SetString(ErrorObject,mes);
   return NULL;
}

/* So far not used here.... */
static int NumPy_CheckObject(PyObject *ob, int t, char *obname,
        char *tname, char *funname)
{       char buf[255];
        if (! PyArray_Check(ob))
        {       sprintf(buf,"Expected an array for parameter %s in fkepler.%s",obname, funname);
        ErrorReturn(buf);
        return 0;
    }
    if (!(((PyArrayObject *)ob)->flags & NPY_CONTIGUOUS)) 
    {    sprintf(buf,"Parameter %s is not contiguous in fkepler.%s",obname, funname);
        ErrorReturn(buf);
        return 0;
    }
    if (!(((PyArrayObject *)ob)->descr->type_num == t))
    {    sprintf(buf,"Parameter %s is not of type %s in fkepler.%s",obname,tname, funname);
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

/* These define the orbit.  Either Tp (time of periapsis passage or 
   M0 (mean anomaly at t=0) may be used as the phase origin parameter
   (Mp = -M0 is deprecated). */
static double tau = 1., e = 0., Tp = 0., M0 = 0., Mp = 0.;  

/* Tolerance for eccentric anomaly solution (rad). */
static double err = 1.e-8;

static double pi = 3.1415926535897931;
static double twopi = 6.2831853071795862;

/*------------------------------------------------------------------------------
* Set the globals defining the orbit.  There are two versions, for
* different choices of the phase parameter.
------------------------------------------------------------------------------*/
static char fkepler_setup_Tp_doc[] =
"setup_Tp(tau, e, Tp):  Define the Keplerian orbit.";

static PyObject *
fkepler_setup_Tp (PyObject *self, PyObject *args) {

/*--- Parse the 3 double arguments. */
    Py_TRY(PyArg_ParseTuple(args, "ddd", &tau, &e, &Tp));
    Mp = twopi*Tp/tau;
    M0 = -Mp;

    /* Return None, counting the new reference. */
    Py_INCREF(Py_None);
    return Py_None;

/** Misbehavior ends up here! */
FAIL:
    return NULL;
}


static char fkepler_setup_M0_doc[] =
"setup_M0(tau, e, M0):  Define the Keplerian orbit.";

static PyObject *
fkepler_setup_M0 (PyObject *self, PyObject *args) {

/*--- Parse the 3 double arguments. */
    Py_TRY(PyArg_ParseTuple(args, "ddd", &tau, &e, &M0));
    Tp = -tau*M0/twopi;

    /* Return None, counting the new reference. */
    Py_INCREF(Py_None);
    return Py_None;

/** Misbehavior ends up here! */
FAIL:
    return NULL;
}


static char fkepler_setup_Mp_doc[] =
"setup_Mp(tau, e, Mp):  Define the Keplerian orbit.  Deprecated; use setup_M0.";

static PyObject *
fkepler_setup_Mp (PyObject *self, PyObject *args) {

/*--- Parse the 3 double arguments. */
    Py_TRY(PyArg_ParseTuple(args, "ddd", &tau, &e, &Mp));
    Tp = tau*Mp/twopi;
    M0 = -Mp;

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
    PySys_WriteStdout("Set fkepler tolerance:  %e\n", err);

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
    double f, cf, hi, lo, EA, W, EA_p, dE;
    int i;
    double imodf;

    const int maxit = 20;

/*--- Parse the double argument. */
    Py_TRY(PyArg_ParseTuple(args, "d", &t));

/*--- Convert time to phase mod 2*pi. */
    f = modf((t-Tp)/tau, &imodf);
    /* printf("phase, f, imodf: %g %g %g\n", (t-Tp)/tau, f, imodf); */
    if (f > 0.5) f = f - 1.;
    else if (f < -0.5) f = f + 1.;
    f = twopi*f;

/*--- Initial guess and bracket for root. */
/* O(e**3) guess: */
    cf = cos(f);
    EA = f + e*sin(f)*(1. + e*(cf +0.5*e*(3*cf*cf - 1.)));
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
    double f, cf, hi, lo, EA, W, EA_p, dE;
    int i;
    double imodf;
    double c, s, th, th2;

    const int maxit = 20;

/*** Parse the double argument. */
    Py_TRY(PyArg_ParseTuple(args, "d", &t));

/*** First, convert time to eccentric anomaly.  We just copy the t2EA
**** code here to avoid overhead. */

/*--- Convert time to phase mod 2*pi. */
    f = modf((t-Tp)/tau, &imodf);
    if (f > 0.5) f = f - 1.;
    else if (f < -0.5) f = f + 1.;
    f = twopi*f;

/*--- Initial guess and bracket for root. */
/* O(e**3) guess: */
    cf = cos(f);
    EA = f + e*sin(f)*(1. + e*(cf +0.5*e*(3*cf*cf - 1.)));
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
    double c, s;
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

/** Create a n x 2 array for the result. */
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

/** Clean up and return a SciPy array. */
    Py_DECREF(t_array);
    return PyArray_Return(EA_array);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(t_array);
    Py_XDECREF(EA_array);
    return NULL;
}

/*------------------------------------------------------------------------------
* Calculate eccentric anomaly for a vector of times; store cosines & sines in
* a provided array.
------------------------------------------------------------------------------*/
static char fkepler_t2EA_store_doc[] =
"t2EA_store(t,store):  Convert times to eccentric anomaly; store cosines & sines.";

/*..............................................................................
. Here is the method itself.
..............................................................................*/
static PyObject *
fkepler_t2EA_store (PyObject *self, PyObject *args) {
    PyObject *input1, *input2;
    PyArrayObject *t_array;
    PyArrayObject *EA_array;
    int n, i;
    double c, s;
    double *t_vec, *EA_vec;

/** Parse the arguments; make them contiguous if necessary. */
    Py_TRY(PyArg_ParseTuple(args, "OO", &input1, &input2));
    t_array = (PyArrayObject *)
            PyArray_ContiguousFromAny(input1, PyArray_DOUBLE, 1, 1);
    if (t_array == NULL) goto FAIL;
    EA_array = (PyArrayObject *)
            PyArray_ContiguousFromAny(input2, PyArray_DOUBLE, 2, 2);
    if (EA_array == NULL) goto FAIL;

/** Check consistency of arrays. */
    n = t_array->dimensions[0];
    if (EA_array->dimensions[0] != n) goto FAIL;
    if (EA_array->dimensions[1] != 2) goto FAIL;

/** Step through the times. */
    t_vec = DDATA(t_array);
    EA_vec = DDATA(EA_array);
    for (i=0; i<n; i++) {
        t2EA(*(t_vec+i), &c, &s);
        *(EA_vec + 2*i) = c;
        *(EA_vec + 2*i + 1) = s;
    }

/** Clean up and return. */
    Py_DECREF(t_array);
    Py_DECREF(EA_array);
    Py_INCREF(Py_None);
    return Py_None;

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
    double f, cf, hi, lo, EA, W, EA_p, dE;
    int i;
    double imodf;

    const int maxit = 20;

/*** First, convert time to eccentric anomaly.  We just copy the t2EA
**** code here to avoid overhead. */

/*--- Convert time to phase mod 2*pi. */
    f = modf((t-Tp)/tau, &imodf);
    if (f > 0.5) f = f - 1.;
    else if (f < -0.5) f = f + 1.;
    f = twopi*f;

/*--- Initial guess and bracket for root. */
/* O(e**3) guess: */
    cf = cos(f);
    EA = f + e*sin(f)*(1. + e*(cf +0.5*e*(3*cf*cf - 1.)));
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
    double c, s;
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

/** Clean up and return a SciPy array. */
    Py_DECREF(t_array);
    return PyArray_Return(TA_array);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(t_array);
    Py_XDECREF(TA_array);
    return NULL;
}

/*------------------------------------------------------------------------------
* Calculate true anomaly for a vector of times; store cosines & sines in
* a provided array.
------------------------------------------------------------------------------*/
static char fkepler_t2TA_store_doc[] =
"t2TA_store(t,store):  Convert times to true anomaly; store cosines & sines.";

/*..............................................................................
. Here is the method itself.
..............................................................................*/
static PyObject *
fkepler_t2TA_store (PyObject *self, PyObject *args) {
    PyObject *input1, *input2;
    PyArrayObject *t_array;
    PyArrayObject *TA_array;
    int n, i;
    double c, s;
    double *t_vec, *TA_vec;

/** Parse the arguments; make them contiguous if necessary. */
    Py_TRY(PyArg_ParseTuple(args, "OO", &input1, &input2));
    t_array = (PyArrayObject *)
            PyArray_ContiguousFromAny(input1, PyArray_DOUBLE, 1, 1);
    if (t_array == NULL) goto FAIL;
    TA_array = (PyArrayObject *)
            PyArray_ContiguousFromAny(input2, PyArray_DOUBLE, 2, 2);
    if (TA_array == NULL) goto FAIL;

/** Check consistency of arrays. */
    n = t_array->dimensions[0];
    if (TA_array->dimensions[0] != n) goto FAIL;
    if (TA_array->dimensions[1] != 2) goto FAIL;

/** Step through the times. */
    t_vec = DDATA(t_array);
    TA_vec = DDATA(TA_array);
    for (i=0; i<n; i++) {
        t2TA(*(t_vec+i), &c, &s);
        *(TA_vec + 2*i) = c;
        *(TA_vec + 2*i + 1) = s;
    }

/** Clean up and return. */
    Py_DECREF(t_array);
    Py_DECREF(TA_array);
    Py_INCREF(Py_None);
    return Py_None;

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
    double cf, th, th2;

    const int maxit = 20;

/*** First, convert time to eccentric anomaly.  We just copy the t2EA
**** code here to avoid overhead. */

/*--- Convert time to phase mod 2*pi. */
    f = modf((t-Tp)/tau, &imodf);
    if (f > 0.5) f = f - 1.;
    else if (f < -0.5) f = f + 1.;
    f = twopi*f;

/*--- Initial guess and bracket for root. */
/* There was only ~10% speed difference among these on a 733MHz G4. */

/* O(0) guess: */
/*     EA = f; */

/* O(e) guess: */
/*     EA = f + e*sin(f); */

/* O(e**2) guess: */
/*     EA = f + e*(sin(f) + 0.5*e*sin(2*f)); */

/* O(e**3) guess: */
    cf = cos(f);
    EA = f + e*sin(f)*(1. + e*(cf +0.5*e*(3*cf*cf - 1.)));

/* O(e**3) guess: */
/*     rf = e*sin(f)/(1.-e*cos(f)); */
/*     EA = f + rf*(1. - 0.5*rf*rf); */

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
    PySys_WriteStdout("t2EA did not converge in t2TA!"); /* *** change to an exception! *** */
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


/*------------------------------------------------------------------------------
* Calculate eccentric anomaly for a vector of times; store cosines & sines in
* a provided array.
------------------------------------------------------------------------------*/
static char fkepler_t2anoms_doc[] =
"t2anoms(t,ea,ta):  Calculate eccentric & true anomaly for times; store cosines & sines.";

/*..............................................................................
. Here is the method itself.
..............................................................................*/
static PyObject *
fkepler_t2anoms (PyObject *self, PyObject *args) {
    PyObject *input1, *input2, *input3;
    PyArrayObject *t_array;
    PyArrayObject *EA_array;
    PyArrayObject *TA_array;
    int n, i;
    double c, s, cta, sta;
    double *t_vec, *EA_vec, *TA_vec;

/** Parse the arguments; make them contiguous if necessary. */
    Py_TRY(PyArg_ParseTuple(args, "OOO", &input1, &input2, &input3));
    t_array = (PyArrayObject *)
            PyArray_ContiguousFromAny(input1, PyArray_DOUBLE, 1, 1);
    if (t_array == NULL) goto FAIL;
    EA_array = (PyArrayObject *)
            PyArray_ContiguousFromAny(input2, PyArray_DOUBLE, 2, 2);
    if (EA_array == NULL) goto FAIL;
    TA_array = (PyArrayObject *)
            PyArray_ContiguousFromAny(input3, PyArray_DOUBLE, 2, 2);
    if (TA_array == NULL) goto FAIL;

/** Check consistency of arrays. */
    n = t_array->dimensions[0];
    if (EA_array->dimensions[0] != n) goto FAIL;
    if (EA_array->dimensions[1] != 2) goto FAIL;
    if (TA_array->dimensions[0] != n) goto FAIL;
    if (TA_array->dimensions[1] != 2) goto FAIL;

/** Step through the times. */
    t_vec = DDATA(t_array);
    EA_vec = DDATA(EA_array);
    TA_vec = DDATA(TA_array);
    for (i=0; i<n; i++) {
        /* PySys_WriteStdout("t2anoms: %i %f\n", i, *(t_vec+i)); */
        t2anoms(*(t_vec+i), &c, &s, &cta, &sta);
        *(EA_vec + 2*i) = c;
        *(EA_vec + 2*i + 1) = s;
        *(TA_vec + 2*i) = cta;
        *(TA_vec + 2*i + 1) = sta;
    }

/** Clean up and return. */
    Py_DECREF(t_array);
    Py_DECREF(EA_array);
    Py_DECREF(TA_array);
    Py_INCREF(Py_None);
    return Py_None;

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(t_array);
    Py_XDECREF(EA_array);
    Py_XDECREF(TA_array);
    return NULL;
}

/*..............................................................................
. This "helper" function does the actual computations; it isn't directly
. accessed by Python.
..............................................................................*/
static void t2anoms (double t, double *cea, double *sea, 
                     double *cta, double *sta) {
    double f, hi, lo, EA, W, EA_p, dE;
    int i;
    double imodf;
    double cf, th, th2;

    const int maxit = 20;

/*** First, convert time to eccentric anomaly.  We just copy the t2EA
**** code here to avoid overhead. */

/*--- Convert time to phase mod 2*pi. */
    f = modf((t-Tp)/tau, &imodf);
    if (f > 0.5) f = f - 1.;
    else if (f < -0.5) f = f + 1.;
    f = twopi*f;

/*--- Initial guess and bracket for root. */
/* There was only ~10% speed difference among these on a 733MHz G4. */

/* O(0) guess: */
/*     EA = f; */

/* O(e) guess: */
/*     EA = f + e*sin(f); */

/* O(e**2) guess: */
/*     EA = f + e*(sin(f) + 0.5*e*sin(2*f)); */

/* O(e**3) guess: */
    cf = cos(f);
    EA = f + e*sin(f)*(1. + e*(cf +0.5*e*(3*cf*cf - 1.)));

/* O(e**3) guess: */
/*     rf = e*sin(f)/(1.-e*cos(f)); */
/*     EA = f + rf*(1. - 0.5*rf*rf); */

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
    PySys_WriteStdout("t2anoms did not converge in t2TA!"); /* *** change to an exception! *** */
    *cea = -1.;
    *sea = -1.;
    *cta = -1.;
    *sta = -1.;
    return;

/*** Now get cos & sin of TA directly from tan(phi/2). */
TACalc:
	*cea = cos(EA);
	*sea = sin(EA);
    if (EA == 0.) {
        *cta = 1.;
        *sta = 0.;
    }
    else if (abs(EA) == pi) {
        *cta = -1.;
        *sta = 0.;
    }
    else {
        th = sqrt((1.+e)/(1.-e)) * tan(0.5*EA);
        th2 = th*th;
        *cta = (1.-th2)/(1.+th2);
        *sta = (1.-*cta)/th;
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
  {"setup_Tp", fkepler_setup_Tp, METH_VARARGS, fkepler_setup_Tp_doc},
  {"setup_M0", fkepler_setup_M0, METH_VARARGS, fkepler_setup_M0_doc},
  {"setup_Mp", fkepler_setup_Mp, METH_VARARGS, fkepler_setup_Mp_doc},
  {"settol", fkepler_settol, METH_VARARGS, fkepler_settol_doc},
  {"t2EA", fkepler_t2EA, METH_VARARGS, fkepler_t2EA_doc},
  {"t2TA", fkepler_t2TA, METH_VARARGS, fkepler_t2TA_doc},
  {"vt2TA", fkepler_vt2TA, METH_VARARGS, fkepler_vt2TA_doc},
  {"vt2EA", fkepler_vt2EA, METH_VARARGS, fkepler_vt2EA_doc},
  {"t2EA_store", fkepler_t2EA_store, METH_VARARGS, fkepler_t2EA_store_doc},
  {"t2TA_store", fkepler_t2TA_store, METH_VARARGS, fkepler_t2TA_store_doc},
  {"t2anoms", fkepler_t2anoms, METH_VARARGS, fkepler_t2anoms_doc},
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

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module fkepler");
}

