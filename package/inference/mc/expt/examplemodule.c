/*******************************************************************************
* examplemodule
*
* Simple examples of a Python extensions, defining a Python method
* that returns the square of a Python double, and functions that
* do some simple interactions with scipy arrays.
*
* Created:  01 Apr 2002 by Tom Loredo
* Modified: 01 Apr 2002
*           20 Dec 2005 converted to scipy
*******************************************************************************/

/*******************************************************************************
* Headers, prototypes, and macros.
*******************************************************************************/
/* Include Python headers; needed for *all* extensions. */
#include "Python.h"

/* Other includes needed for this extension. */
#include "numpy/arrayobject.h"  /* Header needed to use scipy array types. */
#include <math.h>

/* Some notes about the headers above:
* math.h:
*  Mac users note:  math.h will not work on very old (CFM) Macs; in that
*  case you should include "mymath.h" (this could be automated by checking
*  to see if the macro "SYMANTEC__CFM68K__" is defined).
*
* arrayobject.h:
*  Make sure to use the "scipy/" part of the path. */

/* Prototypes for functions defined and used in this extension file that aren't 
   taken care of in an include file. */
/* None are needed for the basic examples here. */

/* Error handler macro for use in any extension; adapted from Dubois & Yang. 
*  This requires that the function using it have the FAIL label defined. */
#define Py_TRY(BOOLEAN) {if (!(BOOLEAN)) goto FAIL;}

/* Macro to cast scipy array data to a 1-d double array; borrowed
   from lapack_litemodule.c in numpy. */
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))

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
* Calculate the square and cube of a nonnegative double (a Python float).
* Negative args raise an exception (just as an example of how to do this!).
------------------------------------------------------------------------------*/
static char example_sqr_cube_doc[] =
"sqr_cube(x):  Return the square and cube of nonnegative x." ;

static PyObject *
example_sqr_cube (PyObject *self, PyObject *args) {
    double x, xx, xxx;
    int i;

/* Parse the double argument. */
    Py_TRY(PyArg_ParseTuple(args, "i", &i));
    PySys_WriteStdout("Int input: %i\n",i);
    x = (double) i;
    if (x<0.) {
    	PyErr_SetString(PyExc_ValueError, "Argument must be >=0!");
    	goto FAIL;
    }

/* Do the calculation, possibly by calling a function defined elsewhere. */
    xx = x*x;
    xxx = x*x*x;

/* Return a PyObject. */
    return Py_BuildValue("dd", xx, xxx);	/* Returns a 2-tuple of floats. */

/** Misbehavior ends up here! */
FAIL:
    return NULL;
}


/*------------------------------------------------------------------------------
* Return a 3-vector of doubles (a 1-d scipy array of type "Float").
------------------------------------------------------------------------------*/
static char example_vec3_doc[] =
"vec3():  Return a filled 3-vector.";

static PyObject *
example_vec3 (PyObject *self, PyObject *args) {
    int dims[1];
    PyArrayObject *result;
    double *data;

/* No arguments to parse.... */

/* Create a 3-vector for the result. */
    dims[0] = 3;
    /* printf("creating vector of size %d...\n",dims[0]); */
    result = (PyArrayObject *)PyArray_FromDims(1, dims, PyArray_DOUBLE);
    if (result == NULL) {
    	PyErr_SetString(PyExc_RuntimeError, "Could not make result array!");
    	goto FAIL;
    }

/* Assign values to its 3 elements. */
    if (SciPy_CheckObject((PyObject *)result, PyArray_DOUBLE, "result", "PyArray_DOUBLE", "vec3") == 0) goto FAIL;
    data = DDATA(result);
    data[0] = 1.;
    data[1] = 2.;
    data[2] = 3.;

/* Sometimes in debugging an extension it is handy to write an intermediate
 * result to the console.  Use PySys_WriteStdout or PySys_WriteStderr for
 * this; they work like printf.  Directly calling printf is problematic, particularly 
 * on the Mac where access to a console is not standard.  The output should be
 * limited to 1000 bytes.  See Python/sysmodule.c if details of these 
 * functions are needed.*/
    PySys_WriteStdout("Example debug string: Element %i is %e\n",1,data[1]);


/** Return a SciPy array. */
    return PyArray_Return(result);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(result);		/* Prevents a memory leak if we don't return result. */
    return NULL;
}


/*------------------------------------------------------------------------------
* Return a 3x3 matrix of doubles (a 2-d scipy array of type "Float").
------------------------------------------------------------------------------*/
static char example_mat3_doc[] =
"mat3():  Return a filled 3x3 matrix.";

static PyObject *
example_mat3 (PyObject *self, PyObject *args) {
    int dims[2] = {3, 3};
    int i, j;
    PyArrayObject *result;
    double *data;

/* No arguments to parse.... */

/** Create a 3x3 matrix for the result. */
    result = (PyArrayObject *)PyArray_FromDims(2, dims, PyArray_DOUBLE);
    if (result == NULL) goto FAIL;

/** Assign values to its 3 elements. */
    data = DDATA(result);
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
             *(data+i*3+j) = 10*i + j;	/* (data+i*3+j) points to element [i][j]. */
        }
    }

/** Return a SciPy array. */
    return PyArray_Return(result);

/** Misbehavior ends up here! */
FAIL:
    Py_XDECREF(result);	/* Prevents a memory leak if we don't return result. */
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
  {"sqr_cube", example_sqr_cube, METH_VARARGS, example_sqr_cube_doc},
  {"vec3", example_vec3, METH_VARARGS, example_vec3_doc},
  {"mat3", example_mat3, METH_VARARGS, example_mat3_doc},
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
void initexample(void) {
  PyObject *m;
  
  /* Create the module and add the functions */
  m = Py_InitModule("example", methods);
  import_array();

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module example");
}

