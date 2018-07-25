////////////////////////////////////////////////////////////////////////////////
// File: rombergs_method.c                                                    //
// Routines:                                                                  //
//    Rombergs_Integration_Method                                             //
////////////////////////////////////////////////////////////////////////////////
#include < math.h >                                     // required for fabs()

static const double richardson[] = {  
  3.333333333333333333e-01, 6.666666666666666667e-02, 1.587301587301587302e-02,
  3.921568627450980392e-03, 9.775171065493646139e-04, 2.442002442002442002e-04,
  6.103888176768601599e-05, 1.525902189669642176e-05, 3.814711817595739730e-06,
  9.536752259018191355e-07, 2.384186359449949133e-07, 5.960464832810451556e-08,
  1.490116141589226448e-08, 3.725290312339701922e-09, 9.313225754828402544e-10,
  2.328306437080797376e-10, 5.820766091685553902e-11, 1.455191522857861004e-11,
  3.637978807104947841e-12, 9.094947017737554185e-13, 2.273736754432837583e-13,
  5.684341886081124604e-14, 1.421085471520220567e-14, 3.552713678800513551e-15,
  8.881784197001260212e-16, 2.220446049250313574e-16
};

#define MAX_COLUMNS 1+sizeof(richardson)/sizeof(richardson[0])

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

////////////////////////////////////////////////////////////////////////////////
//  double Rombergs_Integration_Method( double a, double h, double tolerance, //
//                             int max_cols, double (*f)(double), int *err ); //
//                                                                            //
//  Description:                                                              //
//    If T(f,h,a,b) is the result of applying the trapezoidal rule to approx- //
//    imating the integral of f(x) on [a,b] using subintervals of length h,   //
//    then if I(f,a,b) is the integral of f(x) on [a,b], then                 //
//                           I(f,a,b) = lim T(f,h,a,b)                        //
//    where the limit is taken as h approaches 0.                             //
//    The classical Romberg method applies Richardson Extrapolation to the    //
//    limit of the sequence T(f,h,a,b), T(f,h/2,a,b), T(f,h/4,a,b), ... ,     //
//    in which the limit is approached by successively deleting error terms   //
//    in the Euler-MacLaurin summation formula.                               //
//                                                                            //
//  Arguments:                                                                //
//     double a          The lower limit of the integration interval.         //
//     double h          The length of the interval of integration, h > 0.    //
//                       The upper limit of integration is a + h.             //
//     double tolerance  The acceptable error estimate of the integral.       //
//                       Iteration stops when the magnitude of the change of  //
//                       the extrapolated estimate falls below the tolerance. //
//     int    max_cols   The maximum number of columns to be used in the      //
//                       Romberg method.  This corresponds to a minimum       //
//                       integration subinterval of length 1/2^max_cols * h.  //
//     double *f         Pointer to the integrand, a function of a single     //
//                       variable of type double.                             //
//     int    *err       0 if the extrapolated error estimate falls below the //
//                       tolerance; -1 if the extrapolated error estimate is  //
//                       greater than the tolerance and the number of columns //
//                       is max_cols.                                         //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to a +  h.                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Rombergs_Integration_Method( double a, double h, double tolerance,
                               int max_cols, double (*f)(double), int *err ) { 
   
   double upper_limit = a + h;     // upper limit of integration
   double dt[MAX_COLUMNS];         // dt[i] is the last element in column i.
   double integral = 0.5 * ( (*f)(a) + (*f)(a+h) );
   double x, old_h, delta;
   int j,k;

// Initialize err and the first column, dt[0], to the numerical estimate //
// of the integral using the trapezoidal rule with a step size of h.     //
 
   *err = 0;
   dt[0] = 0.5 * h *  ( (*f)(a) + (*f)(a+h) );

// For each possible succeeding column, halve the step size, calculate  //
// the composite trapezoidal rule using the new step size, and up date  //
// preceeding columns using Richardson extrapolation.                   //

   max_cols = min(max(max_cols,0),MAX_COLUMNS);
   for (k = 1; k < max_cols; k++) {
      old_h = h;

                 // Calculate T(f,h/2,a,b) using T(f,h,a,b) //
 
      h *= 0.5;
      integral = 0.0;
      for (x = a + h; x < upper_limit; x += old_h) integral +=  (*f)(x);
      integral = h * integral + 0.5 * dt[0];

         //  Calculate the Richardson Extrapolation to the limit //

      for (j = 0; j < k; j++) {
         delta =  integral - dt[j];
         dt[j] = integral;
         integral += richardson[j] * delta;
      } 

      //  If the magnitude of the change in the extrapolated estimate //
      //  for the integral is less than the preassigned tolerance,    //
      //  return the estimate with err = 0.                           //

      if ( fabs( delta ) < tolerance ) {
         return integral;
      }
      
             //  Store the current esimate in the kth column. //

      dt[k] = integral;
   }

     // The process didn't converge within the preassigned tolerance //
     // using the maximum number of columns designated.              //
     // Return the current estimate of integral and set err = -1.    //
   
   *err = -1;
   return integral;
}