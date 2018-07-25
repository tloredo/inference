/* Header file for discrete random variable generation with Marsaglia's compact
 * 5-table algorithm. 
 *
 * Copyright: Ed Schofield, 2005-6
 * License: BSD-style (see LICENSE.txt at root of scipy tree)
 *
 * Altered for ppssample by TJL:  Removed handling of PRNG state; added
 * more explicit state info for table allowing simpler save/restore.
 */

#include <stdint.h>
#include "randomkit.h"

typedef struct sampler_t
{
    int32_t t1,t2,t3,t4;           /* limits for table lookups */
    int32_t *AA,*BB,*CC,*DD,*EE;   /* tables for condensed table-lookup */
    int32_t sizeAA, sizeBB, sizeCC, sizeDD, sizeEE; /* table sizes */
    long prob1event;     /* If this is >= 0, this sample point x=prob1event has
                          * p(x) = 1, so the lookup table sampling method won't
                          * work.  (We just let Dran() spit out this value
                          * instead.)  */
} Sampler;

/* Represent probabilities as 30-bit integers and create the 5 tables.
 * The prob mass function is specified as n values p(j) = weights[j]. */
Sampler* init_sampler5tbl(double* weights, unsigned long n);

/* Deallocate it */
void destroy_sampler5tbl(Sampler*);

/* Discrete random variable generating functions using 5 compact tables */
unsigned long Dran(Sampler* sampler, rk_state* rng_state);  /* for a single variate */
void Dran_array(Sampler*, unsigned long* output, unsigned long samplesize,
    rk_state* rng_state);  /* array version */
