/* Sampford's rejection algorithm for finite population sampling with
   probability proportional to size (PPS), *without* replacement, 
   implemented for Python/numpy.

    Created:  June 2006 by Tom Loredo
    License: BSD-style (see LICENSE.txt at root of scipy tree)

*/

#include <stdlib.h>
#include "sampford.h"
#include "randomkit.h"

/* TODO:
    - Use "indices" field, or leave this to Python?
    - Change sample return values from int to long?
    - A more sensible way to handle weight/ratio ordering? 
*/

/* Create object & store the cumulative weights.  We'll let numpy
   do the work of sorting and accumulating before calling this. */
Sampford* 
init_sampford(int npopn, double* cumwts, int* indices) {
    Sampford* sampford = NULL;
    int i;
    
    /* Allocate space.  Note we don't have a sample size yet. */
    sampford = malloc(sizeof(Sampford));
    sampford->nsamp = 0;
    sampford->ratios = NULL;
    sampford->cumratios = NULL;
    
    /* Copy the weight and index info. */
    sampford->npopn = npopn;
    sampford->cumwts = malloc(npopn*sizeof(double));
    sampford->indices = malloc(npopn*sizeof(int));
    for (i=0; i<npopn; i++) {
        sampford->cumwts[i] = cumwts[i];
        sampford->indices[i] = indices[i];
    }
    return sampford;
}


/* Store the weight ratios.  We'll let numpy do the calculations and
   just store the results here. */
void 
set_sample(Sampford* sampford, int nsamp, double* ratios, double* cumratios) {
    int i;
    
    /* Allocate space if needed. */
    if (sampford->nsamp != nsamp) {
        if (sampford->nsamp != 0) {
            free(sampford->ratios);
            free(sampford->cumratios);
        }
        sampford->nsamp = nsamp;
        sampford->ratios = malloc(sampford->npopn*sizeof(double));
        sampford->cumratios = malloc(sampford->npopn*sizeof(double));
    }
    
    /* Copy the info. */
    for (i=0; i<sampford->npopn; i++) {
        sampford->ratios[i] = ratios[i];
        sampford->cumratios[i] = cumratios[i];
    }
}


/* Free memory and destroy. */
void 
destroy_sampford(Sampford* sampford) {
    if (sampford->npopn != 0) {
        free(sampford->cumwts);
        free(sampford->indices);
    }
    if (sampford->nsamp != 0) {
        free(sampford->ratios);
        free(sampford->cumratios);
    }
    free(sampford);
}


/* Obtain a sample via Sampford's rejection method. */
int 
sampford_pps(Sampford* sampford, int* sample, rk_state* rng_state_ptr) {
    int i, j, k, duplicate, ntry;
    int npopn, nsamp;
    double val, tot;
    
    npopn = sampford->npopn;
    nsamp = sampford->nsamp;

/* Cycle the rejection loop until nsamp distinct members get selected. */
    ntry = 0;
    do {
        ntry++;

/* Pick the first value based on the raw weights. */
        /* val = sampford->cumwts[npopn-1] * (1.*rand())/(RAND_MAX+1.); */
        val = sampford->cumwts[npopn-1] * rk_double(rng_state_ptr);
        for (i=0; i<npopn; i++) {
            if (val <= sampford->cumwts[i]) break;
        }
        sample[0] = i;
        /* PySys_WriteStdout("First sample unit:  %i  %i\n", i, ntry); */
    
/* Pick subsequent samples based on the ratios, rejecting sequences with duplicates. */
        duplicate = 0;
        tot = sampford->cumratios[npopn-1];
        for (j=1; j<nsamp; j++) {
            val = tot * (1.*rand())/(RAND_MAX+1.);
            for (i=0; i<npopn; i++) {
                if (val <= sampford->cumratios[i]) break;
            }
            /* If i was already selected, note a duplicate and end the ratio loop. */
            for (k=0; k<j; k++) {
                if (sample[k] == i) {
                    duplicate = 1;
                    break;
                }
            }
            /* PySys_WriteStdout("j, choice, dup:  %i  %i  %i\n", j, i, duplicate); */
            if (duplicate) break;
            sample[j] = i;
        }

/* End the rejection loop if there are no duplicates; else start over. */
    } while (duplicate);

    return ntry;
}


/* Obtain a sample via Sampford's rejection method, using 5-table lookups.
   This should be lots faster for large populations. 
   This assumes the input wt_sampler is in *unsorted* order, but
   the ratio sampler is in *sorted* order. */
int 
sampford_pps_tables(Sampford* sampford, Sampler* wt_sampler, 
                    Sampler* ratio_sampler, int* sample, 
                    rk_state* rng_state_ptr) {
    int i, j, k, duplicate, ntry;
    int npopn, nsamp;

/* Consider making these arguments i/o sampford, unless use is
   made of sampford->indices. */
    npopn = sampford->npopn;
    nsamp = sampford->nsamp;

/* Cycle the rejection loop until nsamp distinct members get selected. */
    ntry = 0;
    do {
        ntry++;

/* Pick the first value based on the raw weights. */
        sample[0] = (int) Dran(wt_sampler, rng_state_ptr);
        /* PySys_WriteStdout("First sample unit:  %i  %i\n", sample[0], ntry); */
    
/* Pick subsequent samples based on the ratios, rejecting sequences with duplicates. */
        duplicate = 0;
        for (j=1; j<nsamp; j++) {
            i = (int) Dran(ratio_sampler, rng_state_ptr);
            i = sampford->indices[i];
            /* If i was already selected, note a duplicate and end the ratio loop. */
            for (k=0; k<j; k++) {
                if (sample[k] == i) {
                    duplicate = 1;
                    break;
                }
            }
            /* PySys_WriteStdout("j, choice, dup:  %i  %i  %i\n", j, i, duplicate); */
            if (duplicate) break;
            sample[j] = i;
        }

/* End the rejection loop if there are no duplicates; else start over. */
    } while (duplicate);

    return ntry;
}
