#include <stdint.h>
#include "randomkit.h"
#include "compact5table.h"

typedef struct sampford_t
{
    /* Data that are fixed for the population */
    int npopn;              /* Population size */
    double *cumwts;         /* Cumulative weights, largest to smallest */
    int* indices; /* Indices for item corresponding to each slot */
    
    /* Data that are fixed for a given sample size. */
    int nsamp;              /* Sample size */
    double *ratios;         /* Weight ratios (used for 5-table acceleration) */
    double *cumratios;      /* Cumulative weight ratios */
} Sampford;

/* Create object & store the cumulative weights */
Sampford* init_sampford(int npopn, double* cumwts, int* indices);

/* Store the weight ratios */
void set_sample(Sampford* sampford, int nsamp, double* ratios, double* cumratios);

/* Free memory and destroy */
void destroy_sampford(Sampford* sampford);

/* Obtain a sample. */
int sampford_pps(Sampford* sampford, int* output, rk_state* rng_state_ptr);

/* Obtain a sample, using 5-table lookups. */
int sampford_pps_tables(Sampford* sampford, Sampler* wt_sampler, 
                    Sampler* ratio_sampler, int* sample, 
                    rk_state* rng_state_ptr);
