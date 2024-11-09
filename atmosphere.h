#ifndef _ATMOSPHERE_H
#define _ATMOSPHERE_H

#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include <fftw3.h>
#include <gsl/gsl_rng.h>

#include "my_thr_lib.h"
#include "dataio.h"

#include "histogram.h"
#include "correlator.h"
#include "sighandler.h"

#define _VERBOSE
#define _SUPER_VERBOSE
#define TIMING

/* What is done first refraction or diffraction */
#define _DIFR_REFR

/*Для 1го процессора*/
#define THREADS_NUMBER	4
#define BLOCK_SIDE_SIZE 64


/* #define MY_RNG_TYPE	gsl_rng_mt19937 */
/* #define MY_RNG_TYPE	gsl_rng_taus2 */
 #define MY_RNG_TYPE	gsl_rng_ranlxd1 
/* #define MY_RNG_TYPE	gsl_rng_ranlxd2 */
/* #define MY_RNG_TYPE	gsl_rng_ranlxs0 */

#define RNG_INIT_BY_TIME 

#define SCALES_FILTER

#define _SAVE_INITIAL_K
#define _SAVE_INITIAL_XY
#define _SAVE_FINAL_K
#define _SAVE_FINAL_XY

#define HISTOGRAM
#define CORRELATOR
#define _MAX_INTENSITY_STAT

/*** definitions for arrays_init.c ***/
void arrays_init (void);
void arrays_free (void);
void plans_init (void);
void plans_destroy (void);
void rng_init (int argc, char** argv);
void rng_destroy (int argc, char** argv);
void arrays_fill (void);

#endif /* _ATMOSPHERE_H */
