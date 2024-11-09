#ifndef correlator_h
#define correlator_h

#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <complex.h>


#include "atmosphere.h"
#include "dataio.h"
#include "filter.h" 

struct correlator_t;

int correlator_isempty(struct correlator_t *correlator);
struct correlator_t* correlator_new(void);
void correlator_ctor(struct correlator_t *correlator, double dx, double R0, unsigned long int N_grid);
void correlator_dtor(struct correlator_t *correlator);

void correlator_collect_statistics(struct correlator_t *correlator, fftw_complex *psi, unsigned long int N_grid);
void correlator_calculate_acf(struct correlator_t *correlator, unsigned long int N_points, unsigned long int realizations_number);
void correlator_save_GNU(struct correlator_t *correlator, char* f_name, unsigned long int NX, unsigned long int NY, double LX, double LY, int step);


#endif // correlator_h
