#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <complex.h>

#include "dataio.h"
#include "filter.h" 


typedef struct
{
    double *psi_acf;
    double *gauss_filter;
	int isempty;
} correlator_t;

int correlator_isempty(correlator_t *correlator)
{
	return correlator->isempty;
}

correlator_t* correlator_new()
{
    return (correlator_t*) malloc(sizeof(correlator_t));
}


void correlator_ctor(correlator_t *correlator, double dx, double R0, unsigned long int N_grid)
{
    correlator->psi_acf = (double*) fftw_malloc(N_grid * N_grid * sizeof(double));
    correlator->gauss_filter = (double*) fftw_malloc(N_grid * N_grid * sizeof(double));
	correlator->isempty = 1;

    double *gauss_filter = correlator->gauss_filter, *psi_acf = correlator->psi_acf;
    unsigned long int R;
    double exp_factor = -dx * dx / (2 * R0 * R0);

    for (size_t i_y = 0; i_y < N_grid; ++i_y){
        for (size_t i_x = 0; i_x < N_grid; ++i_x){
            R = i_x + N_grid * i_y;
            psi_acf[R] = 0,0;
            gauss_filter[R] = exp((i_x * i_x + i_y * i_y) * exp_factor);
        }
    }
}

void correlator_collect_statistics(correlator_t *correlator, fftw_complex *psi, unsigned long int N_grid)
{
    
    double *gauss_filter = correlator->gauss_filter, *psi_acf = correlator->psi_acf;
    unsigned long int N_points = N_grid * N_grid, R, r, R_r;
    double term;

    for (size_t k_y = 0; k_y < N_grid; ++k_y){
        for (size_t k_x = 0; k_x < N_grid; ++k_x){
            r = k_x + N_grid * k_y;
            term = 0.0;
            for (size_t i_y = 0; i_y < N_grid; ++i_y){
                for (size_t i_x = 0; i_x < N_grid; ++i_x){
                    R = i_x + N_grid * i_y;
                    R_r = R + r;
                    if (R_r < N_points)
                        term += (psi[R_r][0] * psi[R][0] + psi[R_r][1] * psi[R][1]) * gauss_filter[R];
                }
            }
            psi_acf[r] += term;
        }
    }
	correlator->isempty = 0;
}


void correlator_calculate_acf(correlator_t *correlator, unsigned long int N_points, unsigned long int realizations_number)
{

    double multiplier = 1.0 / (realizations_number);
    
    for(size_t k = 0; k < N_points; ++k)
    {
		(correlator->psi_acf)[k] *= multiplier;
    }

}

void correlator_save_GNU(correlator_t *correlator, char* f_name, unsigned long int NX, unsigned long int NY, double LX, double LY, int step)
{
    save_GNU_XY_r(f_name, correlator->psi_acf, NX, NY, LX, LY, step);
}

void correlator_dtor(correlator_t *correlator)
{
    fftw_free(correlator->psi_acf);
    fftw_free(correlator->gauss_filter);
}
