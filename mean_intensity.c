#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <complex.h>
#include "utils.h"
#include "filter.h" 


void mean_intensity_collect_statistics(double *mean_I, double *psi_abs_sqr, unsigned long int N_points)
{
    for (size_t r = 0; r < N_points; ++r)
        mean_I[r] += psi_abs_sqr[r];
}


void mean_intensity_scalarize(double *mean_I, double *psi_abs_sqr, unsigned long int N_points, unsigned long int realizations_number)
{
    double multiplier = 1.0 / (realizations_number);

    for (size_t r = 0; r < N_points; ++r)
        mean_I[r] *= multiplier;
}

void mean_intensity_save_GNU(double *mean_I, char* f_name, unsigned long int NX, unsigned long int NY, double LX, double LY, int step)
{
    save_GNU_XY_r(f_name, mean_I, NX, NY, LX, LY, step);
}
