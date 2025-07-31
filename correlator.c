#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <complex.h>
#include "utils.h"
#include "filter.h" 


typedef struct
{
    unsigned long int ind_initial;
    unsigned long int ind_final;
    double *psi_acf;
} correlator_t;


correlator_t* correlator_new()
{
    return (correlator_t*) malloc(sizeof(correlator_t));
}


void correlator_ctor(correlator_t *correlator, unsigned long int NX)
{
    unsigned long int N_points = NX * NX;

    correlator->ind_initial = NX / 4;
    correlator->ind_final = 3 * NX / 4;

    correlator->psi_acf = (double*) fftw_malloc(N_points * sizeof(double));

    for(size_t r = 0; r < N_points; ++r)
    {
		(correlator->psi_acf)[r] = 0.0;
    }
}

void correlator_collect_statistics(correlator_t *correlator, fftw_complex *psi, unsigned long int NX)
{
    unsigned long int ind_initial, ind_final, N_2;
    unsigned long int R, r, R_r_x, R_r_y, R_r;
    double term;
    double *psi_acf = correlator->psi_acf;

    N_2 = NX / 2;
    ind_initial = correlator->ind_initial;
    ind_final   = correlator->ind_final;


    for (size_t r_y = ind_initial; r_y < ind_final; ++r_y)
    {
        for (size_t r_x = ind_initial; r_x < ind_final; ++r_x)
        {
            r = r_x + NX * r_y;
            term = 0.0;
            
            R_r_y = (ind_initial - r_y + N_2);
            for (size_t R_y = ind_initial; R_y < ind_final; ++R_y)
            {
                R = ind_initial + NX * R_y;

                R_r_x = (ind_initial - r_x + N_2);
                for (size_t R_x = ind_initial; R_x < ind_final; ++R_x)
                {
                    R_r = R_r_x + NX * R_r_y;

                    term += psi[R_r][0] * psi[R][0] + psi[R_r][1] * psi[R][1];

                    ++R;
                    ++R_r_x;
                }
                ++R_r_y;
            }
            psi_acf[r] += term;
        }
    }
}

/*
void correlator_collect_statistics(correlator_t *correlator, fftw_complex *psi, unsigned long int NX)
{
    unsigned long int ind_initial, ind_final, N_points = NX * NX, N_2 = NX / 2;
    unsigned long int R, r, R_r_x, R_r_y, R_r;

    ind_initial = correlator->ind_initial;
    ind_final   = correlator->ind_final;

    double term;

    for (size_t r_y = ind_initial; r_y < ind_final; ++r_y){
        for (size_t r_x = ind_initial; r_x < ind_final; ++r_x){
            r = r_x + NX * r_y;
            term = 0.0;
            for (size_t R_y = ind_initial; R_y < ind_final; ++R_y){
                R_r_y = (R_y - r_y + N_2);
                for (size_t R_x = ind_initial; R_x < ind_final; ++R_x){
                    R_r_x = (R_x - r_x + N_2);
                    R = R_x + NX * R_y;

                    R_r = R_r_x + NX * R_r_y;
                    term += psi[R_r][0] * psi[R][0] + psi[R_r][1] * psi[R][1];
                }
            }
            (correlator->psi_acf)[r] += term;
        }
    }
}
*/

void correlator_calculate_acf(correlator_t *correlator, unsigned long int N_points, unsigned long int realizations_number)
{

    double multiplier = 1.0 / (N_points / 4 * realizations_number);
    double *psi_acf = correlator->psi_acf;
    
    for(size_t r = 0; r < N_points; ++r)
    {
		psi_acf[r] *= multiplier;
    }

}

void correlator_save_GNU(correlator_t *correlator, char* f_name, unsigned long int NX, unsigned long int NY, double LX, double LY, int step)
{
    save_GNU_XY_r(f_name, correlator->psi_acf, NX, NY, LX, LY, step);
}
/*
void correlator_save_GNU(correlator_t *correlator, char* f_name, unsigned long int NX, unsigned long int NY, double LX, double LY, int step){
    unsigned long int ind_initial, ind_final, k_xk_y, flag = 0;
    double *input = correlator->psi_acf;

    ind_initial = correlator->ind_initial;
    ind_final   = correlator->ind_final;;
    double delta_x, delta_y, x, y, LX_2, LY_2;
    FILE* f_p;

    delta_x = LX/(NX);
    delta_y = LY/(NY);
    LX_2 = 0.5*(LX);
    LY_2 = 0.5*(LY);

    if (f_name == NULL) {
            f_p = stdout;
    } else {

            if ((f_p = fopen(f_name, "w+")) == NULL) {
                    printf ("Can't open/create file %s!\n", f_name);
        exit (1);
            }
    }

    printf ("Saving file %s in GNUPLOT mode.\n", f_name);
    for (size_t k_y = ind_initial; k_y < ind_final; ++k_y){
        for (size_t k_x = ind_initial; k_x < ind_final; ++k_x){
            if ((k_x % step == 0) && (k_y % step == 0)) {
                flag = 1;
                k_xk_y = k_x + NX * k_y;
                x = delta_x*k_x - LX_2;
                y = delta_y*k_y - LY_2;

                fprintf (f_p, "%.15e	%.15e	%.15e\n", x, y, input[k_xk_y]);
            }
        }
        if (flag) {
            fprintf (f_p, "\n");
            flag = 0;
        }
    }

    if (f_name != NULL) fclose(f_p);
}
*/

void correlator_dtor(correlator_t *correlator)
{
    fftw_free(correlator->psi_acf);
}
