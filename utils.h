#ifndef _utils_H
#define _utils_H
/*
This is library for reading initial and saving data.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <unistd.h>

typedef fftw_complex complex_t;
typedef unsigned long int index_t;


typedef struct {
	unsigned long int NX;
	unsigned long int NY;
	unsigned long int noise_NX;
	unsigned long int noise_NY;
	unsigned long int realizations_number;
	unsigned long int realizations_save_step;
	unsigned long int N_BINS;
	double LX;
	double LY;
	double delta_z;
	double lambda_0;
	double n_0;
	double C_n_sqr;
	double w_0;
	double l_0;
	double L_0;
    double max_intensity_value;
} ProblemConfig;

void read_data_params (char* f_name, unsigned long int *N_X, unsigned long int *N_Y, double *L_X, double *L_Y, double *z);
void read_data_complex (char* f_name, complex_t *output, unsigned long int *N_X, unsigned long int *N_Y, double *L_X, double *L_Y, double *z);
void save_data_complex (char* f_name, complex_t input, unsigned long int N_X, unsigned long int N_Y, double L_X, double L_Y, double z);

void read_data_cr (char* f_name, complex_t *output, unsigned long int *N_X, unsigned long int *N_Y, double *L_X, double *L_Y, double *z);
void save_data_cr (char* f_name, complex_t *input, unsigned long int *N_X, unsigned long int *N_Y, double *L_X, double *L_Y, double *z);

void save_GNU_k_c (char* f_name, complex_t *input, unsigned long int NX, unsigned long int NY, double LX, double LY, int step);
void save_GNU_k_cr (char* f_name, complex_t *input, unsigned long int NX, unsigned long int NY, double LX, double LY, int step);
void save_GNU_real_k_cr (char* f_name, double *input, unsigned long int NX, unsigned long int NY, double LX, double LY, int step);
void save_GNU_XY_c (char* f_name, complex_t *input, unsigned long int NX, unsigned long int NY, double LX, double LY, int step);
void save_GNU_XY_r (char* f_name, double *input, unsigned long int NX, unsigned long int NY, double LX, double LY, int step);
void save_GNU_X_line_c (char* f_name, complex_t *input, unsigned long int NX, unsigned long int NY, double LX, double LY, double Y);

void angle_sqr_avrg_GNU_k_cr (char* f_name, complex_t *input, unsigned long int NX, unsigned long int NY, double LX, double LY);

void save_point_complex_GNU (char* f_name, complex_t input, double index_value1, double index_value2);
void save_point_real_GNU (char* f_name, double input, double index_value1);

void read_conf_file(char* f_name);

void read_rng_state (char *f_name, gsl_rng *rng);
void save_rng_state (char *f_name, gsl_rng *rng);

void append_two_column_ascii_data (char* f_name, unsigned long int index, double value);
void save_GNU_real_array_with_index (char* f_name, double *input, unsigned long int points_number);

void truncate_file (char* f_name);

#endif
