#ifndef _ARRAYS_OPS_H
#define _ARRAYS_OPS_H

void mul_arrays_third_equal_first_by_second (fftw_complex *first, fftw_complex *second, fftw_complex *third, unsigned long int points_number);
void mul_arrays_inplace_first_by_second (fftw_complex *first, fftw_complex *second, unsigned long int points_number);
void mul_arrays_inplace_first_by_real_second (fftw_complex *first, double *second, unsigned long int points_number);

void mul_array_inplace_by_exp_abs_squared (fftw_complex *array, unsigned long int points_number);
void mul_array_inplace_by_real_number (fftw_complex *array, double coeff, unsigned long int points_number);

double L_2_norm (fftw_complex *input, unsigned long int points_number);
void arrays_transpose (fftw_complex *input, fftw_complex *output, unsigned long int row_elements_N, unsigned long int column_elements_M);
void unit_transpose (fftw_complex *input, fftw_complex *output, unsigned long int current_block, unsigned long int row_elements_N, unsigned long int column_elements_M, unsigned long int row_blocks_number, unsigned long int column_blocks_number, unsigned long int row_block_skip, unsigned long int column_block_skip);

void abs_sqr_cdata (fftw_complex *input, double *output, unsigned long int points_number);
void cexp_I_real_calc (double *input, fftw_complex *output, unsigned long int points_number);

void S_k_calc (fftw_complex *rng_pairs, fftw_complex *output);
void rng_pairs_calc (fftw_complex *output, unsigned long int pairs_number);

unsigned long int find_max_index (double *positive_input, unsigned long int points_number);

struct mul_array_third_equal_first_by_second_thread_input {
	fftw_complex *first, *second, *third;
	unsigned long int start_index;
	unsigned long int elements_number;
};

void mul_arrays_third_equal_first_by_second_thr (void *thr_input);

struct mul_arrays_inplace_first_by_second_thread_input {
	fftw_complex *first, *second;
	unsigned long int start_index;
	unsigned long int elements_number;
};

void mul_arrays_inplace_first_by_second_thr (void *thr_input);

struct mul_arrays_inplace_first_by_real_second_thread_input {
	fftw_complex *first;
	double *second;
	unsigned long int start_index;
	unsigned long int elements_number;
};

void mul_arrays_inplace_first_by_real_second_thr (void *thr_input);

struct mul_array_inplace_by_real_number_thread_input {
	fftw_complex *input_array;
	unsigned long int index;
	unsigned long int elements_number;
	double factor;
};

void mul_array_inplace_by_real_number_thr (void *thr_input);

struct L_2_norm_thread_input {
	fftw_complex *input_array;
	unsigned long int index, elements_number;
	double sum;
};

void L_2_norm_thr (void * input_thr);

struct arrays_transpose_thread_input {
	fftw_complex *input, *output;
	unsigned long int index, elements_number, row_elements_N, column_elements_M, row_blocks_number, column_blocks_number, row_block_skip, column_block_skip;
};

void arrays_transpose_thr (void * input_thr);

struct abs_sqr_cdata_thread_input {
	fftw_complex *input_array;
	double *output_array;
	unsigned long int index, elements_number;
};

void abs_sqr_cdata_thr (void * input_thr);

struct cexp_I_real_calc_thread_input {
	fftw_complex *output_array;
	double *input_array;
	unsigned long int index, elements_number;
};

void cexp_I_real_calc_thr (void * input_thr);

/* Atmposphere specific files */

void S_k_calc (fftw_complex *rng_pairs, fftw_complex *output);

struct S_k_calc_thread_input {
	fftw_complex *rng_pairs, *output;
	unsigned long int index, elements_number;
};

void S_k_calc_thr (void *thr_input);

struct rng_pairs_calc_thread_input {
	fftw_complex *output;
	unsigned long int index, elements_number;
	gsl_rng *thread_rng;
};

void rng_pairs_calc_thr (void *thr_input);

#endif /* _ARRAYS_OPS_H */
