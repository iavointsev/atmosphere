/*** in this file we placed simple operations with arrays ***/

struct mul_array_third_equal_first_by_second_thread_input *thr_mul_array_third_equal_first_by_second;
struct mul_arrays_inplace_first_by_second_thread_input *thr_array_mul_arrays_inplace_first_by_second;
struct mul_arrays_inplace_first_by_real_second_thread_input *thr_array_mul_arrays_inplace_first_by_real_second;
struct mul_array_inplace_by_real_number_thread_input *thr_mul_array_inplace_by_real_number;
struct L_2_norm_thread_input *thr_L_2_norm;
struct arrays_transpose_thread_input *thr_arrays_transpose;
struct abs_sqr_cdata_thread_input *thr_abs_sqr_cdata;
struct cexp_I_real_calc_thread_input *thr_cexp_I_real_calc;

struct S_k_calc_thread_input *thr_S_k_calc;
struct rng_pairs_calc_thread_input *thr_rng_pairs_calc;

void mul_arrays_third_equal_first_by_second (complex_t *first, complex_t *second, complex_t *third, unsigned long int points_number) {
	complex_t *temp_ptr1, *temp_ptr2, *temp_ptr3;
	unsigned long int i;

	temp_ptr1 = first;
	temp_ptr2 = second;
	temp_ptr3 = third;
#ifndef THREADS_NUMBER
	for (i = 0; i < points_number; ++i) {
		(*temp_ptr3)[0] = ((*temp_ptr1)[0] * (*temp_ptr2)[0] - (*temp_ptr1)[1] * (*temp_ptr2)[1]);
		(*temp_ptr3)[1] = ((*temp_ptr1)[1] * (*temp_ptr2)[0] + (*temp_ptr1)[0] * (*temp_ptr2)[1]);

		++temp_ptr1; ++temp_ptr2; ++temp_ptr3;
	}
#else
#ifdef SUPER_VERBOSE
	printf ("mul_arrays_third_equal_first_by_second with threads.\n");
#endif
	each_thread_elements = (unsigned long int) (0.5 + (1.0*points_number)/(THREADS_NUMBER));
	last_thread_elements = points_number - each_thread_elements*((THREADS_NUMBER)-1);

	for (i = 0; i < ((THREADS_NUMBER)-1); ++i) {
		thr_mul_array_third_equal_first_by_second [i].first = first;
		thr_mul_array_third_equal_first_by_second [i].second = second;
		thr_mul_array_third_equal_first_by_second [i].third = third;
		thr_mul_array_third_equal_first_by_second [i].start_index = i*each_thread_elements;
		thr_mul_array_third_equal_first_by_second [i].elements_number = each_thread_elements;

		my_thr_data_assign (i,  (void *) &thr_mul_array_third_equal_first_by_second[i]);
	}

	thr_mul_array_third_equal_first_by_second [((THREADS_NUMBER)-1)].first = first;
	thr_mul_array_third_equal_first_by_second [((THREADS_NUMBER)-1)].second = second;
	thr_mul_array_third_equal_first_by_second [((THREADS_NUMBER)-1)].third = third;
	thr_mul_array_third_equal_first_by_second [((THREADS_NUMBER)-1)].start_index = ((THREADS_NUMBER)-1)*each_thread_elements;
	thr_mul_array_third_equal_first_by_second [((THREADS_NUMBER)-1)].elements_number = last_thread_elements;

	my_thr_data_assign (((THREADS_NUMBER)-1),  (void *) &thr_mul_array_third_equal_first_by_second[((THREADS_NUMBER)-1)]);

	my_thr_manager (mul_arrays_third_equal_first_by_second_thr);
#endif
}

void mul_arrays_third_equal_first_by_second_thr (void *thr_input) {
	unsigned long int i, my_start_index;
	complex_t *temp_ptr1, *temp_ptr2, *temp_ptr3;

	struct mul_array_third_equal_first_by_second_thread_input *input;

	input = (struct mul_array_third_equal_first_by_second_thread_input*) thr_input;

	temp_ptr1 = input->first + input->start_index;
	temp_ptr2 = input->second + input->start_index;
	temp_ptr3 = input->third + input->start_index;
	for (i = 0; i < input->elements_number; ++i) {
		(*temp_ptr3)[0] = ((*temp_ptr1)[0] * (*temp_ptr2)[0] - (*temp_ptr1)[1] * (*temp_ptr2)[1]);
		(*temp_ptr3)[1] = ((*temp_ptr1)[1] * (*temp_ptr2)[0] + (*temp_ptr1)[0] * (*temp_ptr2)[1]);

		++temp_ptr1; ++temp_ptr2; ++temp_ptr3;
	}
}

void mul_arrays_inplace_first_by_second (complex_t *first, complex_t *second, unsigned long int points_number) {
	complex_t *temp_ptr1, *temp_ptr2;
	double temp_re;
	unsigned long int i;

#ifndef THREADS_NUMBER
	temp_ptr1 = first;
	temp_ptr2 = second;
	for (i = 0; i < points_number; ++i) {
		temp_re = (*temp_ptr1)[0];
		(*temp_ptr1)[0] = ((*temp_ptr1)[0] * (*temp_ptr2)[0] - (*temp_ptr1)[1] * (*temp_ptr2)[1]);
		(*temp_ptr1)[1] = ((*temp_ptr1)[1] * (*temp_ptr2)[0] + temp_re * (*temp_ptr2)[1]);

		++temp_ptr1; ++temp_ptr2;
	}
#else
#ifdef SUPER_VERBOSE
	printf ("mul_arrays_inplace_first_by_second with threads.\n");
#endif
	each_thread_elements = (unsigned long int) (0.5 + (1.0*points_number)/(THREADS_NUMBER));
	last_thread_elements = points_number - each_thread_elements*((THREADS_NUMBER)-1);

	for (i = 0; i < ((THREADS_NUMBER)-1); ++i) {
		thr_array_mul_arrays_inplace_first_by_second [i].first = first;
		thr_array_mul_arrays_inplace_first_by_second [i].second = second;
		thr_array_mul_arrays_inplace_first_by_second [i].start_index = i*each_thread_elements;
		thr_array_mul_arrays_inplace_first_by_second [i].elements_number = each_thread_elements;

		my_thr_data_assign (i,  (void *) &thr_array_mul_arrays_inplace_first_by_second[i]);
	}

	thr_array_mul_arrays_inplace_first_by_second [((THREADS_NUMBER)-1)].first = first;
	thr_array_mul_arrays_inplace_first_by_second [((THREADS_NUMBER)-1)].second = second;
	thr_array_mul_arrays_inplace_first_by_second [((THREADS_NUMBER)-1)].start_index = ((THREADS_NUMBER)-1)*each_thread_elements;
	thr_array_mul_arrays_inplace_first_by_second [((THREADS_NUMBER)-1)].elements_number = last_thread_elements;

	my_thr_data_assign (i,  (void *) &thr_array_mul_arrays_inplace_first_by_second[i]);

	my_thr_manager (mul_arrays_inplace_first_by_second_thr);
#endif
}

void mul_arrays_inplace_first_by_second_thr (void *thr_input) {
	unsigned long int i, my_start_index;
	double temp_re;
	complex_t *temp_ptr1, *temp_ptr2;
	struct mul_arrays_inplace_first_by_second_thread_input *input;

	input = (struct mul_arrays_inplace_first_by_second_thread_input*) thr_input;

	temp_ptr1 = input->first + input->start_index;
	temp_ptr2 = input->second + input->start_index;

	for (i = 0; i < input->elements_number; ++i) {
		temp_re = (*temp_ptr1)[0];
		(*temp_ptr1)[0] = (temp_re * (*temp_ptr2)[0] - (*temp_ptr1)[1] * (*temp_ptr2)[1]);
		(*temp_ptr1)[1] = ((*temp_ptr1)[1] * (*temp_ptr2)[0] + temp_re * (*temp_ptr2)[1]);

		++temp_ptr1; ++temp_ptr2;
	}
}

void mul_arrays_inplace_first_by_real_second (complex_t *first, double *second, unsigned long int points_number) {
	complex_t *temp_ptr1;
	double *temp_ptr2;
	unsigned long int i;

#ifndef THREADS_NUMBER
	temp_ptr1 = first;
	temp_ptr2 = second;
	for (i = 0; i < points_number; ++i) {
		((*temp_ptr1)[0]) *= (*temp_ptr2);
		((*temp_ptr1)[1]) *= (*temp_ptr2);

		++temp_ptr1; ++temp_ptr2;
	}
#else
#ifdef SUPER_VERBOSE
	printf ("mul_arrays_inplace_first_by_second with threads.\n");
#endif
	each_thread_elements = (unsigned long int) (0.5 + (1.0*points_number)/(THREADS_NUMBER));
	last_thread_elements = points_number - each_thread_elements*((THREADS_NUMBER)-1);

	for (i = 0; i < ((THREADS_NUMBER)-1); ++i) {
		thr_array_mul_arrays_inplace_first_by_real_second [i].first = first;
		thr_array_mul_arrays_inplace_first_by_real_second [i].second = second;
		thr_array_mul_arrays_inplace_first_by_real_second [i].start_index = i*each_thread_elements;
		thr_array_mul_arrays_inplace_first_by_real_second [i].elements_number = each_thread_elements;

		my_thr_data_assign (i,  (void *) &thr_array_mul_arrays_inplace_first_by_real_second[i]);
	}

	thr_array_mul_arrays_inplace_first_by_real_second [((THREADS_NUMBER)-1)].first = first;
	thr_array_mul_arrays_inplace_first_by_real_second [((THREADS_NUMBER)-1)].second = second;
	thr_array_mul_arrays_inplace_first_by_real_second [((THREADS_NUMBER)-1)].start_index = ((THREADS_NUMBER)-1)*each_thread_elements;
	thr_array_mul_arrays_inplace_first_by_real_second [((THREADS_NUMBER)-1)].elements_number = last_thread_elements;

	my_thr_data_assign (i,  (void *) &thr_array_mul_arrays_inplace_first_by_real_second[i]);

	my_thr_manager (mul_arrays_inplace_first_by_real_second_thr);
#endif
}

void mul_arrays_inplace_first_by_real_second_thr (void *thr_input) {
	unsigned long int i, my_start_index;
	complex_t *temp_ptr1;
	double *temp_ptr2;
	struct mul_arrays_inplace_first_by_real_second_thread_input *input;

	input = (struct mul_arrays_inplace_first_by_real_second_thread_input*) thr_input;

	temp_ptr1 = input->first + input->start_index;
	temp_ptr2 = input->second + input->start_index;

	for (i = 0; i < input->elements_number; ++i) {
		((*temp_ptr1)[0]) *= (*temp_ptr2);
		((*temp_ptr1)[1]) *= (*temp_ptr2);

		++temp_ptr1; ++temp_ptr2;
	}
}

void mul_array_inplace_by_real_number (complex_t *array, double coeff, unsigned long int points_number) {
	complex_t *temp_ptr1;
	unsigned long int i;

#ifndef THREADS_NUMBER
	temp_ptr1 = array;
	for (i = 0; i < points_number; ++i) {
		(*temp_ptr1)[0] *= coeff;
		(*temp_ptr1)[1] *= coeff;

		++temp_ptr1;
	}
#else
#ifdef SUPER_VERBOSE
	printf ("mul_array_inplace_by_real_number with threads.\n");
#endif
	unsigned long int each_thread_elements = (unsigned long int) (0.5 + (1.0*points_number)/(THREADS_NUMBER));
	unsigned long int last_thread_elements = points_number - each_thread_elements*((THREADS_NUMBER)-1);

	for (i = 0; i < ((THREADS_NUMBER)-1); ++i) {
		thr_mul_array_inplace_by_real_number[i].input_array = array;
		thr_mul_array_inplace_by_real_number[i].index = i*each_thread_elements;
		thr_mul_array_inplace_by_real_number[i].elements_number = each_thread_elements;
		thr_mul_array_inplace_by_real_number[i].factor = coeff;

		my_thr_data_assign (i,  (void *) &thr_mul_array_inplace_by_real_number[i]);
	}

	thr_mul_array_inplace_by_real_number [((THREADS_NUMBER)-1)].input_array = array;
	thr_mul_array_inplace_by_real_number [((THREADS_NUMBER)-1)].index = ((THREADS_NUMBER)-1)*each_thread_elements;
	thr_mul_array_inplace_by_real_number [((THREADS_NUMBER)-1)].elements_number = last_thread_elements;
	thr_mul_array_inplace_by_real_number [((THREADS_NUMBER)-1)].factor = coeff;

	my_thr_data_assign (((THREADS_NUMBER)-1),  (void *) &thr_mul_array_inplace_by_real_number[((THREADS_NUMBER)-1)]);

	my_thr_manager (mul_array_inplace_by_real_number_thr);
#endif /* THREADS_NUMBER */
}

void mul_array_inplace_by_real_number_thr (void *thr_input) {
	unsigned long int i, my_start_index;
	complex_t *temp_ptr1;
	double factor;
	struct mul_array_inplace_by_real_number_thread_input *input;

	input = (struct mul_array_inplace_by_real_number_thread_input*) thr_input;

	factor = input->factor;
	temp_ptr1 = input->input_array + input->index;

	for (i = 0; i < input->elements_number; ++i) {
		(*temp_ptr1)[0] = ((*temp_ptr1)[0])*factor;
		(*temp_ptr1)[1] = ((*temp_ptr1)[1])*factor;

		++temp_ptr1;
	}
}

double L_2_norm (complex_t *input, unsigned long int points_number) {
	unsigned long int i;
	complex_t *temp_ptr1;
	double result=0.0;

#ifndef THREADS_NUMBER
	temp_ptr1 = input;
	for (i = 0; i < points_number; ++i) {
		result += (*temp_ptr1)[0]*(*temp_ptr1)[0] + (*temp_ptr1)[1]*(*temp_ptr1)[1];

		++temp_ptr1;
	}
#else
#ifdef SUPER_VERBOSE
	printf ("L_2_norm with threads.\n");
#endif
	each_thread_elements = (unsigned long int) (0.5 + (1.0*points_number)/(THREADS_NUMBER));
	last_thread_elements = points_number - each_thread_elements*((THREADS_NUMBER)-1);

	for (i = 0; i < ((THREADS_NUMBER)-1); ++i) {
		thr_L_2_norm[i].input_array = input;
		thr_L_2_norm[i].index = i*each_thread_elements;
		thr_L_2_norm[i].elements_number = each_thread_elements;

		my_thr_data_assign (i,  (void *) &thr_L_2_norm[i]);
	}

	thr_L_2_norm [((THREADS_NUMBER)-1)].input_array = input;
	thr_L_2_norm [((THREADS_NUMBER)-1)].index = ((THREADS_NUMBER)-1)*each_thread_elements;
	thr_L_2_norm [((THREADS_NUMBER)-1)].elements_number = last_thread_elements;

	my_thr_data_assign (((THREADS_NUMBER)-1),  (void *) &thr_L_2_norm[((THREADS_NUMBER)-1)]);

	my_thr_manager (L_2_norm_thr);

	for (i = 0; i < (THREADS_NUMBER); ++i) {
		result += thr_L_2_norm [i].sum;
	}
#endif
	return result;
}

void L_2_norm_thr (void * input_thr) {
	struct L_2_norm_thread_input *input;
	complex_t *temp_ptr1;
	unsigned long int i;

	input = (struct L_2_norm_thread_input *) input_thr;

	input->sum = 0.0;
	temp_ptr1 = input->input_array + input->index;
	for (i = 0; i < input->elements_number; ++i) {
		 input->sum += (*temp_ptr1)[0]*(*temp_ptr1)[0] + (*temp_ptr1)[1]*(*temp_ptr1)[1];

		++temp_ptr1;
	}
}

void arrays_transpose (complex_t *input, complex_t *output, unsigned long int row_elements_N, unsigned long int column_elements_M) {
	unsigned long int row_blocks_number, column_blocks_number, total_number_of_blocks, i,j, row_block_skip, column_block_skip;
	complex_t *block_ptr_in, *block_ptr_out;

/*** We suppose that row_elements_N and column_elements_M are multiples of BLOCK_SIDE_SIZE ***/
	row_blocks_number = row_elements_N/(BLOCK_SIDE_SIZE);
	column_blocks_number = column_elements_M/(BLOCK_SIDE_SIZE);

	total_number_of_blocks = row_blocks_number*column_blocks_number;
	row_block_skip = row_elements_N*(BLOCK_SIDE_SIZE);
	column_block_skip = column_elements_M*(BLOCK_SIDE_SIZE);

#ifndef THREADS_NUMBER
	for (i = 0; i < total_number_of_blocks; ++i) {
		unit_transpose (input, output, i, row_elements_N, column_elements_M, row_blocks_number, column_blocks_number, row_block_skip, column_block_skip);
	}
#else
#ifdef SUPER_VERBOSE
	printf ("arrays_transpose with threads.\n");
#endif
	each_thread_elements = (unsigned long int) (0.5 + (1.0*total_number_of_blocks)/(THREADS_NUMBER));
	last_thread_elements = total_number_of_blocks - each_thread_elements*((THREADS_NUMBER)-1);

	for (i = 0; i < ((THREADS_NUMBER)-1); ++i) {
		thr_arrays_transpose [i].input = input;
		thr_arrays_transpose [i].output = output;
		thr_arrays_transpose [i].index = i*each_thread_elements;
		thr_arrays_transpose [i].elements_number = each_thread_elements;
		thr_arrays_transpose [i].row_elements_N = row_elements_N;
		thr_arrays_transpose [i].column_elements_M = column_elements_M;
		thr_arrays_transpose [i].row_blocks_number = row_blocks_number;
		thr_arrays_transpose [i].column_blocks_number = column_blocks_number;
		thr_arrays_transpose [i].row_block_skip = row_block_skip;
		thr_arrays_transpose [i].column_block_skip = column_block_skip;

		my_thr_data_assign (i,  (void *) &thr_arrays_transpose[i]);
	}

	thr_arrays_transpose [((THREADS_NUMBER)-1)].input = input;
	thr_arrays_transpose [((THREADS_NUMBER)-1)].output = output;
	thr_arrays_transpose [((THREADS_NUMBER)-1)].index = ((THREADS_NUMBER)-1)*each_thread_elements;
	thr_arrays_transpose [((THREADS_NUMBER)-1)].elements_number = last_thread_elements;
	thr_arrays_transpose [((THREADS_NUMBER)-1)].row_elements_N = row_elements_N;
	thr_arrays_transpose [((THREADS_NUMBER)-1)].column_elements_M = column_elements_M;
	thr_arrays_transpose [((THREADS_NUMBER)-1)].row_blocks_number = row_blocks_number;
	thr_arrays_transpose [((THREADS_NUMBER)-1)].column_blocks_number = column_blocks_number;
	thr_arrays_transpose [((THREADS_NUMBER)-1)].row_block_skip = row_block_skip;
	thr_arrays_transpose [((THREADS_NUMBER)-1)].column_block_skip = column_block_skip;

	my_thr_data_assign (((THREADS_NUMBER)-1),  (void *) &thr_arrays_transpose[((THREADS_NUMBER)-1)]);

	my_thr_manager (arrays_transpose_thr);
#endif
}

void arrays_transpose_thr (void * input_thr) {
	struct arrays_transpose_thread_input *input;
	unsigned long int i, final_index;

	input = (struct arrays_transpose_thread_input *) input_thr;

	final_index = input->index + input->elements_number;

	for (i = input->index; i < final_index; ++i) {
		unit_transpose (input->input, input->output, i, input->row_elements_N, input->column_elements_M, input->row_blocks_number, input->column_blocks_number, input->row_block_skip, input->column_block_skip);
	}
}

void unit_transpose (complex_t *input, complex_t *output, unsigned long int current_block, unsigned long int row_elements_N, unsigned long int column_elements_M, unsigned long int row_blocks_number, unsigned long int column_blocks_number, unsigned long int row_block_skip, unsigned long int column_block_skip) {
	complex_t *block_input, *block_output, *temp_ptr_in, *temp_ptr_out;
	unsigned long int i, j;

	j = (unsigned long int) (current_block/row_blocks_number);
	i = current_block - (row_blocks_number*j);

	block_input = input + i*(BLOCK_SIDE_SIZE) + j*row_block_skip;
	block_output = output + i*column_block_skip + j*(BLOCK_SIDE_SIZE);

	for (j = 0; j < (BLOCK_SIDE_SIZE); ++j) {
		temp_ptr_in = block_input + row_elements_N*j;
		temp_ptr_out = block_output + j;
		for (i = 0; i < (BLOCK_SIDE_SIZE); ++i) {
			(*temp_ptr_out)[0] = (*temp_ptr_in)[0];
			(*temp_ptr_out)[1] = (*temp_ptr_in)[1];

			++temp_ptr_in;
			temp_ptr_out += column_elements_M;
		}
	}
}

void abs_sqr_cdata (complex_t *input, double *output, unsigned long int points_number) {
	unsigned long int i;
	complex_t *temp_ptr1;
	double *temp_ptr_r;

#ifndef THREADS_NUMBER
	temp_ptr1 = input;
	temp_ptr_r = output;
	for (i = 0; i < points_number; ++i) {
		(*temp_ptr_r) = (*temp_ptr1)[0]*(*temp_ptr1)[0] + (*temp_ptr1)[1]*(*temp_ptr1)[1];

		++temp_ptr1; ++temp_ptr_r;
	}
#else
#ifdef SUPER_VERBOSE
	printf ("abs_sqr_cdata with threads.\n");
#endif
	each_thread_elements = (unsigned long int) (0.5 + (1.0*points_number)/(THREADS_NUMBER));
	last_thread_elements = points_number - each_thread_elements*((THREADS_NUMBER)-1);

	for (i = 0; i < ((THREADS_NUMBER)-1); ++i) {
		thr_abs_sqr_cdata[i].input_array = input;
		thr_abs_sqr_cdata[i].output_array = output;
		thr_abs_sqr_cdata[i].index = i*each_thread_elements;
		thr_abs_sqr_cdata[i].elements_number = each_thread_elements;

		my_thr_data_assign (i,  (void *) &thr_abs_sqr_cdata[i]);
	}

	thr_abs_sqr_cdata [((THREADS_NUMBER)-1)].input_array = input;
	thr_abs_sqr_cdata [((THREADS_NUMBER)-1)].output_array = output;
	thr_abs_sqr_cdata [((THREADS_NUMBER)-1)].index = ((THREADS_NUMBER)-1)*each_thread_elements;
	thr_abs_sqr_cdata [((THREADS_NUMBER)-1)].elements_number = last_thread_elements;

	my_thr_data_assign (((THREADS_NUMBER)-1),  (void *) &thr_abs_sqr_cdata[((THREADS_NUMBER)-1)]);

	my_thr_manager (abs_sqr_cdata_thr);
#endif /* THREADS_NUMBER */
}

void abs_sqr_cdata_thr (void * input_thr) {
	struct abs_sqr_cdata_thread_input *input;
	complex_t *temp_ptr1;
	double *temp_ptr_r;
	unsigned long int i;

	input = (struct abs_sqr_cdata_thread_input *) input_thr;

	temp_ptr1 = input->input_array + input->index;
	temp_ptr_r = input->output_array + input->index;
	for (i = 0; i < input->elements_number; ++i) {
		(*temp_ptr_r) = (*temp_ptr1)[0]*(*temp_ptr1)[0] + (*temp_ptr1)[1]*(*temp_ptr1)[1];

		++temp_ptr1; ++temp_ptr_r;
	}
}

void cexp_I_real_calc (double *input, complex_t *output, unsigned long int points_number) {
	unsigned long int i;
	complex_t *temp_ptr1;
	double *temp_ptr_r;
	complex double exp_i_S;

#ifndef THREADS_NUMBER
	temp_ptr1 = output;
	temp_ptr_r = input;
	for (i = 0; i < points_number; ++i) {
		exp_i_S = cexp(I*(*temp_ptr_r));
		(*temp_ptr1)[0] = creal(exp_i_S);
		(*temp_ptr1)[1] = cimag(exp_i_S);

		++temp_ptr1; ++temp_ptr_r;
	}
#else
#ifdef SUPER_VERBOSE
	printf ("cexp_I_real_calc with threads.\n");
#endif
	each_thread_elements = (unsigned long int) (0.5 + (1.0*points_number)/(THREADS_NUMBER));
	last_thread_elements = points_number - each_thread_elements*((THREADS_NUMBER)-1);

	for (i = 0; i < ((THREADS_NUMBER)-1); ++i) {
		thr_cexp_I_real_calc[i].input_array = input;
		thr_cexp_I_real_calc[i].output_array = output;
		thr_cexp_I_real_calc[i].index = i*each_thread_elements;
		thr_cexp_I_real_calc[i].elements_number = each_thread_elements;

		my_thr_data_assign (i,  (void *) &thr_cexp_I_real_calc[i]);
	}

	thr_cexp_I_real_calc [((THREADS_NUMBER)-1)].input_array = input;
	thr_cexp_I_real_calc [((THREADS_NUMBER)-1)].output_array = output;
	thr_cexp_I_real_calc [((THREADS_NUMBER)-1)].index = ((THREADS_NUMBER)-1)*each_thread_elements;
	thr_cexp_I_real_calc [((THREADS_NUMBER)-1)].elements_number = last_thread_elements;

	my_thr_data_assign (((THREADS_NUMBER)-1),  (void *) &thr_cexp_I_real_calc[((THREADS_NUMBER)-1)]);

	my_thr_manager (cexp_I_real_calc_thr);
#endif /* THREADS_NUMBER */
}

void cexp_I_real_calc_thr (void * input_thr) {
	struct cexp_I_real_calc_thread_input *input;
	complex_t *temp_ptr1;
	double *temp_ptr_r;
	complex double exp_i_S;
	unsigned long int i;

	input = (struct cexp_I_real_calc_thread_input *) input_thr;

	temp_ptr_r = input->input_array + input->index;
	temp_ptr1 = input->output_array + input->index;
	for (i = 0; i < input->elements_number; ++i) {
		exp_i_S = cexp(I*(*temp_ptr_r));
		(*temp_ptr1)[0] = creal(exp_i_S);
		(*temp_ptr1)[1] = cimag(exp_i_S);

		++temp_ptr1; ++temp_ptr_r;
	}
}

/* Some specific operations for our atmosphere code */

void S_k_calc (complex_t *rng_pairs, complex_t *output) {
	complex_t *temp_ptr1, *temp_ptr2;
	double *temp_ptr_r;
	unsigned long int i;
	complex double xi;

#ifndef THREADS_NUMBER
	/* Calculating noise-size S_k */
	temp_ptr1 = rng_pairs;
	temp_ptr2 = output;
	temp_ptr_r = xi_k_mul;
	for (i=0; i < noise_NX_NY_cr; ++i) {
		xi = (*temp_ptr_r)*sqrt(-2.0*log((*temp_ptr1)[0]))*cexp(I*double_M_PI*((*temp_ptr1)[1]));
		(*temp_ptr2)[0] = creal(xi);
		(*temp_ptr2)[1] = cimag(xi);

		++temp_ptr1; ++temp_ptr2; ++temp_ptr_r;
	}
#else
#ifdef SUPER_VERBOSE
	printf ("S_k_calc with threads.\n");
#endif
	unsigned long int each_thread_elements = (unsigned long int) (0.5 + (1.0*noise_NX_NY_cr)/(THREADS_NUMBER));
	unsigned long int last_thread_elements = noise_NX_NY_cr - each_thread_elements*((THREADS_NUMBER)-1);

	for (i = 0; i < ((THREADS_NUMBER)-1); ++i) {
		thr_S_k_calc[i].rng_pairs = rng_pairs;
		thr_S_k_calc[i].output = output;
		thr_S_k_calc[i].index = i*each_thread_elements;
		thr_S_k_calc[i].elements_number = each_thread_elements;

		my_thr_data_assign (i,  (void *) &thr_S_k_calc[i]);
	}

	thr_S_k_calc[((THREADS_NUMBER)-1)].rng_pairs = rng_pairs;
	thr_S_k_calc[((THREADS_NUMBER)-1)].output = output;
	thr_S_k_calc[((THREADS_NUMBER)-1)].index = ((THREADS_NUMBER)-1)*each_thread_elements;
	thr_S_k_calc[((THREADS_NUMBER)-1)].elements_number = last_thread_elements;

	my_thr_data_assign (((THREADS_NUMBER)-1),  (void *) &thr_S_k_calc[((THREADS_NUMBER)-1)]);

	my_thr_manager (S_k_calc_thr);
#endif /* THREADS_NUMBER */
}

void S_k_calc_thr (void *thr_input) {
	complex_t *temp_ptr1, *temp_ptr2;
	double *temp_ptr_r;
	unsigned long int i;
	complex double xi;

	struct S_k_calc_thread_input *input;

	input = (struct S_k_calc_thread_input *) thr_input;

	temp_ptr1 = input->rng_pairs + input->index;
	temp_ptr2 = input->output + input->index;
	temp_ptr_r = xi_k_mul + input->index;
	for (i=0; i < input->elements_number; ++i) {
		xi = (*temp_ptr_r)*sqrt(-2.0*log((*temp_ptr1)[0]))*cexp(I*double_M_PI*((*temp_ptr1)[1]));
		/* xi = Over_sqrt_delta_kx_delta_ky*sqrt(-2.0*log((*temp_ptr1)[0]))*cexp(I*double_M_PI*((*temp_ptr1)[1])); */
		/* xi = sqrt(-2.0*log((*temp_ptr1)[0]))*cexp(I*double_M_PI*((*temp_ptr1)[1])); */
		(*temp_ptr2)[0] = creal(xi);
		(*temp_ptr2)[1] = cimag(xi);

		++temp_ptr1; ++temp_ptr2; ++temp_ptr_r;
	}
}


void rng_pairs_calc (complex_t *output, unsigned long int pairs_number) {
	unsigned long int i;
	complex_t *temp_ptr1;

	gsl_rng_memcpy (my_rng_current_realization, my_rng);
#ifndef THREADS_NUMBER
	temp_ptr1 = output;
	for (i=0; i < pairs_number; ++i) {
		(*temp_ptr1)[0] = gsl_rng_uniform_pos(my_rng);
		(*temp_ptr1)[1] = gsl_rng_uniform_pos(my_rng);

		++temp_ptr1;
	}
#else
#ifdef SUPER_VERBOSE
	printf ("rng_pairs_calc with threads.\n");
#endif
	unsigned long int each_thread_elements = (unsigned long int) (0.5 + (1.0*pairs_number)/(THREADS_NUMBER));
	unsigned long int last_thread_elements = pairs_number - each_thread_elements*((THREADS_NUMBER)-1);

	for (i = 0; i < ((THREADS_NUMBER)-1); ++i) {
		thr_rng_pairs_calc[i].output = output;
		thr_rng_pairs_calc[i].index = i*each_thread_elements;
		thr_rng_pairs_calc[i].elements_number = each_thread_elements;

		my_thr_data_assign (i,  (void *) &thr_rng_pairs_calc[i]);
	}

	thr_rng_pairs_calc[((THREADS_NUMBER)-1)].output = output;
	thr_rng_pairs_calc[((THREADS_NUMBER)-1)].index = ((THREADS_NUMBER)-1)*each_thread_elements;
	thr_rng_pairs_calc[((THREADS_NUMBER)-1)].elements_number = last_thread_elements;

	my_thr_data_assign (((THREADS_NUMBER)-1),  (void *) &thr_rng_pairs_calc[((THREADS_NUMBER)-1)]);

	my_thr_manager (rng_pairs_calc_thr);
#endif /* THREADS_NUMBER */
}

void rng_pairs_calc_thr (void *thr_input) {
	complex_t *temp_ptr1;
	unsigned long int i;

	struct rng_pairs_calc_thread_input *input;

	input = (struct rng_pairs_calc_thread_input *) thr_input;

	temp_ptr1 = input->output + input->index;
	for (i=0; i < input->elements_number; ++i) {
		(*temp_ptr1)[0] = gsl_rng_uniform_pos(input->thread_rng);;
		(*temp_ptr1)[1] = gsl_rng_uniform_pos(input->thread_rng);;

		++temp_ptr1;
	}
}

unsigned long int find_max_index (double *positive_input, unsigned long int points_number) {
	unsigned long int i, max_index=0;
	double max=0, *temp_ptr_r;

	temp_ptr_r = positive_input;
	for (i=0; i<points_number; ++i) {
		if ((*temp_ptr_r)>max) {
			max = (*temp_ptr_r);
			max_index = i;
		}
		++temp_ptr_r;
	}

	return max_index;
}
