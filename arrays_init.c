/* This file contains special for our programme operations with arrays */

void arrays_init (void) {
#ifdef THREADS_NUMBER
	thr_mul_array_third_equal_first_by_second = (struct mul_array_third_equal_first_by_second_thread_input *) malloc ((THREADS_NUMBER)*sizeof(struct mul_array_third_equal_first_by_second_thread_input));
	thr_array_mul_arrays_inplace_first_by_second = (struct mul_arrays_inplace_first_by_second_thread_input *) malloc ((THREADS_NUMBER)*sizeof(struct mul_arrays_inplace_first_by_second_thread_input));
	thr_array_mul_arrays_inplace_first_by_real_second = (struct mul_arrays_inplace_first_by_real_second_thread_input *) malloc ((THREADS_NUMBER)*sizeof(struct mul_arrays_inplace_first_by_real_second_thread_input));
	thr_mul_array_inplace_by_real_number = (struct mul_array_inplace_by_real_number_thread_input *) malloc ((THREADS_NUMBER)*sizeof(struct mul_array_inplace_by_real_number_thread_input));
	thr_L_2_norm = (struct L_2_norm_thread_input *) malloc ((THREADS_NUMBER)*sizeof(struct L_2_norm_thread_input));
	thr_arrays_transpose = (struct arrays_transpose_thread_input *) malloc ((THREADS_NUMBER)*sizeof(struct arrays_transpose_thread_input));
	thr_abs_sqr_cdata = (struct abs_sqr_cdata_thread_input *) malloc ((THREADS_NUMBER)*sizeof(struct abs_sqr_cdata_thread_input));
	thr_cexp_I_real_calc = (struct cexp_I_real_calc_thread_input *) malloc ((THREADS_NUMBER)*sizeof(struct cexp_I_real_calc_thread_input));

	thr_S_k_calc = (struct S_k_calc_thread_input *) malloc ((THREADS_NUMBER)*sizeof(struct S_k_calc_thread_input));
	thr_rng_pairs_calc = (struct rng_pairs_calc_thread_input *) malloc ((THREADS_NUMBER)*sizeof(struct rng_pairs_calc_thread_input));
#endif /* THREADS_NUMBER */
	psi = (complex_t*) fftw_malloc (NX_NY*sizeof(complex_t));
	psi0 = (complex_t*) fftw_malloc (NX_NY*sizeof(complex_t));
	psi_k = (complex_t*) fftw_malloc (NX_NY*sizeof(complex_t));
/*	xi = (complex_t*) fftw_malloc (NX_NY*sizeof(complex_t)); */
	S_k = (complex_t*) fftw_malloc (NX_NY_cr*sizeof(complex_t));
	exp_refr = (complex_t*) fftw_malloc (NX_NY*sizeof(complex_t));
	exp_diffr = (complex_t*) fftw_malloc (NX_NY*sizeof(complex_t));

	xi_k_mul = (double*) fftw_malloc (noise_NX_NY_cr*sizeof(double));
	S = (double*) fftw_malloc (NX_NY*sizeof(double));
	psi_abs_sqr = (double*) fftw_malloc (NX_NY*sizeof(double));
#ifdef MAX_INTENSITY_STAT
	intensity_maximums = (double*) fftw_malloc (realizations_number*sizeof(double));
#endif // MAX_INTENSITY_STAT
#ifdef SCALES_FILTER
	scales_filter_mask = (double*) fftw_malloc (noise_NX_NY_cr*sizeof(double));
	memset (scales_filter_mask, 0, noise_NX_NY_cr*sizeof(double));
#endif /* SCALES_FILTER */
}

void arrays_free (void) {
#ifdef THREADS_NUMBER
	free(thr_mul_array_third_equal_first_by_second);
	free(thr_array_mul_arrays_inplace_first_by_second);
	free(thr_mul_array_inplace_by_real_number);
	free(thr_L_2_norm);
	free(thr_arrays_transpose);
	free(thr_abs_sqr_cdata);
	free(thr_cexp_I_real_calc);

	free(thr_S_k_calc);
	free(thr_rng_pairs_calc);
#endif /* THREADS_NUMBER */
	fftw_free ((void*) psi);
	fftw_free ((void*) psi0);
	fftw_free ((void*) psi_k);
/*	fftw_free ((void*) xi); */
	fftw_free ((void*) S_k);
	fftw_free ((void*) exp_refr);
	fftw_free ((void*) exp_diffr);

	fftw_free ((void*) xi_k_mul);
	fftw_free ((void*) S);
	fftw_free ((void*) psi_abs_sqr);
#ifdef MAX_INTENSITY_STAT
	fftw_free ((void*) intensity_maximums);
#endif // MAX_INTENSITY_STAT
#ifdef SCALES_FILTER
	fftw_free ((void*) scales_filter_mask);
#endif /* SCALES_FILTER */
}

void plans_init (void) {
	FILE *f_p;
	char f_name[32];

#ifndef THREADS_NUMBER
	strcpy (f_name, "fftw.wisdom");
#else
	strcpy (f_name, "fftw.threads.wisdom");
	fftw_init_threads();
	fftw_plan_with_nthreads((THREADS_NUMBER));
#endif /* THREADS_NUMBER */
	if ((f_p = fopen(f_name, "rb")) == NULL) {
		printf ("Can't read file %s!\n", f_name);
		printf ("Trying to calculate plans without wisdom...\n");
	} else {
		if (fftw_import_wisdom_from_file (f_p)) {
			printf ("Wisdom successfully imported.\n");
		} else {
			printf ("File %s with wisdom is corrupted.\n", f_name);
			printf ("Trying to calculate plans without wisdom...\n");
		}
		fclose (f_p);
	}
/*
	plan_fwd = fftw_plan_dft_2d ((NY), (NX), psi, psi_k, -1, FFTW_EXHAUSTIVE);
	plan_bwd = fftw_plan_dft_2d ((NY), (NX), psi_k, psi, +1, FFTW_EXHAUSTIVE);
	plan_bwd_S_k = fftw_plan_dft_c2r_2d ((NY), (NX), S_k, S, FFTW_EXHAUSTIVE);
*/
	plan_fwd = fftw_plan_dft_2d ((NY), (NX), psi, psi_k, -1, FFTW_PATIENT);
	plan_bwd = fftw_plan_dft_2d ((NY), (NX), psi_k, psi, +1, FFTW_PATIENT);
	plan_bwd_S_k = fftw_plan_dft_c2r_2d ((NY), (NX), S_k, S, FFTW_PATIENT);

	if ((f_p = fopen(f_name, "w+")) == NULL) {
		printf ("Can't open file %s for writing!\n", f_name);
		printf ("Wisdom has been lost!\n");
	} else {
		printf ("Exporting wisdom to file...\n");
		fftw_export_wisdom_to_file (f_p);
		fclose (f_p);
	}
}

void plans_destroy (void) {
	fftw_destroy_plan(plan_fwd);
	fftw_destroy_plan(plan_bwd);
	fftw_destroy_plan(plan_bwd_S_k);
}

void rng_init (int argc, char** argv) {
	FILE* f_p;
	unsigned long int i=0;

	/* Random number generator initialization */
	my_rng = gsl_rng_alloc (MY_RNG_TYPE);
	my_rng_current_realization = gsl_rng_alloc (MY_RNG_TYPE);
	my_rng_best_realization = gsl_rng_alloc (MY_RNG_TYPE);
	if (argc > 5) {
		if (strncmp(argv[5], "0", 256)) {
			read_rng_state (argv[5], my_rng);
		} else {
			gsl_rng_set(my_rng, 0);
		}
	} else {
#ifdef RNG_INIT_BY_TIME
		gsl_rng_set(my_rng, time(NULL));
#else /*RNG_INIT_BY_TIME*/
		gsl_rng_set(my_rng, 0);
#endif /*RNG_INIT_BY_TIME*/
	}
	printf ("RNG type used: %s\n", gsl_rng_name(my_rng));
	/******************************************/
#ifdef THREADS_NUMBER
#ifdef VERBOSE
	printf ("Initializing RNGs for all threads...\n");
#endif /* VERBOSE */
	for (i=0; i < (THREADS_NUMBER); ++i) {
		thr_rng_pairs_calc[i].thread_rng = gsl_rng_alloc (MY_RNG_TYPE);
		gsl_rng_set(thr_rng_pairs_calc[i].thread_rng, gsl_rng_get (my_rng)); 
	}
#endif /* THREADS_NUMBER */
}

void rng_destroy (int argc, char** argv) {
	FILE* f_p;
	unsigned long int i=0;

	/* Saving random number generator state */
	if (argc > 5) {
		if (strncmp(argv[5], "0", 256)) {
			save_rng_state (argv[5], my_rng);
		}
	}
	gsl_rng_free(my_rng);
	gsl_rng_free(my_rng_current_realization);
	gsl_rng_free(my_rng_best_realization);
	/******************************************/
#ifdef THREADS_NUMBER
#ifdef VERBOSE
	printf ("Freeing RNGs for all threads...\n");
#endif /* VERBOSE */
	for (i=0; i < (THREADS_NUMBER); ++i) {
		gsl_rng_free(thr_rng_pairs_calc[i].thread_rng); 
	}
#endif /* THREADS_NUMBER */
}

void arrays_fill (void) {
	const double power_law = -11.0/3.0;
	const double double_M_PI_delta_z = 2.0*(M_PI)*(delta_z);
	const double delta_k_x_delta_ky = delta_k_x*delta_k_y;
	unsigned long int i, j, ij, NX_2 = (NX/2), NY_2 = (NY/2), noise_NY_2 = (noise_NY/2), noise_NX_cr = (noise_NX/2 + 1);
	double k_x, k_y, abs_k_xy, sqr_k_xy, Phi_k, diffr_phase;

	Over_sqrt_delta_kx_delta_ky =(sqrt(LX*LY)/(2.0*(M_PI)))*(M_SQRT1_2); /* (1/sqrt(2))*1.0/sqrt(delta_k_x*delta_k_y), 1/sqrt(2) comes from requirement for complex Gaussian process x=\xi + I*\zeta, \xi and \zeta are Gaussian processes, <|x|^2> = <\xi^2>+<\zeta^2>=2 */
	/* Creating array to multiply on xi_k, noise-size */
	/* positive k_y's */
	for (i = 0; i < noise_NY_2; ++i) {
		k_y = delta_k_y*i;
		for (j = 0; j < noise_NX_cr; ++j) {
			k_x = delta_k_x*j;
			ij = i*noise_NX_cr + j;
			abs_k_xy = sqrt(k_x*k_x + k_y*k_y);
			if (abs_k_xy!=0) {
				Phi_k = (3.300539063636137e-2)*(C_n_sqr)*pow(abs_k_xy, power_law);
			} else {
				Phi_k = 0.0;
			}
			/* xi_k_mul[ij] =  Over_sqrt_delta_kx_delta_ky*k_0*sqrt(double_M_PI_delta_z*Phi_k)*pow(-1.0, 1*(i+j))*delta_k_x_delta_ky; */
			/* Last two factors come to take into account inverse Fourier transform for:
			a) -k_max/2 to +k_max/2, instead of from 0 to k_max;
			b) integration with real delta_k_x and delta_k_y. */
			xi_k_mul[ij] =  Over_sqrt_delta_kx_delta_ky*k_0*sqrt(double_M_PI_delta_z*Phi_k)*delta_k_x_delta_ky;
		}
	}
	/* negative k_y's */
	for (i = noise_NY_2; i < (noise_NY); ++i) {
		k_y = -delta_k_y*((noise_NY) - i);
		for (j = 0; j < noise_NX_cr; ++j) {
			k_x = delta_k_x*j;
			ij = i*noise_NX_cr + j;
			abs_k_xy = sqrt(k_x*k_x + k_y*k_y);
			if (abs_k_xy!=0) {
				Phi_k = (3.300539063636137e-2)*(C_n_sqr)*pow(abs_k_xy, power_law);
			} else {
				Phi_k = 0.0;
			}
			/* xi_k_mul[ij] = Over_sqrt_delta_kx_delta_ky*k_0*sqrt(double_M_PI_delta_z*Phi_k)*pow(-1.0, 1*(((noise_NY) - i) + j))*delta_k_x_delta_ky; */
			/* Last two factors come to take into account inverse Fourier transform for:
			a) -k_max/2 to +k_max/2, instead of from 0 to k_max;
			b) integration with real delta_k_x and delta_k_y. */
			xi_k_mul[ij] = Over_sqrt_delta_kx_delta_ky*k_0*sqrt(double_M_PI_delta_z*Phi_k)*delta_k_x_delta_ky;
		}
	}

	/* Creating array of exponents to take diffraction into account, include 1/(NX*NY) as well. */
	/* positive k_y's */
	for (i = 0; i < NY_2; ++i) {
		/* positive k_x's and positive k_y's */
		k_y = delta_k_y*i;
		for (j = 0; j < NX_2; ++j) {
			k_x = delta_k_x*j;
			ij = i*NX + j;
			sqr_k_xy = k_x*k_x + k_y*k_y;
			diffr_phase = -0.5*sqr_k_xy*(delta_z)/(k_0);
			exp_diffr[ij][0] = cos(diffr_phase)*Over_NX_NY;
			exp_diffr[ij][1] = sin(diffr_phase)*Over_NX_NY;
		}
		/* negative k_x's and positive k_y's */
		for (j = NX_2; j < (NX); ++j) {
			k_x = -delta_k_x*((NX)-j);
			ij = i*NX + j;
			sqr_k_xy = k_x*k_x + k_y*k_y;
			diffr_phase = -0.5*sqr_k_xy*(delta_z)/(k_0);
			exp_diffr[ij][0] = cos(diffr_phase)*Over_NX_NY;
			exp_diffr[ij][1] = sin(diffr_phase)*Over_NX_NY;
		}
	}
	/* negative k_y's */
	for (i = NY_2; i < (NY); ++i) {
		/* positive k_x's and negative k_y's */
		k_y = -delta_k_y*((NY) - i);
		for (j = 0; j < NX_2; ++j) {
			k_x = delta_k_x*j;
			ij = i*NX + j;
			sqr_k_xy = k_x*k_x + k_y*k_y;
			diffr_phase = -0.5*sqr_k_xy*(delta_z)/(k_0);
			exp_diffr[ij][0] = cos(diffr_phase)*Over_NX_NY;
			exp_diffr[ij][1] = sin(diffr_phase)*Over_NX_NY;
		}
		/* negative k_x's and negative k_y's */
		for (j = NX_2; j < (NX); ++j) {
			k_x = -delta_k_x*((NX)-j);
			ij = i*NX + j;
			sqr_k_xy = k_x*k_x + k_y*k_y;
			diffr_phase = -0.5*sqr_k_xy*(delta_z)/(k_0);
			exp_diffr[ij][0] = cos(diffr_phase)*Over_NX_NY;
			exp_diffr[ij][1] = sin(diffr_phase)*Over_NX_NY;
		}
	}
#ifdef SCALES_FILTER
	double k_sqr;
	/* Creating of a mask array which will leave nonzero only harmonics of the noise with k_mask_smallest_sqr < k_sqr < k_mask_largest_sqr */
	/* positive k_y's */
	for (i = 0; i < noise_NY_2; ++i) {
		k_y = delta_k_y*i;
		for (j = 0; j < noise_NX_cr; ++j) {
			k_x = delta_k_x*j;
			ij = i*noise_NX_cr + j;
			k_sqr = k_x*k_x + k_y*k_y;
			if ((k_sqr > k_mask_smallest_sqr) && (k_sqr < k_mask_largest_sqr)) {
				scales_filter_mask[ij] = 1.0;
			}
		}
	}
	/* negative k_y's */
	for (i = noise_NY_2; i < (noise_NY); ++i) {
		k_y = -delta_k_y*((noise_NY) - i);
		for (j = 0; j < noise_NX_cr; ++j) {
			k_x = delta_k_x*j;
			ij = i*noise_NX_cr + j;
			k_sqr = k_x*k_x + k_y*k_y;
			if ((k_sqr > k_mask_smallest_sqr) && (k_sqr < k_mask_largest_sqr)) {
				scales_filter_mask[ij] = 1.0;
			}
		}
	}
	/* save_GNU_real_k_cr ("GNU_filter_mask.data", scales_filter_mask, noise_NX, noise_NY, LX, LY, 1); */
#endif /* SCALES_FILTER */
}
