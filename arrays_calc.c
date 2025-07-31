
void calc_exp_refr (complex_t *rng_pairs_io) {
	const unsigned long int noise_NY_2 = noise_NY/2;
	const unsigned long int noise_NX_cr = (noise_NX/2 +1);
	double u1, u2;
	complex double xi;
	complex double exp_i_S;
	unsigned long int i, j, ij, ij_negative, ij_new;

	/* instead of S_k using exp_refr as a scratch array */;
	S_k_calc (rng_pairs_io, exp_refr);

	/* Making sure that the Y-line of noise-sized S_k has Hermitian symmetry */
	for (i=1; i < noise_NY_2; ++i) {
		ij = i*noise_NX_cr;
		ij_negative = (noise_NY - i)*noise_NX_cr;
		exp_refr[ij_negative][0] = exp_refr[ij][0];
		exp_refr[ij_negative][1] = -exp_refr[ij][1];
	}
	/* These four corner points due to Hermitian symmetry have to be real */
	exp_refr[0][0] *= (M_SQRT2); /* correcting correlator */
	exp_refr[0][1] = 0.0; /* k_x=0, k_y=0 */
	exp_refr[noise_NX_cr-1][0] *= (M_SQRT2); /* correcting correlator */
	exp_refr[noise_NX_cr-1][1] = 0.0; /* k_x=NX/2, k_y=0 */
	exp_refr[noise_NY_2*noise_NX_cr][0] *= (M_SQRT2); /* correcting correlator */
	exp_refr[noise_NY_2*noise_NX_cr][1] = 0.0; /* k_x=0, k_y=NY/2 */
	exp_refr[noise_NY_2*noise_NX_cr + noise_NX_cr-1][0] *= (M_SQRT2); /* correcting correlator */
	exp_refr[noise_NY_2*noise_NX_cr + noise_NX_cr-1][1] = 0.0; /* k_x=NX/2, k_y=NY/2 */

#ifdef SCALES_FILTER
	/* Applying scales filter to noise */
	/* save_GNU_k_cr ("GNU_unfiltered_noise_spectrum.cdata", exp_refr, noise_NX, noise_NY, LX, LY, 1);
	angle_sqr_avrg_GNU_k_cr ("GNU_averaged_unfiltered_noise_spectrum.cdata", exp_refr, noise_NX, noise_NY, LX, LY); */
	mul_arrays_inplace_first_by_real_second (exp_refr, scales_filter_mask, noise_NX_NY_cr);
	/* save_GNU_k_cr ("GNU_filtered_noise_spectrum.cdata", exp_refr, noise_NX, noise_NY, LX, LY, 1);
	angle_sqr_avrg_GNU_k_cr ("GNU_averaged_filtered_noise_spectrum.cdata", exp_refr, noise_NX, noise_NY, LX, LY); */
#endif /* SCALES_FILTER */

	memset (S_k, 0, NX_cr*NY*sizeof(complex_t));

	for (i=0; i <= noise_NY_2; ++i)
		for (j=0; j < noise_NX_cr; ++j) {
			ij = i*noise_NX_cr + j;
			ij_new = i*NX_cr + j;
			
			S_k[ij_new][0] = exp_refr[ij][0];
			S_k[ij_new][1] = exp_refr[ij][1];
		}

	for (i=(noise_NY_2+1); i < noise_NY; ++i)
		for (j=0; j < noise_NX_cr; ++j) {
			ij = i*noise_NX_cr + j;
			ij_new = (NY - noise_NY + i)*NX_cr + j;

			S_k[ij_new][0] = exp_refr[ij][0];
			S_k[ij_new][1] = exp_refr[ij][1];
		}
	/* Save the S_k noise spectrum *
	angle_sqr_avrg_GNU_k_cr ("GNU_averaged_filtered_spectrum.cdata", S_k, NX, NY, LX, LY);
	*******************************/
	/* Switching back to XY-space */
	fftw_execute (plan_bwd_S_k);
	
	/* save_GNU_XY_r ("GNU.S_XY_with_sign_change2.rdata", S, NX, NY, LX, LY, 1);
	exit(1); */
	
	/* Computing exp_refr = exp(I*S) */
	/*
	for (i=0; i < NX_NY; ++i) {
		exp_i_S = cexp(I*S[i]);
		exp_refr[i][0] = creal(exp_i_S);
		exp_refr[i][1] = cimag(exp_i_S);
	}
	*/
	cexp_I_real_calc (S, exp_refr, NX_NY);
}
