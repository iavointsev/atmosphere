/* Program for simulation of optic pulses propagation in turbulent atmosphere. By Korotkevich A.O., alexkor@math.unm.edu*/
#include <math.h>

#include <fftw3.h>
#include <gsl/gsl_rng.h>

#include <time.h>

#include <complex.h>
#include "utils.h" 
#include "atmosphere.h"

/* NX and NY should be even */
unsigned long int NX=0, NY=0, noise_NX=0, noise_NY=0, realizations_number=0;
double LX=0.0, LY=0.0, delta_z=0.0, lambda_0=0.0, n_0=0.0, C_n_sqr=0.0, w_0=0.0, l_0=0.0, L_0=0.0;

unsigned long int N_BINS = 0;
double max_intensity_value = 0.0;

double delta_k_x, delta_k_y, Over_NX_NY, k_0;
unsigned long int NX_NY, NX_NY_cr, noise_NX_NY_cr, NX_cr;

complex_t *psi, *psi0, *psi_k, *xi, *exp_refr, *S_k, *exp_diffr;
double *abs_k_xy_cr, *Phi_k, *xi_k_mul, *S, *psi_abs_sqr, *intensity_maximums, Over_sqrt_delta_kx_delta_ky;
fftw_plan plan_fwd, plan_bwd, plan_bwd_S_k;

#ifdef SCALES_FILTER
double *scales_filter_mask;
double k_mask_smallest_sqr=0.0, k_mask_largest_sqr=0.0;
#endif /* SCALES_FILTER */

const double double_M_PI = 2.0*(M_PI);

gsl_rng *my_rng, *my_rng_current_realization, *my_rng_best_realization;

unsigned long int each_thread_elements, last_thread_elements;

char file_name[256];

#include "arrays_ops.h"
#include "arrays_ops.c"
#include "arrays_calc.c"
#include "arrays_init.c"


int main (int argc, char** argv) {
	char name[256], name2[256], ext[256], temp_char[256];
	double distance_initial=0.0, distance_current, LX_read=0.0, LY_read=0.0, seconds, max_intensity_overall=0.0, zero_avrg=0.0;
	unsigned long int z_step, step_number=0, NX_read=0, NY_read=0, realizations_save_step=0, realization_current, max_intensity_index, zero_index;
	struct timespec start, finish;
	
	/*Parsing command line.*/
	if (argc < 5)
	{
		printf ("Using is the following:\n");
		printf("%s conf_file_name initial_data_file_name step_number output_data_file_name [saved_rng_state or 0]\n", argv[0]);
		exit (1);
	}

	read_conf_file (argv[1], &NX, &NY, &noise_NX, &noise_NY, &LX, &LY, &delta_z, &lambda_0, &n_0, &C_n_sqr, &w_0, &l_0, &L_0, &realizations_number, &realizations_save_step, &max_intensity_value, &N_BINS);
#ifdef SCALES_FILTER
	k_mask_smallest_sqr = double_M_PI/L_0;
	k_mask_smallest_sqr *= k_mask_smallest_sqr; /* To make squared k */
	k_mask_largest_sqr = double_M_PI/l_0;
	k_mask_largest_sqr *= k_mask_largest_sqr; /* To make squared k */
#endif /* SCALES_FILTER */

	delta_k_x = 2.0*(M_PI)/LX;
	delta_k_y = 2.0*(M_PI)/LY;
	Over_NX_NY = 1.0/(NX*NY);
	k_0 = 2.0*(M_PI)*n_0/lambda_0;

	NX_NY = NX*NY;
	NX_cr = (NX/2+1);
	NX_NY_cr = NX_cr*NY;
	noise_NX_NY_cr = (noise_NX/2+1)*noise_NY;

	zero_index = (NY/2)*NX + (NX/2);

	sscanf(argv[3], "%lu", &step_number);
/*** initialization and calculation of arrays ***/
#ifdef THREADS_NUMBER
	my_thr_pool_init ((THREADS_NUMBER));
#endif /* THREADS_NUMBER */

	arrays_init ();
	plans_init ();
	rng_init (argc, argv);
	arrays_fill ();

	read_data_complex (argv[2], psi, &NX_read, &NY_read, &LX_read, &LY_read, &distance_initial); /* Reading data in XY-space */

	if ((NX_read != NX) || (NY_read != NY) || (LX_read != LX) || (LY_read != LY)) {
		printf ("File %s has parameters incompatible with conf-file %s! Exiting...\n", argv[2], argv[1]);
		exit(1);
	}

	memcpy (psi0, psi, NX_NY*sizeof(complex_t));
#ifdef SAVE_INITIAL_XY
	save_GNU_XY_c ("GNU.initial_data_XY.cdata", psi, (NX), (NY), (LX), (LY), 10);
#endif /* SAVE_INITIAL_XY */

#ifdef SAVE_INITIAL_K
	fftw_execute (plan_fwd); /* Out of place forward Fourier psi->psi_k */
	mul_array_inplace_by_real_number (psi_k, Over_NX_NY, NX_NY); /* Norming Fourier image by 1/(NX*NY) */
	save_GNU_k_c ("GNU.initial_data_k.cdata", psi_k, (NX), (NY), (LX), (LY), 10);
#endif /* SAVE_INITIAL_K */

	sprintf (ext, "%05lux%05lu.noise_%05lux%05lu_%s.", NX, NY, noise_NX, noise_NY, gsl_rng_name(my_rng));
#ifdef THREADS_NUMBER
	sprintf (temp_char, "%lu_threads.data", (unsigned long int)(THREADS_NUMBER));
#else
	sprintf (temp_char, "no_threads.data");
#endif /* THREADS_NUMBER */
	strcat (ext, temp_char);
	strcpy (name, "zero_avrg.");
	strcat (name, ext);
	strcpy (name2, "zero_value.");
	strcat (name2, ext);
	distance_current = distance_initial;
#ifdef TIMING
	clock_gettime(CLOCK_MONOTONIC, &start);
#endif /*TIMING */

#ifdef HISTOGRAM
    if (!(max_intensity_value)) {
        printf("max_intensity parameter is zero!");
        exit(1);
    }
    struct histogram_t *histogram = histogram_new();
    histogram_ctor(histogram, NX, realizations_number, max_intensity_value, N_BINS);
#endif /* HISTOGRAM */

#ifdef CORRELATOR
    struct correlator_t *correlator = correlator_new();
    correlator_ctor(correlator, NX);
#endif /* CORRELATOR */

#ifdef MEAN_INTENSITY
    double *mean_I;
    mean_I = (double*) fftw_malloc (NX_NY*sizeof(double));
    for(size_t r = 0; r < NX_NY; ++r)
    {
		mean_I[r] = 0.0;
    }
#endif /* MEAN_INTENSITY */

	for (realization_current=0; realization_current < realizations_number; ++realization_current) {
		memcpy (psi, psi0, NX_NY*sizeof(complex_t));
/*** calculations begin ***/
		for (z_step = 0; z_step < step_number; ++z_step) {
#ifdef VERBOSE
			printf ("Step number %lu, initial distance z= %.15e\n", z_step+1, distance_current);
#endif /* VERBOSE */

#ifdef DIFR_REFR
			/* We suppose that array psi_k contains data in k-space */
			/* Taking into account diffraction */
			/* Switch to the k-space */
			fftw_execute (plan_fwd);
			mul_arrays_inplace_first_by_second (psi_k, exp_diffr, NX_NY); /* exp_diffr already include factor 1/(NX*NY) */ 
			/* Switch back to XY-space */
			fftw_execute (plan_bwd); /* Now result is in psi-array */
			/* Take into account refraction */
			rng_pairs_calc (exp_refr, noise_NX_NY_cr); /* exp_refr = S_k noise-sized */
			calc_exp_refr (exp_refr); /* exp_refr = S_k noise-sized */
			mul_arrays_inplace_first_by_second (psi, exp_refr, NX_NY);
			/* Full step completed */
#else /* REFR->DIFR */
			/* We suppose that array psi_k contains data in k-space */
			/* Take into account refraction */
			rng_pairs_calc (exp_refr, noise_NX_NY_cr); /* exp_refr = S_k noise-sized */
			calc_exp_refr (exp_refr); /* exp_refr = S_k noise-sized */
			mul_arrays_inplace_first_by_second (psi, exp_refr, NX_NY);
			/* Taking into account diffraction */
			/* Switch to the k-space */
			fftw_execute (plan_fwd); /* Now result is in psi_k-array, multiplied by NX*NY */
			mul_arrays_inplace_first_by_second (psi_k, exp_diffr, NX_NY); /* exp_diffr already include factor 1/(NX*NY) */
			/* Switch back to XY-space */
			fftw_execute (plan_bwd); /* Now result is in psi-array */
			/* Full step completed */
#endif /* DIFR_REFR */
			distance_current += delta_z;
		}
        abs_sqr_cdata (psi, psi_abs_sqr, NX_NY);
		zero_avrg += psi_abs_sqr[zero_index];
		append_two_column_ascii_data (name, (realization_current+1), zero_avrg/(realization_current+1));
		append_two_column_ascii_data (name2, (realization_current+1), psi_abs_sqr[zero_index]);

#ifdef MAX_INTENSITY_STAT
		max_intensity_index = find_max_index (psi_abs_sqr, NX_NY);
		intensity_maximums[realization_current] = psi_abs_sqr[max_intensity_index];
		if (intensity_maximums[realization_current] > max_intensity_overall) {
			max_intensity_overall = intensity_maximums[realization_current];
			gsl_rng_memcpy(my_rng_best_realization, my_rng_current_realization);
		}
        if (!((realization_current+1)%realizations_save_step)) {
            sprintf (temp_char, "intensity_maximums_%06lu_realizations.", (realization_current+1));
            strcat (temp_char, ext);
            save_GNU_real_array_with_index (temp_char, intensity_maximums, (realization_current+1));
            sprintf (temp_char, "extreme_case_rng_%06lu_realizations.", (realization_current+1));
            strcat (temp_char, ext);
            save_rng_state (temp_char, my_rng_best_realization);
            }
#endif // MAX_INTENSITY_STAT
       //
#ifdef HISTOGRAM
        histogram_fill (histogram, psi_abs_sqr, NX);
#endif /* HISTOGRAM */

#ifdef CORRELATOR
    correlator_collect_statistics(correlator, psi, NX);
#endif /* CORRELATOR */
#ifdef MEAN_INTENSITY
    mean_intensity_collect_statistics(mean_I, psi_abs_sqr, NX_NY);
#endif /* MEAN_INTENSITY */
    }

#ifdef HISTOGRAM
    histogram_scale (histogram);
    sprintf (name, "histogram_%06lu_realizations.", realizations_number);
	strcat (name, ext);
    histogram_save_to_file (histogram, name);
#endif /* HISTOGRAM */

#ifdef TIMING
	clock_gettime(CLOCK_MONOTONIC, &finish);
	seconds = (finish.tv_sec - start.tv_sec);
	seconds += (finish.tv_nsec - start.tv_nsec) / (1.0e+9);
	printf ("Time for %lu realization each of %lu steps of atmosphere code is %.15e seconds.\n", realizations_number, step_number, seconds);
#endif /*TIMING */

#ifdef MAX_INTENSITY_STAT
	sprintf (name, "intensity_maximums_%06lu_realizations.", realizations_number);
	strcat (name, ext);
	save_GNU_real_array_with_index (name, intensity_maximums, realizations_number);
	sprintf (name, "extreme_case_rng_%06lu_realizations.", realizations_number);
	strcat (name, ext);
	save_rng_state (name, my_rng_best_realization);
#endif // MAX_INTENSITY_STAT

	save_data_complex (argv[4], psi, NX, NY, LX, LY, distance_current);

#ifdef SAVE_FINAL_K
	save_GNU_k_c ("GNU.final_data_1024_k.cdata", psi_k, NX, NY, LX, LY, 1);
#endif /* SAVE_FINAL_K */

#ifdef SAVE_FINAL_XY
	save_GNU_XY_c ("GNU.final_data_1024_XY.cdata", psi, NX, NY, LX, LY, 1);
#endif /* SAVE_INITIAL_XY */

#ifdef CORRELATOR
    correlator_calculate_acf(correlator, NX_NY, realizations_number);
    correlator_save_GNU(correlator, "GNU.correlator.cdata", NX, NY, LX, LY, 2);
#endif /* CORRELATOR */

#ifdef MEAN_INTENSITY
    mean_intensity_scalarize(mean_I, psi_abs_sqr, NX_NY, realizations_number);
    mean_intensity_save_GNU(mean_I, "GNU.mean_intensity.cdata", NX, NY, LX, LY, 2);
#endif /* MEAN_INTENSITY */

/*** freeing memory ***/
#ifdef THREADS_NUMBER
	my_thr_pool_clear ();
#endif

#ifdef HISTOGRAM
    histogram_dtor(histogram);
    free(histogram); 
#endif /* HISTOGRAM */

#ifdef CORRELATOR
    correlator_dtor(correlator);
    free(correlator); 
#endif /* CORRELATOR */

#ifdef MEAN_INTENSITY
	fftw_free ((void*) mean_I);
#endif /* MEAN_INTENSITY */

	plans_destroy ();
	rng_destroy (argc, argv);
	arrays_free ();

	return 0;
}
