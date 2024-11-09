/*** Program for generation of initial conditions ***/

#include "dataio.h"

/*** #include <unistd.h>
#include <fenv.h>
#include <time.h> ***/

#define _VERBOSE
#define _SAVE_FILE_GNU


int main (int argc, char **argv) {
	FILE *f_p;
	fftw_complex *complex_array1, *temp_ptr1;
	double x, y, distance=0.0;
	unsigned long int i, j, ij;

	/*Parsing command line.*/
	if (argc < 3)
	{
		printf ("Using is the following:\n");
		printf("%s conf_file_name output_data_file_name\n", argv[0]);
		exit (1);
	}

	unsigned long int NX=0, NY=0, noise_NX=0, noise_NY=0, realizations_number=0, realizations_save_step=0;
	double LX=0.0, LY=0.0, delta_z=0.0, lambda_0=0.0, n_0=0.0, C_n_sqr=0.0, w_0=0.0, l_0=0.0, L_0=0.0;

    unsigned long int N_BINS;
    double max_intensity_value = 0.0;

	read_conf_file (argv[1], &NX, &NY, &noise_NX, &noise_NY, &LX, &LY, &delta_z, &lambda_0, &n_0, &C_n_sqr, &w_0, &l_0, &L_0, &realizations_number, &realizations_save_step, &max_intensity_value, &N_BINS);

#ifdef VERBOSE
	printf ("NX=%lu	NY=%lu	LX=%.6e	LY=%.6e	w_0=%.6e\n", NX, NY, LX, LY, w_0);
#endif /* VERBOSE */
	double delta_x = LX/NX;
	double delta_y = LY/NY;
	double LX_2 = 0.5*LX;
	double LY_2 = 0.5*LY;
	double Over_sqr_w_0 = 1.0/(w_0*w_0);
	double sum=0.0;

	complex_array1 = (fftw_complex*) fftw_malloc (NX*NY*sizeof(fftw_complex));

	for (i=0; i < NY; ++i) {
		y = delta_y*i - LY_2;
		for (j=0; j < NX; ++j) {
			x = delta_x*j - LX_2;
			ij = NX*i + j;
			complex_array1[ij][0] = exp(-Over_sqr_w_0*(x*x + y*y));
			complex_array1[ij][1] = 0.0;
			sum += (complex_array1[ij][0]);
		}
	}
	printf ("Integral of |psi(r)|=%.15e\n", sum*(LX*LY/(NX*NY)));

#ifdef SAVE_FILE_GNU
	save_GNU_XY_c ("GNU.really_initial_data_XY.cdata", complex_array1, NX, NY, LX, LY, 4);
#endif /* SAVE_FILE_GNU */

	save_data_complex (argv[2], complex_array1, &NX, &NY, &LX, &LY, &distance);

	return 0;
}
