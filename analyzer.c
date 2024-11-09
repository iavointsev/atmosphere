#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <unistd.h>
#include <fenv.h>

/*** Program for analysis of saved data ***/

#include "dataio.h"

#define _SAVE_PSI_K
#define SAVE_X_LINE

int main (int argc, char** argv) {
	
	char file_name[256], name[256], ext[256];
	fftw_complex *input, *temp_ptr1;
	double LX=0.0, LY=0.0, z=0.0, over_N, delta_x;
	unsigned long int i, NX=0, NY=0, save_grid_step=0.0;
	
	/*Parsing command line.*/
	if (argc < 3)
	{
		printf ("Using is the following:\n");
		printf("%s input_file_name save_grid_step\n", argv[0]);
		exit (1);
	}

	sscanf(argv[2],"%lu", &save_grid_step);
	strcpy(file_name,argv[1]);
	strcpy(name,strtok(file_name,"."));
	strcpy(ext,strtok(NULL,"."));

	read_data_params (argv[1], &NX, &NY, &LX, &LY, &z);
	input = (fftw_complex*) fftw_malloc ((NX)*(NY)*sizeof(fftw_complex));

	if ( strcmp(ext, "cdata") == 0) {
#ifdef VERBOSE
		printf("File extension is %s thus we suppose complex data, XY-plane...\n", ext);
#endif
		read_data_complex (argv[1], input, &NX, &NY, &LX, &LY, &z);

		printf ("Distance = %.15e\n", z);
		strcpy (file_name,"GNU.");
		strcat (file_name, name);
		strcat (file_name, "_XY.");
		strcat (file_name, ext);

		save_GNU_XY_c (file_name, input, NX, NY, LX, LY, save_grid_step);
#ifdef SAVE_X_LINE
		strcpy (file_name,"GNU.");
		strcat (file_name, name);
		strcat (file_name, "_X_line.");
		strcat (file_name, ext);
		save_GNU_X_line_c (file_name, input, NX, NY, LX, LY, 0.0);
#endif /* SAVE_X_LINE */
#ifdef SAVE_PSI_K
		fftw_execute(fftw_plan_dft_2d (NY, NX, input, input, -1, FFTW_ESTIMATE));
		over_N = 1.0/(NX*NY);
		temp_ptr1 = input;
		for (i = 0; i < (NX*NY); ++i) {
			(*temp_ptr1)[0] *= over_N;
			(*temp_ptr1)[1] *= over_N;

			++temp_ptr1;
		}
		strcpy (file_name,"GNU.");
		strcat (file_name, name);
		strcat (file_name, "_k.");
		strcat (file_name, ext);
		save_GNU_k_c (file_name, input, NX, NY, LX, LY, save_grid_step);
#endif /* SAVE_PSI_K */
	} else {
		printf("Unknown file extension %s. Exiting...\n", ext);
	}

	free (input);
	return 0;
}
