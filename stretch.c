/*
Scaling of the c2r Fourier array.
*/
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <fftw3.h>
#include <complex.h>
#include "dataio.h"

int main(int argc, char** argv) {
	double z=0.0, LX=0.0, LY=0.0;
	fftw_complex *initial, *output;
	unsigned long int factor=0, i, j, ij, ij_new, new_NX=0, new_NY=0, NX=0, NY=0, NX_cr, new_NX_cr;
	
	/*Parsing command line.*/
	if (argc < 4)
	{
		printf ("Using is the following:\n");
		printf("%s initial_file_name factor new_file_name\n", argv[0]);
		exit (1);
	}

	read_data_params (argv[1], &NX, &NY, &LX, &LY, &z);
	printf ("Parameters read: NX=%lu, NY=%lu, LX=%.15e, LY=%.15e, z=%.15e\n", NX, NY, LX, LY, z);
	NX_cr = NX/2 + 1;
	initial = (fftw_complex*) malloc (NX_cr*NY*sizeof(fftw_complex));
	read_data_cr (argv[1], initial, &NX, &NY, &LX, &LY, &z);

	sscanf (argv [2], "%lu", &factor);
	
	new_NX = NX*factor;
	new_NY = NY*factor;
	new_NX_cr = new_NX/2 + 1;
	printf ("New parameters: NX=%lu, NY=%lu, LX=%.15e, LY=%.15e, z=%.15e\n", new_NX, new_NY, LX, LY, z);
	output = (fftw_complex*) malloc (new_NX_cr*new_NY*sizeof(fftw_complex));

	memset (output, 0, new_NX_cr*NY*sizeof(fftw_complex));

	for (i=0; i <= NY/2; ++i)
		for (j=0; j < NX_cr; ++j) {
			ij = i*NX_cr + j;
			ij_new = i*new_NX_cr + j;
			
			output[ij_new][0] = initial[ij][0];
			output[ij_new][1] = initial[ij][1];
		}

	for (i=(NY/2+1); i < NY; ++i)
		for (j=0; j < NX_cr; ++j) {
			ij = i*NX_cr + j;
			ij_new = (new_NY - NY + i)*new_NX_cr + j;

			output[ij_new][0] = initial[ij][0];
			output[ij_new][1] = initial[ij][1];
		}

	save_data_cr (argv[3], output, &new_NX, &new_NY, &LX, &LY, &z);
/**********************************************************************/
	free (initial);
	free (output);

	exit (0);
}
