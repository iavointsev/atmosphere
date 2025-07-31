/*
compare two complex PSI arrays in XY plane.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fftw3.h>
#include <complex.h>
#include "utils.h"

double L_infty_norm (complex_t *input, unsigned long int points_number) {
	double max=0.0, current_abs_sqr;
	unsigned long int i;

	for (i=0; i < points_number; ++i) {
		current_abs_sqr = (input[i][0])*(input[i][0]) + (input[i][1])*(input[i][1]);
		if (current_abs_sqr > max) max = current_abs_sqr;
	}

	return sqrt(max);
}

unsigned long int L_infty_norm_index (complex_t *input, unsigned long int points_number) {
	double max=0.0, current_abs_sqr;
	unsigned long int i, max_index=0;

	for (i=0; i < points_number; ++i) {
		current_abs_sqr = (input[i][0])*(input[i][0]) + (input[i][1])*(input[i][1]);
		if (current_abs_sqr > max) {
			max = current_abs_sqr;
			max_index = i;
		}
	}

	return max_index;
}

int main(int argc, char** argv) {
	ldiv_t	division;
	double z=0.0, LX=0.0, LY=0.0, abs_diff, abs_value;
	complex_t *input1, *input2;
	unsigned long int factor=0, i, j, ij, ij_big, big_NX, big_NY, NX=0, NY=0, step_x=0, step_y=0;
	
	/*Parsing command line.*/
	if (argc < 6)
	{
		printf ("Using is the following (we suppose file1<=file2):\n");
		printf("%s file1_name step_x step_y file2_name diff_file_name\n", argv[0]);
		exit (1);
	}

	read_data_params (argv[1], &NX, &NY, &LX, &LY, &z);
	printf ("Parameters read: NX=%lu, NY=%lu, LX=%.15e, LY=%.15e, z=%.15e\n", NX, NY, LX, LY, z);
	sscanf (argv [2], "%lu", &step_x);
	sscanf (argv [3], "%lu", &step_y);
	big_NX = NX*step_x;
	big_NY = NY*step_y;
	input1 = (complex_t*) malloc (NX*NY*sizeof(complex_t));
	input2 = (complex_t*) malloc (big_NX*big_NY*sizeof(complex_t));
	read_data_complex (argv[1], input1, &NX, &NY, &LX, &LY, &z);
	read_data_complex (argv[4], input2, &big_NX, &big_NY, &LX, &LY, &z);

	for (i=0; i < NY; ++i)
		for (j=0; j < NX; ++j) {
			ij = i*NX + j;
			ij_big = (i*step_y)*big_NX + j*step_x;
			
			input1[ij][0] -= input2[ij_big][0];
			input1[ij][1] -= input2[ij_big][1];
		}
	printf ("Old L_infty_norm(file1-file2)=%.15e\n", L_infty_norm(input1, NX*NY));

	ij = L_infty_norm_index(input1, NX*NY);
	/*
	j = (ij%NX);
	i = ((ij-j)/NX);
	*/
	division = ldiv(ij, NX);
	i = division.quot;
	j = division.rem;

	ij_big = (i*step_y)*big_NX + j*step_x;
	abs_diff = sqrt((input1[ij][0])*(input1[ij][0]) + (input1[ij][1])*(input1[ij][1]));
	abs_value = sqrt((input2[ij_big][0])*(input2[ij_big][0]) + (input2[ij_big][1])*(input2[ij_big][1]));

	printf ("L_infty_norm(file1-file2)=%.15e, relative error at that point =%.15e\n", abs_diff, abs_diff/abs_value);
	save_data_complex (argv[5], input1, &NX, &NY, &LX, &LY, &z);
/**********************************************************************/
	free (input1);
	free (input2);

	exit (0);
}
