/*
This is library for reading initial and saving final data.
*/

#include "dataio.h"

void read_data_params (char* f_name, unsigned long int *N_X, unsigned long int *N_Y, double *L_X, double *L_Y, double *z) {
	unsigned long int elements_number;
	FILE* f_p;

	if ((f_p = fopen(f_name, "rb")) == NULL) {
		printf ("Can't read file %s!\n", f_name);
		exit (1);
	}
#ifdef VERBOSE
	printf ("Reading parameters of file %s...\n", f_name);
#endif
	elements_number = fread (N_X, sizeof (unsigned long int), 1, f_p);
	elements_number = fread (N_Y, sizeof (unsigned long int), 1, f_p);
	elements_number = fread (L_X, sizeof (double), 1, f_p);
	elements_number = fread (L_Y, sizeof (double), 1, f_p);
	elements_number = fread (z, sizeof (double), 1, f_p);

	if (elements_number != 1) {
		printf ("Alarm! Alarm! File %s wasn't read succesfully!\n", f_name);
		exit (1);
	}
}

void read_data_complex (char* f_name, fftw_complex *output, unsigned long int *N_X, unsigned long int *N_Y, double *L_X, double *L_Y, double *z) {
	unsigned long int elements_number;
	FILE* f_p;

	if ((f_p = fopen(f_name, "rb")) == NULL) {
		printf ("Can't read file %s!\n", f_name);
		exit (1);
	}
#ifdef VERBOSE
	printf ("Reading file %s...\n", f_name);
#endif
	elements_number = fread (N_X, sizeof (unsigned long int), 1, f_p);
	elements_number = fread (N_Y, sizeof (unsigned long int), 1, f_p);
	elements_number = fread (L_X, sizeof (double), 1, f_p);
	elements_number = fread (L_Y, sizeof (double), 1, f_p);
	elements_number = fread (z, sizeof (double), 1, f_p);
	elements_number = fread (output, sizeof (fftw_complex), (*N_X)*(*N_Y), f_p);

	fclose (f_p);

	if (elements_number != (*N_X)*(*N_Y)) {
		printf ("Alarm! Alarm! File %s wasn't read succesfully!\n", f_name);
		exit (1);
	}
}

void save_data_complex (char* f_name, fftw_complex *input, unsigned long int *N_X, unsigned long int *N_Y, double *L_X, double *L_Y, double *z) {
	unsigned long int elements_number;
	FILE* f_p;

	if ((f_p = fopen(f_name, "w+b")) == NULL) {
		printf ("Can't open/create file %s!\n", f_name);
		exit (1);
	}
#ifdef VERBOSE
	printf ("Saving file %s\n", f_name);
#endif
	elements_number = fwrite (N_X, sizeof (unsigned long int), 1, f_p);
	elements_number = fwrite (N_Y, sizeof (unsigned long int), 1, f_p);
	elements_number = fwrite (L_X, sizeof (unsigned long int), 1, f_p);
	elements_number = fwrite (L_Y, sizeof (unsigned long int), 1, f_p);
	elements_number = fwrite (z, sizeof (double), 1, f_p);
	elements_number = fwrite (input, sizeof (fftw_complex), (*N_X)*(*N_Y), f_p);

	fclose(f_p);

	if (elements_number != (*N_X)*(*N_Y)) {
		printf ("Alarm! Alarm! File %s write has failed!\n", f_name);
		exit (1);
	}
}

void read_data_cr (char* f_name, fftw_complex *output, unsigned long int *N_X, unsigned long int *N_Y, double *L_X, double *L_Y, double *z) {
	unsigned long int elements_number, NX_cr;
	FILE* f_p;

	if ((f_p = fopen(f_name, "rb")) == NULL) {
		printf ("Can't read file %s!\n", f_name);
		exit (1);
	}
#ifdef VERBOSE
	printf ("Reading file %s...\n", f_name);
#endif
	elements_number = fread (N_X, sizeof (unsigned long int), 1, f_p);
	elements_number = fread (N_Y, sizeof (unsigned long int), 1, f_p);
	elements_number = fread (L_X, sizeof (double), 1, f_p);
	elements_number = fread (L_Y, sizeof (double), 1, f_p);
	elements_number = fread (z, sizeof (double), 1, f_p);
	NX_cr = (*N_X)/2+1;
	elements_number = fread (output, sizeof (fftw_complex), NX_cr*(*N_Y), f_p);

	fclose (f_p);

	if (elements_number != NX_cr*(*N_Y)) {
		printf ("Alarm! Alarm! File %s wasn't read succesfully!\n", f_name);
		exit (1);
	}
}

void save_data_cr (char* f_name, fftw_complex *input, unsigned long int *N_X, unsigned long int *N_Y, double *L_X, double *L_Y, double *z) {
	unsigned long int elements_number, NX_cr;
	FILE* f_p;

	if ((f_p = fopen(f_name, "w+b")) == NULL) {
		printf ("Can't open/create file %s!\n", f_name);
		exit (1);
	}
#ifdef VERBOSE
	printf ("Saving file %s\n", f_name);
#endif
	elements_number = fwrite (N_X, sizeof (unsigned long int), 1, f_p);
	elements_number = fwrite (N_Y, sizeof (unsigned long int), 1, f_p);
	elements_number = fwrite (L_X, sizeof (unsigned long int), 1, f_p);
	elements_number = fwrite (L_Y, sizeof (unsigned long int), 1, f_p);
	elements_number = fwrite (z, sizeof (double), 1, f_p);
	NX_cr = (*N_X)/2+1;
	elements_number = fwrite (input, sizeof (fftw_complex), NX_cr*(*N_Y), f_p);

	fclose(f_p);

	if (elements_number != NX_cr*(*N_Y)) {
		printf ("Alarm! Alarm! File %s write has failed!\n", f_name);
		exit (1);
	}
}

void save_GNU_k_c (char* f_name, fftw_complex *input, unsigned long int NX, unsigned long int NY, double LX, double LY, int step) {
	/* Saving harmonics for complex2complex (no symmetry) transform in GNU format */
        long int i, j, ij, flag = 0;
	double delta_kx, delta_ky, kx, ky;
	FILE* f_p;
	/* 
	 * Function for output in GNUPLOT format 
	 * of the complex plain in k-space.
	 */
	delta_kx = 2.0*(M_PI)/LX;
	delta_ky = 2.0*(M_PI)/LY;
	
	if (f_name == NULL) {
                f_p = stdout;
        } else {

                if ((f_p = fopen(f_name, "w+")) == NULL) {
                        printf ("Can't open/create file %s!\n", f_name);
			exit (1);
                }
        }

        printf ("Saving file %s in GNUPLOT mode.\n", f_name);

	for (i = NY/2; i >= 0; --i) {
		for (j = (NX/2+1); j < NX; ++j) {
			if ((i%step == 0) && (j%step == 0)) {
				flag = 1;
				ij = i * NX + j;
				kx = -delta_kx*(NX - j);
				ky = delta_ky*i;

				fprintf (f_p, "%.15e	%.15e	%.15e	%.15e\n", kx, ky, input[ij][0], input[ij][1]);
			}
                }
		
		for (j = 0; j < (NX/2+1); ++j) {
			if ((i%step == 0) && (j%step == 0)) {
				flag = 1;
				ij = i * NX + j;
				kx = delta_kx*j;
				ky = delta_ky*i;

				fprintf (f_p, "%.15e	%.15e	%.15e	%.15e\n", kx, ky, input[ij][0], input[ij][1]);
			}
                }
		
		if (flag) {
			fprintf (f_p, "\n");
			flag = 0;
		}
	}

	for (i = NY-1; i > (NY/2); --i) {
		for (j = (NX/2+1); j < NX; ++j) {
			if ((i%step == 0) && (j%step == 0)) {
				flag = 1;
				ij = i * NX + j;
				kx = -delta_kx*(NX - j);
				ky = -delta_ky*(NY - i);

				fprintf (f_p, "%.15e	%.15e	%.15e	%.15e\n", kx, ky, input[ij][0], input[ij][1]);
			}
                }
		
		for (j = 0; j < (NX/2+1); ++j) {
			if ((i%step == 0) && (j%step == 0)) {
				flag = 1;
				ij = i * NX + j;
				kx = delta_kx*j;
				ky = -delta_ky*(NY - i);

				fprintf (f_p, "%.15e	%.15e	%.15e	%.15e\n", kx, ky, input[ij][0], input[ij][1]);
			}
                }
		
		if (flag) {
			fprintf (f_p, "\n");
			flag = 0;
		}
	}
	
        if (f_name != NULL) fclose(f_p);
}

void save_GNU_k_cr (char* f_name, fftw_complex *input, unsigned long int NX, unsigned long int NY, double LX, double LY, int step) {
	/* Saving harmonics for complex2real transform, i.e. with Hermitian symmetry, in GNU format */
        long int i, j, ij;
	unsigned long int flag = 0, NX_cr = (NX/2+1);
	double delta_kx, delta_ky, kx, ky;
	FILE* f_p;
	/* 
	 * Function for output in GNUPLOT format 
	 * of the complex plain in k-space.
	 */
	delta_kx = 2.0*(M_PI)/LX;
	delta_ky = 2.0*(M_PI)/LY;
	
	if (f_name == NULL) {
                f_p = stdout;
        } else {

                if ((f_p = fopen(f_name, "w+")) == NULL) {
                        printf ("Can't open/create file %s!\n", f_name);
			exit (1);
                }
        }

        printf ("Saving file %s in GNUPLOT mode.\n", f_name);

	for (i = NY/2; i >= 0; --i) {
		for (j = 0; j < NX_cr; ++j) {
			if ((i%step == 0) && (j%step == 0)) {
				flag = 1;
				ij = i * NX_cr + j;
				kx = delta_kx*j;
				ky = delta_ky*i;

				fprintf (f_p, "%.15e	%.15e	%.15e	%.15e\n", kx, ky, input[ij][0], input[ij][1]);
			}
                }
		if (flag) {
			fprintf (f_p, "\n");
			flag = 0;
		}
	}

	for (i = NY-1; i > (NY/2); --i) {
		for (j = 0; j < NX_cr; ++j) {
			if ((i%step == 0) && (j%step == 0)) {
				flag = 1;
				ij = i * NX_cr + j;
				kx = delta_kx*j;
				ky = -delta_ky*(NY - i);

				fprintf (f_p, "%.15e	%.15e	%.15e	%.15e\n", kx, ky, input[ij][0], input[ij][1]);
			}
                }
		if (flag) {
			fprintf (f_p, "\n");
			flag = 0;
		}
	}
	
        if (f_name != NULL) fclose(f_p);
}

void save_GNU_real_k_cr (char* f_name, double *input, unsigned long int NX, unsigned long int NY, double LX, double LY, int step) {
	/* Saving harmonics for complex2real transform, i.e. with Hermitian symmetry, in GNU format */
        long int i, j, ij;
	unsigned long int flag = 0, NX_cr = (NX/2+1);
	double delta_kx, delta_ky, kx, ky;
	FILE* f_p;
	/* 
	 * Function for output in GNUPLOT format 
	 * of the complex plain in k-space.
	 */
	delta_kx = 2.0*(M_PI)/LX;
	delta_ky = 2.0*(M_PI)/LY;
	
	if (f_name == NULL) {
                f_p = stdout;
        } else {

                if ((f_p = fopen(f_name, "w+")) == NULL) {
                        printf ("Can't open/create file %s!\n", f_name);
			exit (1);
                }
        }

        printf ("Saving file %s in GNUPLOT mode.\n", f_name);

	for (i = NY/2; i >= 0; --i) {
		for (j = 0; j < NX_cr; ++j) {
			if ((i%step == 0) && (j%step == 0)) {
				flag = 1;
				ij = i * NX_cr + j;
				kx = delta_kx*j;
				ky = delta_ky*i;

				fprintf (f_p, "%.15e	%.15e	%.15e	%.15e\n", kx, ky, input[ij], input[ij]);
			}
                }
		if (flag) {
			fprintf (f_p, "\n");
			flag = 0;
		}
	}

	for (i = NY-1; i > (NY/2); --i) {
		for (j = 0; j < NX_cr; ++j) {
			if ((i%step == 0) && (j%step == 0)) {
				flag = 1;
				ij = i * NX_cr + j;
				kx = delta_kx*j;
				ky = -delta_ky*(NY - i);

				fprintf (f_p, "%.15e	%.15e	%.15e	%.15e\n", kx, ky, input[ij], input[ij]);
			}
                }
		if (flag) {
			fprintf (f_p, "\n");
			flag = 0;
		}
	}
	
        if (f_name != NULL) fclose(f_p);
}

void angle_sqr_avrg_GNU_k_cr (char* f_name, fftw_complex *input, unsigned long int NX, unsigned long int NY, double LX, double LY) {
	/* Saving histogram for complex2real transform, i.e. with Hermitian symmetry, in GNU format */
        long int i, j, ij;
	unsigned long int NX_cr = (NX/2+1), *norm_numbers, boxindex;
	double delta_kx, delta_ky, kx, ky, k, *values;
	FILE* f_p;
	/* 
	 * Function for averaging of intensity 
	 * in k-space. We will average over NX_cr values of k (boxes).
	 */
	delta_kx = 2.0*(M_PI)/LX;
	delta_ky = 2.0*(M_PI)/LY;
	values = (double *) malloc (NX_cr*sizeof(double));
	norm_numbers = (unsigned long int *) malloc (NX_cr*sizeof(unsigned long int));
	memset (values, 0, NX_cr*sizeof(double));
	memset (norm_numbers, 0, NX_cr*sizeof(unsigned long int));

	if (f_name == NULL) {
                f_p = stdout;
        } else {

                if ((f_p = fopen(f_name, "w+")) == NULL) {
                        printf ("Can't open/create file %s!\n", f_name);
			exit (1);
                }
        }

        printf ("Saving file %s in GNUPLOT mode.\n", f_name);

	for (i = NY/2; i >= 0; --i) {
		for (j = 0; j < NX_cr; ++j) {
			ij = i * NX_cr + j;
			kx = delta_kx*j;
			ky = delta_ky*i;
			k = sqrt(kx*kx + ky*ky);
			boxindex = truncf(k/delta_kx);
			if (boxindex < NX_cr) {
				values[boxindex] += ((input[ij][0])*(input[ij][0]) + (input[ij][1])*(input[ij][1]));
				norm_numbers[boxindex] += 1;
			}
                }
	}

	for (i = NY-1; i > (NY/2); --i) {
		for (j = 0; j < NX_cr; ++j) {
			ij = i * NX_cr + j;
			kx = delta_kx*j;
			ky = -delta_ky*(NY - i);
			k = sqrt(kx*kx + ky*ky);
			boxindex = truncf(k/delta_kx);
			if (boxindex < NX_cr) {
				values[boxindex-1] += ((input[ij][0])*(input[ij][0]) + (input[ij][1])*(input[ij][1]));
				norm_numbers[boxindex-1] += 1;
			}
                }
	}

	for (i = 0; i < NX_cr; ++i) {
		if (values[i] != 0.0) {
			values[i] /= norm_numbers[i];
		}
		fprintf (f_p, "%.15e	%.15e\n", (i+1)*delta_kx, values[i]);
	}

	free (values); free (norm_numbers);

        if (f_name != NULL) fclose(f_p);
}

void save_GNU_XY_c (char* f_name, fftw_complex *input, unsigned long int NX, unsigned long int NY, double LX, double LY, int step) {
        long int i, j, ij, flag = 0;
	double delta_x, delta_y, x, y, LX_2, LY_2;
	FILE* f_p;
	/* 
	 * Function for output in GNUPLOT format 
	 * of the complex plain in XY-space.
	 */
	delta_x = LX/(NX);
	delta_y = LY/(NY);
	LX_2 = 0.5*(LX);
	LY_2 = 0.5*(LY);
	
	if (f_name == NULL) {
                f_p = stdout;
        } else {

                if ((f_p = fopen(f_name, "w+")) == NULL) {
                        printf ("Can't open/create file %s!\n", f_name);
			exit (1);
                }
        }

        printf ("Saving file %s in GNUPLOT mode.\n", f_name);

	for (i = 0; i < NY; ++i) {
		for (j = 0; j < NX; ++j) {
			if ((i%step == 0) && (j%step == 0)) {
				flag = 1;
				ij = i * NX + j;
				x = delta_x*j - LX_2;
				y = delta_y*i - LY_2;

				fprintf (f_p, "%.15e	%.15e	%.15e	%.15e\n", x, y, input[ij][0], input[ij][1]);
			}
                }
		if (flag) {
			fprintf (f_p, "\n");
			flag = 0;
		}
	}
	
        if (f_name != NULL) fclose(f_p);
}

void save_GNU_XY_r (char* f_name, double *input, unsigned long int NX, unsigned long int NY, double LX, double LY, int step) {
        long int i, j, ij, flag = 0;
	double delta_x, delta_y, x, y, LX_2, LY_2;
	FILE* f_p;
	/* 
	 * Function for output in GNUPLOT format 
	 * of the complex plain in XY-space.
	 */
	delta_x = LX/(NX);
	delta_y = LY/(NY);
	LX_2 = 0.5*(LX);
	LY_2 = 0.5*(LY);
	
	if (f_name == NULL) {
                f_p = stdout;
        } else {

                if ((f_p = fopen(f_name, "w+")) == NULL) {
                        printf ("Can't open/create file %s!\n", f_name);
			exit (1);
                }
        }

        printf ("Saving file %s in GNUPLOT mode.\n", f_name);

	for (i = 0; i < NY; ++i) {
		for (j = 0; j < NX; ++j) {
			if ((i%step == 0) && (j%step == 0)) {
				flag = 1;
				ij = i * NX + j;
				x = delta_x*j - LX_2;
				y = delta_y*i - LY_2;

				fprintf (f_p, "%.15e	%.15e	%.15e\n", x, y, input[ij]);
			}
                }
		if (flag) {
			fprintf (f_p, "\n");
			flag = 0;
		}
	}
	
        if (f_name != NULL) fclose(f_p);
}

void save_GNU_X_line_c (char* f_name, fftw_complex *input, unsigned long int NX, unsigned long int NY, double LX, double LY, double Y) {
        unsigned long int j, ij, Y_index, index_offset;
	double delta_x, LX_2, y;
	FILE* f_p;
	/* 
	 * Function for output in GNUPLOT format 
	 * of the line along X-axis corresponding to Y_index in XY-space.
	 */
	delta_x = LX/(NX);
	LX_2 = 0.5*(LX);
	
	if (f_name == NULL) {
                f_p = stdout;
        } else {

                if ((f_p = fopen(f_name, "w+")) == NULL) {
                        printf ("Can't open/create file %s!\n", f_name);
			exit (1);
                }
        }

        printf ("Saving file %s in GNUPLOT mode.\n", f_name);

	Y_index = (unsigned long int)round((Y/LY +0.5)*NY); /* Y_index = (Y+0.5*LY)*NY/LY; */
	index_offset = NX*Y_index;
	y = (LY*Y_index)/NY - 0.5*LY;

	for (j = 0; j < NX; ++j) {
		ij = index_offset + j;

		fprintf (f_p, "%.15e	%.15e	%.15e	%.15e\n", delta_x*j - LX_2, y, input[ij][0], input[ij][1]);
	}

        if (f_name != NULL) fclose(f_p);
}

void save_point_complex_GNU (char* f_name, fftw_complex input, double index_value1, double index_value2) {
	FILE* f_p;
	if ((f_p = fopen(f_name, "a")) == NULL) {
		printf ("Can't open/create file %s!\n", f_name);
		exit (1);
	}
	fprintf (f_p, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", index_value1, index_value2, sqrt(input[0]*input[0] + input[1]*input[1]), atan2(input[1], input[0]), input[0], input[1]);
	fclose(f_p);
}

void save_point_real_GNU (char* f_name, double input, double index_value1) {
	FILE* f_p;
	if ((f_p = fopen(f_name, "a")) == NULL) {
		printf ("Can't open/create file %s!\n", f_name);
		exit (1);
	}
	fprintf (f_p, "%.15e\t%.15e\n", index_value1, input);
	fclose(f_p);
}

void read_conf_file (char* f_name, unsigned long int *NX, unsigned long int *NY, unsigned long int *noise_NX, unsigned long int *noise_NY, double *LX, double *LY, double *delta_z, double *lambda_0, double *n_0, double *C_n_sqr, double *w_0, double *l_0, double *L_0, unsigned long int *realizations_number, unsigned long int *realizations_save_step, double *max_intensity_value, unsigned long int *N_BINS) {
	FILE* f_p;
	char *buffer, *tok, *remainder;
	unsigned long int buffer_size=256;

	buffer = (char *) malloc ((buffer_size+1)*sizeof(char));
	if ((f_p = fopen(f_name, "r")) == NULL) {
		printf ("Can't open/read file %s!\n", f_name);
		exit (1);
	}

	printf ("Parsing configuration file...\n");

	while(fgets(buffer, buffer_size, f_p)!= NULL){
		remainder=buffer;
		if ((tok = strsep(&remainder, "=")) != NULL) {
			if ((tok != NULL) && (remainder != NULL)) {
				if (!strcmp(tok,"NX")) {
					sscanf(remainder, "%lu", NX);
					printf ("NX=%lu\n", *NX);
				} else if (!strcmp(tok,"NY")) {
					sscanf(remainder, "%lu", NY);
					printf ("NY=%lu\n", *NY);
				} else if (!strcmp(tok,"noise_NX")) {
					sscanf(remainder, "%lu", noise_NX);
					printf ("noise_NX=%lu\n", *noise_NX);
				} else if (!strcmp(tok,"noise_NY")) {
					sscanf(remainder, "%lu", noise_NY);
					printf ("noise_NY=%lu\n", *noise_NY);
				} else if (!strcmp(tok,"LX")) {
					sscanf(remainder, "%le", LX);
					printf ("LX=%.6e\n", *LX);
				} else if (!strcmp(tok,"LY")) {
					sscanf(remainder, "%le", LY);
					printf ("LY=%.6e\n", *LY);
				} else if (!strcmp(tok,"delta_z")) {
					sscanf(remainder, "%le", delta_z);
					printf ("delta_z=%.6e\n", *delta_z);
				} else if (!strcmp(tok,"lambda_0")) {
					sscanf(remainder, "%le", lambda_0);
					printf ("lambda_0=%.6e\n", *lambda_0);
				} else if (!strcmp(tok,"n_0")) {
					sscanf(remainder, "%le", n_0);
					printf ("n_0=%.6e\n", *n_0);
				} else if (!strcmp(tok,"C_n_sqr")) {
					sscanf(remainder, "%le", C_n_sqr);
					printf ("C_n_sqr=%.6e\n", *C_n_sqr);
				} else if (!strcmp(tok,"w_0")) {
					sscanf(remainder, "%le", w_0);
					printf ("w_0=%.6e\n", *w_0);
				} else if (!strcmp(tok,"l_0")) {
					sscanf(remainder, "%le", l_0);
					printf ("l_0=%.6e\n", *l_0);
				} else if (!strcmp(tok,"L_0")) {
					sscanf(remainder, "%le", L_0);
					printf ("L_0=%.6e\n", *L_0);
				} else if (!strcmp(tok,"realizations_number")) {
					sscanf(remainder, "%lu", realizations_number);
					printf ("realizations_number=%lu\n", *realizations_number);
				} else if (!strcmp(tok,"realizations_save_step")) {
					sscanf(remainder, "%lu", realizations_save_step);
					printf ("realizations_save_step=%lu\n", *realizations_save_step);
				} else if (!strcmp(tok,"max_intensity_value")) {
					sscanf(remainder, "%le", max_intensity_value);
					printf ("max_intensity_value=%.6e\n", *max_intensity_value);
				} else if (!strcmp(tok,"N_BINS")) {
					sscanf(remainder, "%lu", N_BINS);
					printf ("N_BINS=%lu\n", *N_BINS);
				} else {
					printf ("Variable \"%s\" is not known.\n", tok);
				}
			} else {
				printf ("%s", tok);
			}
		}
	}

	fflush(stdout);
	if (!((*NX)*(*NY)*(*noise_NX)*(*noise_NY)*(*LX)*(*LY)*(*delta_z)*(*lambda_0)*(*n_0)*(*C_n_sqr)*(*w_0)*(*max_intensity_value))) {
		printf ("Something wrong: one or more variables are zeros! Exiting...\n");
		exit(1);
	}

	free(buffer);
	fclose(f_p);
}

void read_rng_state (char *f_name, gsl_rng *rng) {
	FILE *f_p;

	if ((f_p = fopen(f_name, "rb")) == NULL) {
		printf ("Can't read file %s, initializing with default seed!\n", f_name);
		gsl_rng_set(rng, 0);
	} else if (gsl_rng_fread(f_p, rng)) {
		printf ("Error reading RNG state from file %s, initializing with default seed!\n", f_name);
		gsl_rng_set(rng, 0);
	} else {
		printf ("RNG initialized with state from file %s.\n", f_name);
	}
	fclose(f_p);
}

void save_rng_state (char *f_name, gsl_rng *rng) {
	FILE *f_p;

	if ((f_p = fopen(f_name, "w+")) == NULL) {
		printf ("Can't open file %s for writing, \n", f_name);
		exit (1);
	}
	if (gsl_rng_fwrite(f_p, rng)) {
		printf ("Error writing RNG state to file %s, RNG data is lost!\n", f_name);
	} else {
		printf ("RNG state saved to file %s.\n", f_name);
	}
	fclose(f_p);
}

void append_two_column_ascii_data (char* f_name, unsigned long int index, double value) {
	FILE* f_p;

	if ((f_p = fopen(f_name, "a")) == NULL) {
		printf ("Can't open file %s!\n", f_name);
		exit (1);
	}
	
	fprintf (f_p, "\n");
	fprintf (f_p, "%.06lu\t%.15e", index, value);

	fclose (f_p);
}

void save_GNU_real_array_with_index (char* f_name, double *input, unsigned long int points_number) {
        unsigned long int ij;
	FILE* f_p;
	/* 
	 * Function for output in GNUPLOT format 
	 * real array with array index as the first column
	 */
        if (f_name == NULL) {
                f_p = stdout;
        } else {

                if ((f_p = fopen(f_name, "w+")) == NULL) {
                        printf ("Can't open/create file %s!\n", f_name);
			exit (1);
                }
        }

        printf ("Saving file %s in GNUPLOT mode.\n", f_name);
	
	for (ij = 0; ij < points_number; ++ij) {
		fprintf (f_p, "%lu	%.15e\n", ij, input [ij]);
	}

        if (f_name != NULL) fclose(f_p);
}

void truncate_file (char* f_name) {
	FILE* f_p;
	if ((f_p = fopen(f_name, "w+")) == NULL) {
		printf ("Can't open/create file %s!\n", f_name);
		exit (1);
	}
}
