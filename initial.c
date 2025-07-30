/*** Program for generation of initial conditions ***/

#include "dataio.h"

#define _SAVE_FILE_GNU

typedef enum
{
	PLANE,
	GAUSS
} WaveProfile;

// TODO: убрать в initial.h? или вообще убрать?
static void generate_plane_wave(fftw_complex *restrict psi,
								unsigned long int NX, unsigned long int NY,
								double x_start, double y_start,
								double dx, double dy,
								double *psi_integral);

static void generate_gauss_wave(fftw_complex *restrict psi,
								unsigned long int NX, unsigned long int NY,
								double x_start, double y_start,
								double dx, double dy,
								double over_width_sqr,
								double *psi_integral);

fftw_complex *generate_wave(unsigned long int NX, unsigned long int NY,
							double LX, double LY, double w_0,
							WaveProfile profile)
{
	const double dx = LX / NX;
	const double dy = LY / NY;
	const double x_start = -LX / 2.0;
	const double y_start = -LY / 2.0;
	const double over_width_sqr = (profile == GAUSS) ? 1.0 / (w_0 * w_0) : 0.0;

	fftw_complex *psi = (fftw_complex *)fftw_malloc(NX * NY * sizeof(fftw_complex));
	if (!psi)
	{
		fprintf(stderr, "Allocation (%lu x %lu) array failed\n", NX, NY);
		exit(EXIT_FAILURE);
	}

	double psi_integral = 0.0;

	switch (profile)
	{
	case GAUSS:
		generate_gauss_wave(psi, NX, NY, x_start, y_start, dx, dy,
							over_width_sqr, &psi_integral);
		break;
	case PLANE:
		generate_plane_wave(psi, NX, NY, x_start, y_start, dx, dy,
							&psi_integral);
		break;
	default:
		fprintf(stderr, "Unknown wave profile.\n");
		fftw_free(psi);
		exit(EXIT_FAILURE);
	}

	printf("Estimated integral of |psi(r)| = %.15e\n", psi_integral * dx * dy);
	return psi;
}

static void generate_gauss_wave(fftw_complex *restrict psi,
								unsigned long int NX, unsigned long int NY,
								double x_start, double y_start,
								double dx, double dy,
								double over_width_sqr,
								double *psi_integral)
{
	for (size_t iy = 0; iy < NY; ++iy)
	{
		const double y = y_start + dy * iy;
		const double y2 = y * y;
		for (size_t ix = 0; ix < NX; ++ix)
		{
			const double x = x_start + dx * ix;
			const double value = exp(-over_width_sqr * (x * x + y2));

			const size_t idx = iy * NX + ix;
			psi[idx][0] = value;
			psi[idx][1] = 0.0;
			*psi_integral += value;
		}
	}
}

static void generate_plane_wave(fftw_complex *restrict psi,
								unsigned long int NX, unsigned long int NY,
								double x_start, double y_start,
								double dx, double dy,
								double *psi_integral)
{
	for (size_t iy = 0; iy < NY; ++iy)
	{
		const double y = y_start + dy * iy;
		const double tanh_y = 0.5 * (tanh(-(y - 1.2) / 0.1) - tanh(-(y + 1.2) / 0.1);
        
        for (size_t ix = 0; ix < NX; ++ix) {
			const double x = x_start + dx * ix;
			const double tanh_x = 0.5 * (tanh(-(x - 1.2) / 0.1) - tanh(-(x + 1.2) / 0.1);
            const double value = tanh_x * tanh_y;
            
            const size_t idx = iy * NX + ix;
            psi[idx][0] = value;
            psi[idx][1] = 0.0;
            *psi_integral += value;
        }
	}
}

int main(int argc, char **argv)
{
	FILE *f_p;
	fftw_complex *psi;
	/*Parsing command line.*/
	if (argc < 3)
	{
		fprintf(stderr, "Using is the following:\n");
		fprintf(stderr, "\t%s conf_file_name output_data_file_name [gauss | plane]\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	ProblemConfig config = read_conf_file(argv[1]);
	WaveProfile profile = PLANE;
	if (argc >= 4)
	{
		if (!strcmp(argv[3], "gauss"))
		{
			profile = GAUSS;
		}
		else if (!strcmp(argv[3], "plane"))
		{
			profile = PLANE;
		}
		else
		{
			fprintf(stderr, "Unknown profile type %s\n", argv[3]);
			exit(EXIT_FAILURE);
		}
	}

	unsigned long int NX = config.NX,
					  NY = config.NY;

	double LX = config.LX,
		   LY = config.LY,
		   w_0 = config.w_0;

	double distance = 0.0; // TODO: избавиться
	fftw_complex *psi = generate_wave(NX, NY, LX, LY, w_0, profile);

#ifdef SAVE_FILE_GNU
	save_GNU_XY_c("GNU.reallY_initial_data_XY.cdata", psi, NX, NY, LX, LY, 4);
#endif /* SAVE_FILE_GNU */

	save_data_complex(argv[2], psi, &NX, &NY, &LX, &LY, &distance);

	fftw_free((void *)psi);
	return 0;
}
