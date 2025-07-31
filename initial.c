/*** Program for generation of initial conditions ***/

#include "utils.h"

#define _SAVE_FILE_GNU

typedef enum
{
	PLANE,
	GAUSS
} WaveProfile;

typedef struct {
    const char* conf_file;
    WaveProfile profile;
} CmdArgs;


// TODO: убрать в initial.h? или вообще убрать?
static void generate_plane_wave(complex_t *restrict psi,
								size_t NX, size_t NY,
								double x_start, double y_start,
								double dx, double dy,
								double *psi_integral);


static void generate_gauss_wave(complex_t *restrict psi,
								size_t NX, size_t NY,
								double x_start, double y_start,
								double dx, double dy,
								double over_width_sqr,
								double *psi_integral);


complex_t *generate_wave(size_t NX, size_t NY,
							double LX, double LY, double w_0,
							WaveProfile profile)
{
	const double dx = LX / NX;
	const double dy = LY / NY;
	const double x_start = -LX / 2.0;
	const double y_start = -LY / 2.0;
	const double over_width_sqr = (profile == GAUSS) ? 1.0 / (w_0 * w_0) : 0.0;

	complex_t *psi = (complex_t *)fftw_malloc(NX * NY * sizeof(complex_t));
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


static CmdArgs parse_args(int argc, char **argv) {
    CmdArgs args = {
        .conf_file	= NULL,
        .profile 	= PLANE,
    };

	if (argc < 2) {
        fprintf(stderr, "Usage: %s config_file [OPTIONS]\n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -p, --profile [gauss|plane]  Wave profile (default: plane)\n");
        exit(EXIT_FAILURE);
    }

	args.conf_file = argv[1];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--profile") == 0) {
            if (argc < i + 1) {
                fprintf(stderr, "Missing profile value\n");
                exit(EXIT_FAILURE);
            }
            i++;
            if (strcmp(argv[i], "gauss") == 0) {
                args.profile = GAUSS;
            } else if (strcmp(argv[i], "plane") == 0) {
                args.profile = PLANE;
            } else {
                fprintf(stderr, "Invalid profile: %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    return args;
}

static void generate_inital_fname(buffer, size, config) {
	snprintf(buffer, size, "initial.grid_%lux%lu.noise_%lux%lu.realizations_%lu.cdata",
		config.NX, config.NY, config.noise_NX, config.noise_NY, config.realizations_number);
}

int main(int argc, char **argv)
{
	CmdArgs cmd_args = parse_args(argc, argv);

	ProblemConfig config = read_conf_file(cmd_args.conf_file);
	size_t NX = config.NX, NY = config.NY;
	double LX = config.LX, LY = config.LY, w_0 = config.w_0;

	char inital_fname[128];
	generate_inital_fname(buffer, size, config);

	complex_t *psi = generate_wave(NX, NY, LX, LY, w_0, profile);
	save_data_complex(inital_fname, psi, NX, NY, LX, LY, 0.0);

	fftw_free((void *)psi); // небезопасно, если save_data_complex упадет
	return 0;
}


static void generate_gauss_wave(complex_t *restrict psi,
								size_t NX, size_t NY,
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


static void generate_plane_wave(complex_t *restrict psi,
								size_t NX, size_t NY,
								double x_start, double y_start,
								double dx, double dy,
								double *psi_integral)
{
	for (size_t iy = 0; iy < NY; ++iy)
	{
		const double y = y_start + dy * iy;
		const double tanh_y = 0.5 * (tanh(-(y - 1.2) / 0.1) - tanh(-(y + 1.2) / 0.1));
        
        for (size_t ix = 0; ix < NX; ++ix) {
			const double x = x_start + dx * ix;
			const double tanh_x = 0.5 * (tanh(-(x - 1.2) / 0.1) - tanh(-(x + 1.2) / 0.1));
            const double value = tanh_x * tanh_y;
            
            const size_t idx = iy * NX + ix;
            psi[idx][0] = value;
            psi[idx][1] = 0.0;
            *psi_integral += value;
        }
	}
}
