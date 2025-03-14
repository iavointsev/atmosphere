#ifndef correlator_h
#define correlator_h

struct correlator_t;

struct correlator_t* correlator_new(void);
void correlator_ctor(struct correlator_t *correlator, unsigned long int NX);
void correlator_dtor(struct correlator_t *correlator);

void correlator_collect_statistics(struct correlator_t *correlator, fftw_complex *psi, unsigned long int NX);
void correlator_calculate_acf(struct correlator_t *correlator, unsigned long int N_points, unsigned long int realizations_number);
void correlator_save_GNU(struct correlator_t *correlator, char* f_name, unsigned long int NX, unsigned long int NY, double LX, double LY, int step);


#endif // correlator_h
