#ifndef mean_intensity_h
#define mean_intensity_h

void mean_intensity_collect_statistics(double *mean_I, double *psi_abs_sqr, unsigned long int N_points);
void mean_intensity_scalarize(double *mean_I, double *psi_abs_sqr, unsigned long int N_points, unsigned long int realizations_number);

void mean_intensity_save_GNU(double *mean_I, char* f_name, unsigned long int NX, unsigned long int NY, double LX, double LY, int step);

#endif // mean_intensity_h
