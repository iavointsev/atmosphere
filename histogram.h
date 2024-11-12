/* Header file for histogram functions */
#ifndef histogram_h
#define histogram_h

#include <stdlib.h>

struct histogram_t;

int histogram_isempty(struct histogram_t *histogram);
struct histogram_t* histogram_new();
void histogram_ctor(struct histogram_t *histogram, unsigned long int NX, unsigned long realizations_number, double max_intensity_value, unsigned long int N_BINS);
void histogram_dtor(struct histogram_t *histogram);

void histogram_fill(struct histogram_t *histogram, double *input, size_t NX);
void histogram_scale(struct histogram_t *histogram);
void histogram_save_to_file(struct histogram_t *histogram, char* f_name);


size_t histogram_get_ind_initial(struct histogram_t *histogram);
size_t histogram_get_ind_final(struct histogram_t *histogram);

#endif // histogram_h

