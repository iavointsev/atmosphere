#include <stdlib.h>
#include <math.h>

#define FILTER_FACTOR 2

typedef struct {
    size_t ind_initial;
    size_t ind_final;
} filter_t;

filter_t* filter_new () {
    return (filter_t*) malloc(sizeof(filter_t));
}

void filter_ctor (filter_t *filter, unsigned long int NX) {
    size_t tmp = pow(2, FILTER_FACTOR);
    filter->ind_initial = NX * (tmp - 1) / (2 * tmp);
    filter->ind_final   = NX * (tmp + 1) / (2 * tmp);
}

void filter_dtor (filter_t *filter) {
}

size_t filter_get_ind_initial(filter_t *filter) {
    return filter->ind_initial;
}

size_t filter_get_ind_final(filter_t *filter) {
    return filter->ind_final;
}
