#define FILTER_FACTOR 2

struct filter_t;

struct filter_t* filter_new ();
void filter_ctor (struct filter_t *filter, unsigned long int NX);
void filter_dtor (struct filter_t *histogram);

size_t filter_get_ind_initial(struct filter_t *filter);
size_t filter_get_ind_final(struct filter_t *filter);

