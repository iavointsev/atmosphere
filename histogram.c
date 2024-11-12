/*Functions for creating histogram */
#include <stdio.h>
#include <stdlib.h>

#include "filter.h" 

#define N_BINS_STD 256

typedef struct
{
    unsigned long int _N_BINS;
    double *x;
    unsigned long long int *values;
    double *scaled;
    double index_factor; 
    double scale_factor;
    struct filter_t *filter;
	int isempty;
} histogram_t;

int histogram_isempty(histogram_t *histogram)
{
	return histogram->isempty;
}

histogram_t* histogram_new()
{
    return (histogram_t*) malloc(sizeof(histogram_t));
}

void histogram_ctor(histogram_t *histogram, unsigned long int NX, unsigned long realizations_number, double max_intensity_value, unsigned long N_BINS)
{
    unsigned char verify_N_BINS = (N_BINS > 0) && !(N_BINS & (N_BINS - 1));
    if (!verify_N_BINS) 
	{
        printf("Warning! %lu is invalid value for number of bins (must be power of 2). Using %lu as standard value", N_BINS, (unsigned long int)N_BINS_STD);
        N_BINS = N_BINS_STD;
    }

    double delta = max_intensity_value / N_BINS;

    printf("N_BINS = %lu", N_BINS);

    histogram->_N_BINS = N_BINS;
    histogram->isempty = 1;
    histogram->x = calloc(N_BINS, sizeof(double));
    histogram->values = calloc(N_BINS, sizeof(unsigned long int));
    histogram->scaled = calloc(N_BINS, sizeof(double));

    for (size_t i = 0; i < N_BINS; ++i) 
	{
        histogram->x[i] = i * delta; 
        histogram->values[i] = 0; 
        histogram->scaled[i] = 0.0; 
    }

    size_t histogram_number_of_points = NX * NX / (FILTER_FACTOR * FILTER_FACTOR);
    histogram->index_factor = N_BINS * 1.0 / max_intensity_value;
    histogram->scale_factor = 1.0 / (histogram_number_of_points * realizations_number);

    histogram->filter = filter_new();
    filter_ctor(histogram->filter, NX);
}

void histogram_dtor(histogram_t *histogram)
{
    free(histogram->x);
    free(histogram->values);
    free(histogram->scaled);
    filter_dtor(histogram->filter);
    free(histogram->filter);
}


void histogram_fill(histogram_t *histogram, double *input, size_t NX) 
{
    size_t ind_initial, ind_final;
    size_t ij, output_ind;

    ind_initial = filter_get_ind_initial(histogram->filter);
    ind_final   = filter_get_ind_final(histogram->filter);

    for (size_t indy = ind_initial; indy < ind_final; ++indy)
	{
        for (size_t indx = ind_initial; indx < ind_final; ++indx)
		{
            ij = indx + NX * indy;
            output_ind = (size_t) (input[ij] * histogram->index_factor);
            if (output_ind >= histogram->_N_BINS)
                continue;
            ++(histogram->values[output_ind]);
        }
    }

	histogram->isempty = 0;
}


void histogram_scale(histogram_t *histogram) 
{
    for (size_t i = 0; i < histogram->_N_BINS; ++i) 
        histogram->scaled[i] = histogram->scale_factor * histogram->values[i];    
}


void histogram_save_to_file(histogram_t *histogram, char* f_name) 
{
    FILE* f_p;

    if (f_name == NULL)
        f_p = stdout;
    else
	{
        if ((f_p = fopen(f_name, "w+")) == NULL) 
		{
            printf ("Can't open/create file %s!\n", f_name);
			exit (-1);
        }
    }
        
    for (size_t i = 0; i < histogram->_N_BINS; ++i)
        fprintf(f_p, "%.15e\t%.15e\n", histogram->x[i], histogram->scaled[i]);
}

size_t histogram_get_ind_initial(histogram_t *histogram) 
{
    return filter_get_ind_initial(histogram->filter);
}

size_t histogram_get_ind_final(histogram_t *histogram) 
{
    return filter_get_ind_final(histogram->filter);
}
