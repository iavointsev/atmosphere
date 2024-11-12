/*Handlie SIGTERM on parma*/

#include "sighandler.h"

#define ERROR_MESSAGE_SIZE 256


extern unsigned long int realization_number;
extern unsigned long int realization_current;
extern unsigned long int NX;
extern unsigned long int NY;
extern unsigned long int NX_NY;
extern double LX;
extern double LY;

#ifdef HISTOGRAM
extern struct histogram_t * histogram;
#endif /* HISTOGRAM */

#ifdef CORRELATOR
extern struct correlator_t * correlator;
#endif /* CORRELATOR */

void sighandler_ctor(struct sigaction sa)
{
	sa.sa_handler = handle_sigterm;
    sa.sa_flags = 0;
	sigemptyset(&sa.sa_mask);

    if (sigaction(SIGTERM, &sa, NULL) == -1) 
	{
        perror("Signal handler error while setting up\n");
        exit(EXIT_FAILURE);
    }
}

void handle_sigterm(int signum) 
{
	char message[ERROR_MESSAGE_SIZE];
	snprintf(message, ERROR_MESSAGE_SIZE, "SIGTERM received, saving progress...\n \t%lu / %lu iterations done", 
			realization_current, realization_number);
    write(STDOUT_FILENO, message, strlen(message));
	
#ifdef HISTOGRAM
	if(!histogram_isempty(histogram))
	{
		snprintf(message, ERROR_MESSAGE_SIZE, "\tSaving histogram...\n");
    	write(STDOUT_FILENO, message, strlen(message));
		histogram_scale(histogram);
		histogram_save_to_file(histogram, "histogram_dump.data");
	}
#endif /* HISTOGRAM */

#ifdef CORRELATOR
	if(!correlator_isempty(correlator))
	{
		snprintf(message, ERROR_MESSAGE_SIZE, "\tSaving acf...\n"); 
    	write(STDOUT_FILENO, message, strlen(message));
		correlator_calculate_acf(correlator, NX_NY, realization_current);
		correlator_save_GNU(correlator, "GNU.correlator.cdata", NX, NY, LX, LY, 2);
	}
#endif /* CORRELATOR */

	_exit(-1);
}

