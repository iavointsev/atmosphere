/* This is library for parallelizing of simple linear cycles */

#ifndef _MY_THR_LIB_H
#define _MY_THR_LIB_H

#include <stdio.h>
#include <stdlib.h>

#include <pthread.h>
#include <signal.h>

#define _DEBUG

struct thread_input {
	pthread_t th;
	pthread_mutex_t run;
	pthread_cond_t wait;
	unsigned long int thr_index;
	void *data_in;
	void (*f_ptr) (void*);
};

void my_thr_pool_init (unsigned long int thread_number);
void my_thr_data_assign (unsigned long int thr_index, void * data_input);
void my_thr_pool_clear (void);

void* elementary_thread (void* thr_input);
void my_thr_manager (void (*thr_function_requested) (void*));

#endif /* _MY_THR_LIB_H */
