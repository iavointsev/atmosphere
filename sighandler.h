/*Handlie SIGTERM on parma*/

#ifndef sighandler_h
#define sighandler_h

#include <signal.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include "atmosphere.h"
#include "histogram.h"
#include "correlator.h"

void sighandler_ctor(struct sigaction sa);
void handle_sigterm(int signum);

#endif // sighandler_h

