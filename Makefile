CFLAGS=-march=znver3 -Ofast -flto -std=gnu99 -W -pedantic -Wcast-qual \
	-Wpointer-arith -Wcast-align -fno-schedule-insns -fschedule-insns2 \
	-fstrict-aliasing -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer \
	-I /home/hse/ivointsev/laser/nodes/AMD74F3/fftw/include

LDFLAGS=-march=znver3 -Ofast -flto -lm -lpthread \
	-lgsl -lgslcblas -L /home/hse/ivointsev/laser/nodes/AMD74F3/fftw/include -lfftw3 -lfftw3_threads

CC=gcc

all: atmosphere initial analyzer stretch compare

%.o: %.c 
	$(CC) -c $(CFLAGS) -o $@ $<

atmosphere: atmosphere.o dataio.o my_thr_lib.o histogram.o filter.o correlator.o
	$(CC) -Wall -o $@ atmosphere.o my_thr_lib.o dataio.o histogram.o filter.o correlator.o $(LDFLAGS)

analyzer: analyzer.o dataio.o my_thr_lib.o 
	$(CC) -Wall -o $@ analyzer.o my_thr_lib.o dataio.o $(LDFLAGS)

initial: initial.o dataio.o my_thr_lib.o histogram.o filter.o correlator.o
	$(CC) -Wall -o $@ initial.o my_thr_lib.o dataio.o $(LDFLAGS)

stretch: stretch.o dataio.o 
	$(CC) -Wall -o $@ stretch.o dataio.o $(LDFLAGS)

compare: compare.o dataio.o 
	$(CC) -Wall -o $@ compare.o dataio.o $(LDFLAGS)

clean: 
	rm -f *.o 
distclean: clean 
	rm atmosphere analyzer initial stretch compare; 
	./rm_data.sh; rm fftw.wisdom; rm *.cdata; rm LOG

