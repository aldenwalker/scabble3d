CCC=g++
CC=gcc
CFLAGS= -O3#-g -Wall
IFLAGS=-I/sw/include
LDFLAGS=-L/sw/lib -lglpk -lgmpxx -lgmp
all: scabble

scabcc: scabble.cc word.cc rational.cc lp.cc
	$(CCC) $(CFLAGS) $(IFLAGS) -c scabble.cc word.cc rational.cc lp.cc

scabc: matrix.c
	$(CC) $(CFLAGS) $(IFLAGS) -c matrix.c

scabble: scabcc scabc
	cd exlp-package; make
	$(CCC) $(CFLAGS) $(IFLAGS) -o scabble scabble.o word.o rational.o matrix.o lp.o exlp-package/*.o $(LDFLAGS)


clean: 
	rm scabble
	rm *.o
	cd exlp-package; rm *.o
