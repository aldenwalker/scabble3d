CCC=g++
CC=gcc
CFLAGS= -g -Wall
IFLAGS=-I/sw/include
LDFLAGS=-L/sw/lib -lglpk -lgmp -lpthread
all: scabble3d

scabble3d.o: scabble3d.c *.h
	$(CC) $(CFLAGS) $(IFLAGS) -c scabble3d.c `pkg-config --cflags --libs gtk+-2.0`
	
word.o: word.c word.h
	$(CC) $(CFLAGS) $(IFLAGS) -c word.c 
	
matrix.o: matrix.c matrix.h
	$(CC) $(CFLAGS) $(IFLAGS) -c matrix.c 

scl_backend.o: scl_backend.c *.h
	$(CC) $(CFLAGS) $(IFLAGS) -c scl_backend.c `pkg-config --cflags --libs gtk+-2.0`

triangle_and_vertex.o: triangle_and_vertex.c triangle_and_vertex.h
	$(CC) $(CFLAGS) $(IFLAGS) -c triangle_and_vertex.c


.PHONY: exlp-package
exlp-package: 
	cd exlp-package; make

scabble3d: exlp-package scabble3d.o matrix.o word.o scl_backend.o triangle_and_vertex.o
	$(CCC) $(CFLAGS) $(IFLAGS) -o scabble3d scabble3d.o matrix.o word.o scl_backend.o triangle_and_vertex.o exlp-package/*.o `pkg-config --cflags --libs gtk+-2.0` $(LDFLAGS)



clean: 
	rm scabble3d
	rm *.o
	cd exlp-package; rm *.o
