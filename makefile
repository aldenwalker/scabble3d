CCC=g++
CC=gcc
CFLAGS= -g -Wall
IFLAGS=-I/sw/include
LDFLAGS=-L/sw/lib -lglpk -lgmp
all: scabble3d

scabc: scabble3d.c matrix.c word.c scl_backend.c triangle_and_vertex.c
	$(CC) $(CFLAGS) $(IFLAGS) -c scabble3d.c matrix.c word.c scl_backend.c triangle_and_vertex.c `pkg-config --cflags --libs gtk+-2.0`

scabble3d: scabc
	cd exlp-package; make
	$(CCC) $(CFLAGS) $(IFLAGS) -o scabble3d scabble3d.o matrix.o word.o scl_backend.o triangle_and_vertex.o exlp-package/*.o `pkg-config --cflags --libs gtk+-2.0` $(LDFLAGS)

clean: 
	rm scabble3d
	rm *.o
	cd exlp-package; rm *.o
