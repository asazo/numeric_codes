CC=g++
CFLAGS=-Wall -g

all: main


p1.o: p1.cpp wavesolver.h
	$(CC) $(CFLAGS) -c p1.cpp

p1: p1.o wavesolver.h
	$(CC) $(CFLAGS) -o p1 p1.o wavesolver.h

p2.o: p2.cpp ellipticsolver.h
	$(CC) $(CFLAGS) -c p2.cpp

p2: p2.o ellipticsolver.h
	$(CC) $(CFLAGS) -o p2 p2.o ellipticsolver.h

p3.o: p3.cpp solver.h
	$(CC) $(CFLAGS) -c p3.cpp

p3: p3.o solver.h
	$(CC) $(CFLAGS) -o p3 p3.o solver.h



clean:
	rm -v *.o p3 p2 p1
