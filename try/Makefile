CC = cc
FC = f77
CFLAGS = -O2
FFLAGS = -O2

SRC = src
LIB = ../lib/build/libfilter.a

all: test

test.o: $(SRC)/test.c
	$(CC) $(CFLAGS) -c $(src)/test.c

test: test.o $(LIB)
	$(FC) $(FFLAGS) -o test test.o $(LIB)

playmatrix2: playmatrix2.o $(LIB)
	$(FC) $(FFLAGS) -o playmatrix2 playmatrix2.o $(LIB)

clean:
	rm -f test.o test

