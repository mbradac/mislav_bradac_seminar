CC=g++
SRCDIR=src
BINDIR=bin
TESTDIR=test
SPIKEDIR=spike
RM=rm
INC=include
CFLAGS=-Wall -std=c++11 -march=native -I $(INC)

sequence: $(SRCDIR)/sequence.cc $(INC)/sequence.h
	$(CC) $(CFLAGS) $(SRCDIR)/sequence.cc -c -o $(BINDIR)/sequence.o

sequence_test: $(TESTDIR)/sequence_test.cc sequence
	$(CC) $(CFLAGS) $(TESTDIR)/sequence_test.cc $(BINDIR)/sequence.o -o $(BINDIR)/sequence_test

check_avx2: $(SPIKEDIR)/check_avx2.cc
	$(CC) $(CFLAGS) $(SPIKEDIR)/check_avx2.cc -o $(BINDIR)/check_avx2

test: sequence_test
	$(BINDIR)/sequence_test

clean:
	$(RM) -rf BINDIR
