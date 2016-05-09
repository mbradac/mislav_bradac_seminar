CC=g++
SRCDIR=src
BINDIR=bin
TESTDIR=test
SPIKEDIR=spike
SIMDEXMPSDIR=simd_examples
RM=rm
INC=include
CFLAGS=-Wall -std=c++11 -march=native -I $(INC) -O2

sequence: $(SRCDIR)/sequence.cc $(INC)/sequence.h
	$(CC) $(CFLAGS) $(SRCDIR)/sequence.cc -c -o $(BINDIR)/sequence.o

sequence_test: $(TESTDIR)/sequence_test.cc sequence
	$(CC) $(CFLAGS) $(TESTDIR)/sequence_test.cc $(BINDIR)/sequence.o -o $(BINDIR)/sequence_test

matrix_test: $(TESTDIR)/matrix_test.cc sequence
	$(CC) $(CFLAGS) $(TESTDIR)/matrix_test.cc $(BINDIR)/sequence.o -o $(BINDIR)/matrix_test

matrix_test2: $(TESTDIR)/matrix_test2.cc sequence
	$(CC) $(CFLAGS) $(TESTDIR)/matrix_test2.cc $(BINDIR)/sequence.o -o $(BINDIR)/matrix_test2

matrix_test3: $(TESTDIR)/matrix_test3.cc sequence
	$(CC) $(CFLAGS) $(TESTDIR)/matrix_test3.cc $(BINDIR)/sequence.o -o $(BINDIR)/matrix_test3

check_avx2: $(SPIKEDIR)/check_avx2.cc
	$(CC) $(CFLAGS) $(SPIKEDIR)/check_avx2.cc -o $(BINDIR)/check_avx2

fasta_statistics: $(SPIKEDIR)/fasta_statistics.cc sequence
	$(CC) $(CFLAGS) $(SPIKEDIR)/fasta_statistics.cc $(BINDIR)/sequence.o -o $(BINDIR)/fasta_statistics

string_compare: $(SIMDEXMPSDIR)/string_compare.cc
	$(CC) $(CFLAGS) $(SIMDEXMPSDIR)/string_compare.cc -o $(BINDIR)/string_compare

test: sequence_test matrix_test matrix_test2 matrix_test3
	$(BINDIR)/sequence_test
	$(BINDIR)/matrix_test
	$(BINDIR)/matrix_test2
	$(BINDIR)/matrix_test3

clean:
	$(RM) -rf BINDIR
