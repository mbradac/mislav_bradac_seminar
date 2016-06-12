CC=g++
SRCDIR=src
BINDIR=bin
TESTDIR=test
SPIKEDIR=spike
SIMDEXMPSDIR=simd_examples
RM=rm
INC=include
CFLAGS=-Wall -std=c++11 -march=native -I $(INC) -O2

$(BINDIR):
	mkdir -p $(BINDIR)

sequence: $(SRCDIR)/sequence.cc $(INC)/sequence.h | $(BINDIR)
	$(CC) $(CFLAGS) $(SRCDIR)/sequence.cc -c -o $(BINDIR)/sequence.o

sequence_test: $(TESTDIR)/sequence_test.cc sequence | $(BINDIR)
	$(CC) $(CFLAGS) $(TESTDIR)/sequence_test.cc $(BINDIR)/sequence.o -o $(BINDIR)/sequence_test

matrix_test: $(TESTDIR)/matrix_test.cc sequence | $(BINDIR)
	$(CC) $(CFLAGS) $(TESTDIR)/matrix_test.cc $(BINDIR)/sequence.o -o $(BINDIR)/matrix_test

matrix_test2: $(TESTDIR)/matrix_test2.cc sequence | $(BINDIR)
	$(CC) $(CFLAGS) $(TESTDIR)/matrix_test2.cc $(BINDIR)/sequence.o -o $(BINDIR)/matrix_test2

matrix_test3: $(TESTDIR)/matrix_test3.cc sequence | $(BINDIR)
	$(CC) $(CFLAGS) $(TESTDIR)/matrix_test3.cc $(BINDIR)/sequence.o -o $(BINDIR)/matrix_test3

check_avx2: $(SPIKEDIR)/check_avx2.cc | $(BINDIR)
	$(CC) $(CFLAGS) $(SPIKEDIR)/check_avx2.cc -o $(BINDIR)/check_avx2

fasta_statistics: $(SPIKEDIR)/fasta_statistics.cc sequence | $(BINDIR)
	$(CC) $(CFLAGS) $(SPIKEDIR)/fasta_statistics.cc $(BINDIR)/sequence.o -o $(BINDIR)/fasta_statistics

string_compare: $(SIMDEXMPSDIR)/string_compare.cc | $(BINDIR)
	$(CC) $(CFLAGS) $(SIMDEXMPSDIR)/string_compare.cc -o $(BINDIR)/string_compare

vector_matrix_4x4: $(SIMDEXMPSDIR)/vector_matrix_4x4.cc | $(BINDIR)
	$(CC) $(CFLAGS) $(SIMDEXMPSDIR)/vector_matrix_4x4.cc -o $(BINDIR)/vector_matrix_4x4

smith_waterman_no_simd: $(SPIKEDIR)/smith_waterman_no_simd.cc $(INC)/sequence.h | $(BINDIR)
	$(CC) $(CFLAGS) $(SPIKEDIR)/smith_waterman_no_simd.cc -c -o $(BINDIR)/smith_waterman_no_simd.o

smith_waterman_no_simd_test: $(TESTDIR)/smith_waterman_no_simd_test.cc smith_waterman_no_simd sequence | $(BINDIR)
	$(CC) $(CFLAGS) $(TESTDIR)/smith_waterman_no_simd_test.cc $(BINDIR)/smith_waterman_no_simd.o $(BINDIR)/sequence.o -o $(BINDIR)/smith_waterman_no_simd_test

smith_waterman: $(SRCDIR)/smith_waterman.cc $(INC)/sequence.h $(INC)/specialized_simd.h | $(BINDIR)
	$(CC) $(CFLAGS) $(SRCDIR)/smith_waterman.cc -c -o $(BINDIR)/smith_waterman.o

search: $(SRCDIR)/search.S | $(BINDIR)
	$(CC) $(SRCDIR)/search.S -c -o $(BINDIR)/search.o

smith_waterman_test: $(TESTDIR)/smith_waterman_test.cc smith_waterman sequence search | $(BINDIR)
	$(CC) $(CFLAGS) $(TESTDIR)/smith_waterman_test.cc $(BINDIR)/smith_waterman.o $(BINDIR)/search.o $(BINDIR)/sequence.o -o $(BINDIR)/smith_waterman_test

test: sequence_test matrix_test matrix_test2 matrix_test3 smith_waterman_no_simd_test smith_waterman_test
	$(BINDIR)/sequence_test
	$(BINDIR)/matrix_test
	$(BINDIR)/matrix_test2
	$(BINDIR)/matrix_test3
	$(BINDIR)/smith_waterman_no_simd_test
	$(BINDIR)/smith_waterman_test

clean:
	$(RM) -rf BINDIR
