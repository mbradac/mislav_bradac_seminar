# Smith Waterman using AVX2

This is implementation of Smith Waterman algorithm using AVX2 instruction set.

Core loop of algorithm is very similiar to Rognes' core loop in SWIPE. SWIPE uses only SSE registers (128 bits) while this implementation uses AVX registers (256 bits) thus twice as much sequences can be processed in parallel. You can take a look at SWIPE here: <https://github.com/torognes/swipe>.

Some implementation details were inspired by opal, but opal doesn't implement core loop in assembly nor does any loop unrolling thus this implementation is on some examples ~25% faster. You can take a look at opal here: <https://github.com/Martinsos/opal>.

While opal provides maximal 32 bit score range and SWIPE 64 bit score range, this implementation's currently supported maximal score range is 16 bits. (In the near future I will probably add possibility of maximal 32 bit score range.)

### Requirements

Support for AVX2 and g++ compiler that supports c++11.

### Usage

You can use this implementation of algorithm as a library or as a command line tool. To use it as a command line tool you should `make command_line_smith_waterman` and then call `bin/command_line_smith_waterman`. List of tool's arguments:

```
  -q FILE     FILE is query file
  -d FILE     FILE is database file
  -m FILE     FILE is matrix file
  -r N        N is gap extend penalty        [default: 1]
  -e N        N is gap open_extend penalty   [default: 3]
  -x N        N is either 0, 1 or 2          [default: 0]
                   0 - dynamic score range
                   1 - 8bit score range
                   2 - 16bit score range
  -s          silent - don't print results
```

Query file and database file are standard FASTA files, but they must not contain comments starting with semicolon (';') nor empty lines. Comments starting with greater-than ('>') are expected and treated as a name and description of next sequence. To see format of matrix file take a look at any matrix in `test_data/matrix`.
Example of usage:
```
bin/command_line_smith_waterman -q test_data/query/P18080.fasta \
	-d test_data/database/uniprot_sprot12071.fasta \
	-m test_data/matrix/blosum50.mat -r 1 -e 3 -x 0
```
