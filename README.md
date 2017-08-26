# qgrs-cpp
C++ implementation of QGRS mapping algorithm - described [here](http://bioinformatics.ramapo.edu/QGRS/credits.php).  This program differs from the actual algorithm used by QGRS Mapper recarding overlapping motifs and the maximum length of the motif.  In particular, a maximum of 45nt is applied for all motifs with 3+ tetrads, while a maximum of 30nt is applied for motifs with 2 tetrads.  The details of this, and dealing with overlapping motifs is described in detail in the following papers:

- Frees, S., Crum, M., Menendez, C., Bagga, P.S. (2014) ["QGRS-Conserve: A Computational Method for Discovering Evolutionarily Conserved G-quadruplex Motifs"](https://humgenomics.biomedcentral.com/articles/10.1186/1479-7364-8-8). Human Genetics, BioMed Central, Vol 8 (8)
- Menendez, C., Frees, S., and Bagga, P. (2012) [QGRS-H Predictor: A Web Server for Predicting Homologous Quadruplex forming G-Rich Sequence Motifs in Nucleotide Sequences](https://academic.oup.com/nar/article/40/W1/W96/1074452/QGRS-H-Predictor-a-web-server-for-predicting). Nucleic Acids Res. doi: 10.1093/nar/gks42240: W96-W103.

## Getting the code
You can clone this repository with the following command (you need git installed on your machine!).
```
git clone https://github.com/freezer333/qgrs-cpp
```
Alternatively, you can download a zip file [here](https://github.com/freezer333/qgrs-cpp/archive/master.zip)

## Building the code
Source code is found in the `/src` directory.  You can build the program with any modern (C++11) compiler, on any platform/OS - there are no external dependencies.

If you are using maxOS or any UNIX/Linux-type system, you can build the program using `make`.  If you have not configured your system to with `make`, you can compile the source code like so:

```
g++ -o qgrs src/default.cpp src/qgrs.cpp --std=c++11
```

The code has not been tested on Windows.

## Using the tool
Once you've build the executable, you can access the help segment of the program using the `-h` flag.

```
./qgrs -h
```

The following command line options are supported:

```
-h                    help (this)
-csv                  output in csv mode
-json                 output in JSON mode
-i [input filename]   read input from a file as specified (defaults to stdin)
-o [output filename]  write output to file as specified (defaults to stdout)
-t [n]                filter output to QGRS with number tetrads equal to or greater than n (defaults to 2)
-s [n]                filter output to QGRS with g-score equal to or greater than n (defaults to 17)
-g [c]                replace all G characters within tetrads with given character.
-v                    include overlapping QGRS in output (very large outputs may be generated)
-notitle              omit column headings in output (does not apply for JSON)
```

### Manually entering sequence data
The simplest way to run it from the command line is using all the defaults.  To do this, simply execute `bin/qgrs`.  

You'll need to enter the sequence data - you can enter multiple lines of characters.  To allow the program to start (using the sequence), you must signal EOF, by typing Ctrl+D at the command line.  You will see output generated.

```
bin/qgrs
GGAGGAGGAGG
^D
1           0  3  6  9  2  21  GGAGGAGGAGG
```

### Understanding the output
By default, output is in columns - with the following values:

```
Column 1 - ID (order found in sequence).  
Column 2 - Position of the start of the first tetrad (relative to beginning of input sequence)
Column 3 - Position of the start of the second tetrad (relative to beginning of input sequence)
Column 4 - Position of the start of the third tetrad (relative to beginning of input sequence)
Column 5 - Position of the start of the fourth tetrad (relative to beginning of input sequence)
Column 6 - Number of tetrads
Column 7 - G-Score
Column 8 - Sequence
```

### Using sequence data from a file
On a system that supports input redirection, you can use sequence data from a file as such:

```
bin/qgrs <samples/input.txt
```

The program will map QGRS from the input.txt file provided with the source code distribution.

Alternatively, you can use the `-i` flag and specify the input name - which will give you the same results.

```
bin/qgrs -i samples/input.txt
```

### Saving output to files
You can redirect stdout when using the tool, or you can use the -o flag to have the mapping results stored in a resulting file.

```
bin/qgrs -i samples/input.txt >samples/output.txt
```

```
bin/qgrs -i samples/input.txt -o samples/output.txt
```

### Saving output as a CSV /JSON files
The tool supports output in CSV format instead of the column-based output used as default.

```
bin/qgrs -csv -i samples/input.txt > samples/output.csv
```

Alternative, you can output as a JSON document

```
bin/qgrs -json -i samples/input.txt > samples/output.json
```

### Filtering mapping results
QGRS can be filtered by minimum tetrad and minimum gscore using the `-t` and `-s` flags respectively.

The following will map QGRS from `samples/nras.txt`, filtering such that only QGRS with 3 tetrads and scores above 50 are presented in the output.

```
bin/qgrs -i samples/nras.txt -t 3 -s 50
1             14    18    23    29  3  62  GGGAGGGGCGGGTCTGGG
2             81    97   119   123  3  52  GGGACTCAGGCGCCTGGGGCGCCGACTGATTACGTAGCGGGCGGG
3           6690  6695  6700  6704  3  63  GGGGGGGGTTGGGGGGG
```

These flags can be specified in any order.

### Overlapping QGRS
Most QGRS found have alternative forms found within them - or overlapping them.  More information about this concept can be found [here](http://bioinformatics.ramapo.edu/QGRS/help_overlaps.php).

By default, the tools does not output overlapping sequences.  You may use the `-v` option to enable their output.

```
bin/qgrs -i samples/nras.txt -t 3 -s 50 -v
1             14    18    23    29  3  62  GGGAGGGGCGGGTCTGGG
1.1           14    19    23    29  3  62  GGGAGGGGCGGGTCTGGG
2             81    97   119   123  3  52  GGGACTCAGGCGCCTGGGGCGCCGACTGATTACGTAGCGGGCGGG
2.1           81    96   119   123  3  51  GGGACTCAGGCGCCTGGGGCGCCGACTGATTACGTAGCGGGCGGG
3           6690  6695  6700  6704  3  63  GGGGGGGGTTGGGGGGG
3.1         6690  6694  6700  6703  3  62  GGGGGGGGTTGGGGGG
3.2         6690  6694  6700  6704  3  62  GGGGGGGGTTGGGGGGG
3.3         6690  6694  6701  6704  3  61  GGGGGGGGTTGGGGGGG
3.4         6690  6695  6700  6703  3  62  GGGGGGGGTTGGGGGG
3.5         6690  6695  6701  6704  3  62  GGGGGGGGTTGGGGGGG
3.6         6691  6695  6700  6703  3  62  GGGGGGGTTGGGGGG
3.7         6691  6695  6700  6704  3  63  GGGGGGGTTGGGGGGG
3.8         6691  6695  6701  6704  3  62  GGGGGGGTTGGGGGGG
```

Note the ID numbers - for example, QGRS #3 has 8 overlapping sequence variants - 3.1-3.8, each of which have lower g-scores.  Filters (tetrads, g-scores) are applied to overlaps as well.

### Substituting tetrad characters
To make it easier to read output, it's often nice to see which G characters are actually involved in the QGRS.  You can use the `-s` flag to substitute tetrad G's with another character.

```
bin/qgrs -i samples/nras.txt -t 3 -s 50 -v -g @
1             14    18    23    29  3  62  @@@A@@@GC@@@TCT@@@
1.1           14    19    23    29  3  62  @@@AG@@@C@@@TCT@@@
2             81    97   119   123  3  52  @@@ACTCAGGCGCCTG@@@CGCCGACTGATTACGTAGC@@@C@@@
2.1           81    96   119   123  3  51  @@@ACTCAGGCGCCT@@@GCGCCGACTGATTACGTAGC@@@C@@@
3           6690  6695  6700  6704  3  63  @@@GG@@@TT@@@G@@@
3.1         6690  6694  6700  6703  3  62  @@@G@@@GTT@@@@@@
3.2         6690  6694  6700  6704  3  62  @@@G@@@GTT@@@G@@@
3.3         6690  6694  6701  6704  3  61  @@@G@@@GTTG@@@@@@
3.4         6690  6695  6700  6703  3  62  @@@GG@@@TT@@@@@@
3.5         6690  6695  6701  6704  3  62  @@@GG@@@TTG@@@@@@
3.6         6691  6695  6700  6703  3  62  @@@G@@@TT@@@@@@
3.7         6691  6695  6700  6704  3  63  @@@G@@@TT@@@G@@@
3.8         6691  6695  6701  6704  3  62  @@@G@@@TTG@@@@@@
```

Only a single character should be used as the substitute.
