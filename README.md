# wf-align (CSE185 Project Demo)

(Work in progress!)

This is a demonstration project for CSE185. It implements a smaller, simpler version of bwa-backtrack. See the [BWA](https://bio-bwa.sourceforge.net) page for more details.

# Install instructions

Installation doesn't require any additional libraries.

Once required libraries are installed, you can install `wf-align` with the following command:

```
python setup.py install
```

Note: if you do not have root access, you can run the commands above with additional options to install locally:
```
python setup.py install --user
```

If the install was successful, typing `wf-align --help` should show a useful message.

# Basic usage

The basic usage of `wf-align` is:

```
wf-align reference.fa reads.fq [-o output.sam] [other options]
```

To run `wf-align` on a small test example (using files in this repo):
```
wf-align example-files/test_reference.fa example-files/test_reads.fastq 
```

This should produce the output below:
```
hello
```

To compare to output of `bwa mem`, run:
```
bwa mem example-files/test_reference.fa example-files/test_reads.fastq
```

# wf-align options

There are 2 required inputs to `wf-align`, a reference fasta file and a fastq file full of reads. Users may additionally specify the options below:

* `-o FILE`, `--output FILE`: Write output to file. By default, output is written to stdout.

# File format

The output file format is the same as the bwa mem method, a sam file. See: https://samtools.github.io/hts-specs/SAMv1.pdf

# Testing

To run tests:
```
# Run command line tests
sh tests/cmdline_tests.sh

# Run unit tests
python -m pytest --cov=.
```

# Methodology

I used several different resources to come up with my own hand-stitched methodology:

# Sources

https://www.youtube.com/watch?v=byNR4CbYiPQ&list=PLQ-85lQlPqFP-q0Ig_GdjWohbC9C3f8tz&index=8

https://www.cs.cmu.edu/~15451-f18/lectures/lec26-suffarray.pdf

https://www.cs.cmu.edu/~15451-f18/lectures/lec25-bwt.pdf

# Contributors

This repository was generated by Willard Ford, with (much) inspiration from the [CSE 185 Example Repository](https://github.com/gymreklab/cse185-demo-project#readme) and the work of my fellow students.

Please submit a pull request with any corrections or suggestions.
