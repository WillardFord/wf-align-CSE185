# Benchmarking wf-align against BWA-MEM

For benchmarking I used SARS COV2 reference and fastq files. The files can be found at the following links. I used the ERR10000004 version of the fastq file.

```
fasta:
https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2

fastq:
https://www.ebi.ac.uk/ena/browser/view/PRJEB37886
```

## Running wf-align

To run wf-align using the SARS COV2 data I used the following command:

``
wf-align sequence-4.fasta ERR10000004.fastq -o output.sam -m metric.txt
``

Then by examining the metric file, we can see that (on my macbook air) the alignment took 4.65 seconds and correctly mapped ~18.9% of reads. This low percentage can be attributed to only aligning exact reads.

## Running BWA-MEM

To run BWA-MEM using the same data we have to run two comands listed below. The first generates the neccessary auxiliary structures and the second aligns reads to the reference genome.

```
bwa index sequence-4.fasta
bwa mem sequence-4.fasta  ERR10000004.fastq
```

We can see from the standard output that both runs (on my macbook air) took about 1.964 seconds. Using the following command we can also determine that BWA-MEM only failed to align 3 total reads resulting in 5361/5364 being mapped.

```
cat output_bwa.sam | cut -f4 | grep ^0 | wc -l
```

## Results

BWA-MEM outperforms wf-align in both run time and percentage of reads mapped.

wf-align takes approximately 2.4 times longer than BWA-MEM on my machine, we would expect this to get even worse as we use larger datasets. 

wf-align is significantly beat out in percentage of reads mapped. This can be explained by wf-align only mapping exact matches, whereas BWA-MEM allows close mis-reads, which is relatively common in this dataset.

## Memory usage

Both wf-align and BWA-MEM scale linearly with memory usage. More in-depth analysis still needs to be completed.
