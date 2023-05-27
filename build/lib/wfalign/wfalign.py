#!/usr/bin/env python
"""
Command-line script to perform alignment of fastq reads to a reference genome

Similar to bwa-mem
"""


import argparse
from . import utils as utils
# from mypileup import __version__ How does this work?
import os
import sys

def main():
    parser = argparse.ArgumentParser(
        prog="",
        description=""
    )

    # Inputs

    # TODO What happens if I call this method without providing sample methods.
    parser.add_argument("fasta-ref",\
                        help="This is a fasta file representing the reference genome to "\
                        "align your reads back to. For more "\
                        "information about the fasta file format: "\
                        "https://zhanggroup.org/FASTA/", \
                        metavar="FA FILE", type = str)
    
    parser.add_argument("fastq", help="Fastq file containing reads to be aligned to the "\
                        "reference genome. For more information "\
                        "about the fastq file format: "\
                        "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217/", \
                        metavar="FASTQ FILE", type = str)
    
    # Output
    parser.add_argument("-o", "--output", help="Output file path in the sam format. "\
                        "This file will be overwritten if it already exists. "\
                        "For more information about the SAM file format: "\
                        "https://samtools.github.io/hts-specs/SAMv1.pdf\n"
                        "Default: stdout", metavar="FILE", type = str, required=False)

    # Other options

    # Ideas:
    #   bam output

    # Parse args
    args = parser.parse_args()


    # Set up output file
    if args.out == None:
        outf = sys.stdout
    else: outf = open(args.out, "w")

    # Load Fastq Reads
    if not os.path.exists(args.fastq):
        utils.ERROR("{fastq} does not exist".format(fastq=args.fastq))
    if args.fastq[-6:] != ".fastq" and args.fastq[-3:] != ".fa":
        utils.ERROR("{fastq} has wrong file ending. Is it in the correct format? Refer to: " \
                    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217"\
                    .format(fastq=args.fastq))
    rds_fastq = open(args.fastq, "r")

    # Load FASTA Reference
    if not os.path.exists(args.fasta_ref):
        utils.ERROR("{fasta} does not exist".format(fasta=args.fasta_ref))
    #TODO edit:
    if args.fastq[-6:] != ".fastq" and args.fastq[-3:] != ".fa":
        utils.ERROR("{fasta} has wrong file ending. Is it in the correct format? Refer to: " \
                    "https://zhanggroup.org/FASTA/"\
                    .format(fasta=args.fasta_ref))
    ref_fasta = open(args.fastq, "r")

    # Perform alignment

     # Write header lines to output
    #   TODO this is where I can specify lots of additional options.
    #   Simply grab codes from https://samtools.github.io/hts-specs/SAMv1.pdf
    header = utils.GET_FASTA_HEADER()
    outf.write(header)

    # Create BWT array of ref # TODO
    fasta_header = ref_fasta.readline().strip()
    reference = ref_fasta.read().strip()
    print(reference)
    # Create FM Indeex

    print("hello")

   
   

    # Align genomes 

    # 


if __name__ == '__main__':
    main()
