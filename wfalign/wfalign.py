#!/usr/bin/env python
"""
Command-line script to perform alignment of fastq reads to a reference genome

Similar to bwa alignment
"""

import argparse
from . import utils as utils
import os
import sys

def main():
    parser = argparse.ArgumentParser(
        prog="",
        description=""
    )

    # Inputs

    # TODO What happens if I call this method without providing sample methods.
    parser.add_argument("fasta_ref",\
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
    if args.output == None:
        outf = sys.stdout
    else: outf = open(args.output, "w")

    # Load FASTA Reference
    if not os.path.exists(args.fasta_ref):
        utils.ERROR("{fasta} does not exist".format(fasta=args.fasta_ref))
    if args.fasta_ref[-6:] != ".fasta" and args.fasta_ref[-3:] != ".fa":
        utils.ERROR("{fasta} has wrong file ending. Is it in the correct format? Refer to: " \
                    "https://zhanggroup.org/FASTA/"\
                    .format(fasta=args.fasta_ref))
    ref_fasta = open(args.fasta_ref, "r")

    # Load Fastq Reads
    if not os.path.exists(args.fastq):
        utils.ERROR("{fastq} does not exist".format(fastq=args.fastq))
    if args.fastq[-6:] != ".fastq" and args.fastq[-3:] != ".fa":
        utils.ERROR("{fastq} has wrong file ending. Is it in the correct format? Refer to: " \
                    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217"\
                    .format(fastq=args.fastq))
    fastq = open(args.fastq, "r")


    # Write header lines to output
    # Simply grab codes from https://samtools.github.io/hts-specs/SAMv1.pdf
    sam_header = utils.GET_FASTA_HEADER()
    outf.write(sam_header)
    del sam_header


    # Read in fasta file.

    rname = "*"
    reference = ""

    # list of lists: [rname, reference, sa, bwt, f_bwt]
    # Reference Name, Reference, Suffix Array, Burrows Wheler Transform, 
    # First Column of BWT Matrix, Last to First
    ref_data = []

    for line in ref_fasta.readlines():
        if line[0] == ">":
            if reference != "":
                N = len(reference)
                # Construct Suffix Array
                sa = utils.SUFFIX_ARRAY(reference, N)

                # Create BWT array of ref
                bwt = utils.BURROWS_WHEELER(reference, sa, N)
                # Reconstruct first column of bwt matrix
                f_bwt = utils.RADIX_SORT(bwt, N)
                ltf = utils.LAST2FIRST(bwt, f_bwt, N)
                ref_data.append([rname, reference, sa, bwt, f_bwt, ltf])
            reference = ""
            rname = line[1:-1] # don't grab \n char
        else:
            reference += line.replace("\n","") # Last line might not have \n
    
    # Because we only add aux data structures when we reach the next chromosome,
    # we manually add the last one here which has no subsequent chromosome.
    N = len(reference)
    # Construct Suffix Array
    sa = utils.SUFFIX_ARRAY(reference, N)
    # Create BWT array of ref
    bwt = utils.BURROWS_WHEELER(reference, sa, N)
    # Reconstruct first column of bwt matrix
    f_bwt = utils.RADIX_SORT(bwt, N)
    ltf = utils.LAST2FIRST(bwt, f_bwt, N)
    ref_data.append([rname, reference, sa, bwt, f_bwt, ltf])

    ref_fasta.close()

    # Perform alignment TODO

    qname = ""  # Line 1
    seq = ""    # Line 2
    quals = ""  # Line 4
    
    line_num = -1
    for line in fastq:
        line_num += 1
        if line_num % 4 == 0:
            qname = line[1:-1]  # Don't grab \n character
        elif line_num % 4 == 1:
            seq = line[:-1]     # Don't grab \n character
        elif line_num % 4 == 2:
            continue
        else:
            # can't slice because last line might not have \n character
            quals = line.replace("\n","")   

            # Now use all previous info to align
            loc:int = 0
            for i in range(len(ref_data)):
                # ref_data = [rname, reference, sa, bwt, f_bwt, ltf]
                loc = utils.FIND( # Only finds exact matches
                    REF=ref_data[i][1], SA=ref_data[i][2], 
                    BWT=ref_data[i][3], F_BWT=ref_data[i][4]
                    LTF=ref_data[i][5]
                )
                if loc > 0:
                    outf.write(utils.GET_ALIGNMENT(
                        QNAME=qname, TEMPLATE=seq, QUAL=quals, POS=loc,
                        RNAME= ref_data[0],
                    ))
                    break
            # If no exact match was found
            if loc == 0:
                outf.write(utils.GET_ALIGNMENT(
                        QNAME=qname, TEMPLATE=seq, QUAL=quals, POS=loc,
                        RNAME= "*",
                    ))

    ref_fasta.close()

if __name__ == '__main__':
    main()
