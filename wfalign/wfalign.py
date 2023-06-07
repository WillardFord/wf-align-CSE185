#!/usr/bin/env python
"""
Command-line script to perform alignment of fastq reads to a reference genome

Similar to bwa alignment
"""

import argparse
from . import utils as utils
import os
import sys
import time

def main():
    parser = argparse.ArgumentParser(
        prog="",
        description=""
    )

    # Inputs
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

    # Generate Metics
    total_ref_bps = 0
    time_building = 0
    start_fasta = time.time()

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
                # Generate Auxilury Data Structures
                # ref_data[i] = [Ref_Name, Suffix Array, BurrowsWheelerTransform, \
                #                   LasttoFirst, CountVector]
                build_start = time.time()

                ref_data.append(utils.GENERATE_AUXS(reference, rname))

                build_end = time.time()
                time_building += build_end - build_start
                total_ref_bps += len(reference)
            reference = ""
            rname = line.split()[0][1:] # don't grab > char
        else:
            reference += line.replace("\n","") # Last line might not have \n
    
    # Because we only add aux data structures when we reach the next chromosome,
    # we manually add the last one here which has no subsequent chromosome.
    build_start = time.time()

    ref_data.append(utils.GENERATE_AUXS(reference, rname))

    build_end = time.time()
    time_building += build_end - build_start
    total_ref_bps += len(reference)

    ref_fasta.close()

    # Iterate Benchmarking Results
    end_fasta = time.time()
    time_in_search = 0

    num_seqs = 0
    num_uniquely_aligned = 0
    num_unaligned = 0
    total_seq_bps = 0

    # Perform alignment

    qname = ""  # Line 1 in fasta format
    seq = ""    # Line 2 in fasta format
    quals = ""  # Line 4 in fasta format
    
    line_num = -1
    for line in fastq:
        line_num += 1
        if line_num % 4 == 0:
            qname = line.split()[0][1:]  # Don't grab > character
        elif line_num % 4 == 1:
            seq = line[:-1]     # Don't grab \n character
        elif line_num % 4 == 2:
            continue
        else:
            # can't slice because last line might not have \n character
            quals = line.replace("\n","")   
            m = len(seq)

            
            # Now use all previous info to align
            find_start = time.time()
            loc:int = 0
            for i in range(len(ref_data)):
                # ref_data[i] = [Ref_Name, Suffix Array, BurrowsWheelerTransform, 
                #                   LasttoFirst, CountVector]
                loc, multiple_matches = utils.FIND( 
                    SA=ref_data[i][1], BWT=ref_data[i][2], \
                    LTF=ref_data[i][3], C=ref_data[i][4], \
                    PATTERN=seq, M=m, \
                ) + 1 # sam files are 1 indexed
                if loc > 0:
                    outf.write(utils.GET_ALIGNMENT(
                        QNAME=qname, TEMPLATE=seq, QUAL=quals, POS=loc, \
                        RNAME= ref_data[i][0],
                    ))
                    num_uniquely_aligned += 1 if not multiple_matches else 0
                    break
            find_end = time.time()

            # If no exact match was found, no statement can be made about location
            if loc == 0:
                outf.write(utils.GET_ALIGNMENT(
                        QNAME=qname, TEMPLATE=seq, QUAL=quals, POS=loc, \
                        RNAME= "*",
                ))
                num_unaligned += 1

            # Iterate Benchmarking Results
            
            time_in_search += find_end - find_start 
            total_seq_bps += m
            num_seqs += 1

    outf.close()

    end_align = time.time()


    # Calculate Benchmarking Results

    fasta_time = end_fasta - start_fasta
    aligning = end_align - end_fasta
    '''
    fasta_time      : Time spent loading data from fasta files and building structures
    aligning        : Time spent loading data from fastq files and aligning

    # These two help differentiate IO time from algorithm time
    time_in_search  : Time spent literally searching for matches
    time_building   : Time spent literally building auxilury data structures

    total_ref_bps   : Total length of reference genome in base pairs
    total_seq_bps   : Total length of reads in base pairs

    num_seqs        : Total number of reads
    num_unaligned   : Total number of reads that did not align to any reference chromosome
    num_uniquely_aligned    : Total number of reads that were uniquely aligned
    # This number is slightly inlfated because once we find matches on a single chromosome
    # we will not search the remaining chromosomes. This helps with speed but increases
    # the number of reads uniquely aligned.
    
    '''

    metrics = utils.GET_METRICS(TOT_RDS_LEN=total_seq_bps, TOT_REF_LEN=total_ref_bps, \
                NUM_RDS=num_seqs, NUM_NA=num_unaligned, NUM_UA=num_uniquely_aligned, \
                A_TIME=aligning, AF_TIME=time_in_search, REF_TIME=fasta_time, \
                REFB_TIME=time_building,
                )

    # Output Benchmarking Results

    cur_time = time.strftime("%d:%H:%M")
    outm = open(f"{cur_time}_wf_align_metrics.txt", "w")
    outm.write(metrics)
    outm.close()

if __name__ == '__main__':
    main()
