"""
Utilities for align
"""
import sys

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def ERROR(msg:str) -> None:
    """
    Print an error message and die

    Parameters
    ----------
    msg : str
       Error message to print
    """
    sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + "{msg}\n".format(msg=msg) )
    sys.exit(1)

def GET_FASTA_HEADER(SN:str = None, LN:str = None) -> str:
    """
    Returns the header of the output fasta file in version 1.6 format.
    More details: https://samtools.github.io/hts-specs/SAMv1.pdf

    Parameters
    ----------
    SN : str
        Reference Sequence Name. Not usable without VN.
    LN : str
        Reference Sequence Length. Not usable without SN.
    """
    # fasta version of this output
    version = "1.6"

    # Sorting order of the reads in the fasta output
    #   TODO update this when I have a more concrete idea of my output.
    SO = "unknown"
    header = "@HD VN:{ver} SO:{so}\t".format(ver=version, so=SO)

    if SN != None and LN != None:
        header += "@SQ SN:{sn} LN:{ln}\t".format(sn=SN, ln=LN)

    return header