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
        Reference Sequence Name. Not usable without LN.
    LN : str
        Reference Sequence Length. Not usable without SN.
    """
    # fasta version of this output
    version = "1.6"

    # Sorting order of the reads in the fasta output
    #   TODO update this when I have a more concrete idea of my output.
    SO = "unknown"
    header = "@HD VN:{ver} SO:{so}\n".format(ver=version, so=SO)

    if SN != None and LN != None:
        header += "@SQ SN:{sn} LN:{ln}\n".format(sn=SN, ln=LN)

    return header

def RADIX_SORT_sufs(sufs:list, Q:int, N:int) -> list:
    """
    Return sorted list of suffix objects given list of suffix objects.
    Sorts using radix sort on suffix.rank[0] and suffix.rank[1].

    Adapted from the following link: https://www.geeksforgeeks.org/radix-sort/

    Parameters
    ----------
    sufs : list
        Set of suffix objects to be sorted
    Q : int
        Which entry should we be sorting by in the suffix object
    N : int
        Length of sufs. len(sufs)
    """
    output = [0] * N

    # Find norm values s.t. each value >= 1
    # Find new alphabet length
    min_value = 10
    max_value = -10
    for i in range(N):
        if sufs[i].rank[Q] > max_value:
            max_value = sufs[i].rank[Q]
        if sufs[i].rank[Q] < min_value:
            min_value = sufs[i].rank[Q]
    alphabet_len = max_value - min_value + 1
    norm = -min_value

    # Sort
    count = [0] * alphabet_len

    for i in range(N):
        count[sufs[i].rank[Q]+norm] += 1

    for i in range(1, alphabet_len):
        count[i] += count[i - 1]

    # Build the output array
    i = N - 1
    while i >= 0:
        index = sufs[i].rank[Q]+norm
        output[count[index] - 1] = sufs[i]
        count[index] -= 1
        i -= 1
    return output

def SUFFIX_ARRAY(REF:str, N:int) -> list:
    """
    Builds a suffix array from a given string and length of string.

    Adapted from the following link: 
        https://www.geeksforgeeks.org/suffix-array-set-2-a-nlognlogn-algorithm

    Parameters
    ----------
    REF : str
        Input string to build suffix array from
    N : str
        Length of REF. len(REF)
    """
    dictionary = {"A":"B", "C":"C", "G":"D","T":"E", "$":"$"}
    suffixes = [suffix() for _ in range(N)]

    for i in range(N):
        suffixes[i].index = i
        suffixes[i].rank[0] = (ord(dictionary[REF[i]]) - ord("A"))
        suffixes[i].rank[1] = \
            (ord(dictionary[REF[i+1]]) - ord("A")) if ((i + 1) < N) else 0
 
    # sort suffixes
    suffixes = RADIX_SORT_sufs(suffixes, 1, N)
    suffixes = RADIX_SORT_sufs(suffixes, 0, N)
    ind = [0] * N 
    k = 4
    while (k < 2 * N):
        rank = 0
        prev_rank = suffixes[0].rank[0]
        suffixes[0].rank[0] = rank
        ind[suffixes[0].index] = 0

        for i in range(1, N):

            if (suffixes[i].rank[0] == prev_rank and
                suffixes[i].rank[1] == suffixes[i - 1].rank[1]):
                prev_rank = suffixes[i].rank[0]
                suffixes[i].rank[0] = rank

            # Otherwise increment rank
            else: 
                prev_rank = suffixes[i].rank[0]
                rank += 1
                suffixes[i].rank[0] = rank
            ind[suffixes[i].index] = i

        for i in range(N):
            nextindex = suffixes[i].index + k // 2
            suffixes[i].rank[1] = suffixes[ind[nextindex]].rank[0] if (nextindex < N) else -1

        suffixes = RADIX_SORT_sufs(suffixes, 1, N)
        suffixes = RADIX_SORT_sufs(suffixes, 0, N) 
        k *= 2

    suffixArr = [0] * N
    for i in range(N):
        suffixArr[i] = suffixes[i].index

    return suffixArr

"""
Class to store suffixes and their indices 
"""
class suffix:
    def __init__(self):
        self.index = 0
        self.rank = [0, 0]

def BURROWS_WHEELER(REF:str, SA:list,  N:int) -> list:
    """
    Builds a Burrows Wheeler Transformed version of the string from a Suffix
    Array version of the string.

    Parameters
    ----------
    REF : str
        Input string to build bwt from
    SA : list
        Suffix Array of String
    N : int
        Length of REF. len(REF)
    """
    bwt = [0] * N
    for i in range(N):
        bwt[i] = REF[SA[i]-1]
    return bwt

def RADIX_SORT(BWT:list, N:int) -> list:
    """
    Return sorted list of suffix objects given list of suffix objects.
    Sorts using radix sort on suffix.rank[0] and suffix.rank[1].

    Adapted from the following link: https://www.geeksforgeeks.org/radix-sort/

    Parameters
    ----------
    BWT : list
        Set of characters to be sorted
    N : int
        Length of BWT. len(BWT)
    """
    output = [0] * N

    # Find norm values s.t. each value >= 1
    # Find new alphabet length
    min_value = "!"
    max_value = "Z"
    for i in range(N):
        if ord(BWT[i]) - ord("A") > ord(max_value) - ord("A"):
            max_value = ord(BWT[i]) - ord("A")
        if ord(BWT[i]) - ord("A") < ord(min_value) - ord("A"):
            min_value = ord(BWT[i]) - ord("A")
    alphabet_len = ord(max_value) - ord(min_value) + 1
    norm = -ord(min_value)

    # Sort
    count = [0] * alphabet_len

    for i in range(N):
        count[ord(BWT[i]) - ord("A") +norm] += 1

    for i in range(1, alphabet_len):
        count[i] += count[i - 1]

    # Build the output array
    i = N - 1
    while i >= 0:
        index = ord(BWT[i]) - ord("A") +norm
        output[count[index] - 1] = BWT[i]
        count[index] -= 1
        i -= 1
    return output

def FIND() -> int: #TODO
    """
    Returns int representing first index of first character of
    instance of sequence in reference

    Parameters
    ----------
    BWT : list
        Set of characters to be sorted
    N : int
        Length of BWT. len(BWT)
    """
    # TODO
    pass

def GET_ALIGNMENT(QNAME:str) -> str: #TODO
    """
    Returns tabular data of string alignment information according to 
    section 1.4 of the SAM file format specs:
    https://samtools.github.io/hts-specs/SAMv1.pdf

    Parameters
    ---------- TODO
    QNAME : str
        Query Name that corresponds to the read id in the fastq file input
    N : int
        Length of BWT. len(BWT)
    """
    output = ""
    # 1 QNAME String 
    output += QNAME + "\t"
    # 2 FLAG Int

    # 3 RNAME String
    
    # 4 POS Int
    
    # 5 MAPQ Int
    
    # 6 CIGAR String
    
    # 7 RNEXT String
    
    # 8 PNEXT Int
    
    # 9 TLEN Int
    
    # 10 SEQ String
    
    # 11 QUAL


    return output + "\n"

if __name__ == "__main__":
    string = "TACAGTATCGA"
    
    print(SUFFIX_ARRAY(string, len(string)))