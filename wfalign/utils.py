"""
Utilities for align
"""
import sys
import time

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

def GENERATE_AUXS(REF:str, RNAME:str) -> list:
    """
    Calls other functions to generate all auxilury data structures needed for pattern matching. 
    Those are as follows:


    Parameters
    ----------
    REF : str
        String to construct auxiliury data structures on.
    RNAME : str
        Name of reference string. Often Chr1 or simply 1
    """

    REF += "$"
    N = len(REF)
    # Construct Suffix Array
    sa = SUFFIX_ARRAY(REF, N)

    # Create BWT array of ref
    bwt = BURROWS_WHEELER(REF, sa, N)
    # Reconstruct first column of bwt matrix
    f_bwt, ltf = RADIX_SORT_L2F(bwt, N)
    # Find first occurrences of each character in f_bwt

    c = [0, 0, 0, 0, 0, N] # A, C, G, N, T + 1 extra to get where T's end
    j = 0
    for i in range(N):
        if j == 0 and f_bwt[i] == "A":
            c[j] = i
            j += 1
        if j == 1 and f_bwt[i] == "C":
            c[j] = i
            j += 1
        if j == 2 and f_bwt[i] == "G":
            c[j] = i
            j += 1
        if j == 3 and f_bwt[i] == "N":
            c[j] = i
            j += 1
        if j == 4 and f_bwt[i] == "T":
            c[j] = i
            j += 1
    return [RNAME, sa, bwt, ltf, c]

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
    dictionary = {"A":"B", "C":"C", "G":"D","T":"F", "$":"$", "N":"E"}
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

"""
Class to store original indices and value. Used for generating last to first vector
"""
class l2f_aux:
    def __init__(self, i:int, val:str):
        self.index = i
        self.rank = val

def RADIX_SORT_L2F(BWT:list, N:int,) -> list:
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
    l2fs = [l2f_aux(i, x) for i, x in enumerate(BWT)]
    output = [0] * N

    # Find norm values s.t. each value >= 1
    # Find new alphabet length
    min_value = "!"
    max_value = "Z"
    for i in range(N):
        # We need to generate last to first array as well.
        if ord(l2fs[i].rank) - ord("A") > ord(max_value) - ord("A"):
            max_value = ord(l2fs[i].rank) - ord("A")
        if ord(l2fs[i].rank) - ord("A") < ord(min_value) - ord("A"):
            min_value = ord(l2fs[i].rank) - ord("A")
    alphabet_len = ord(max_value) - ord(min_value) + 1
    norm = ord(min_value)

    # Sort
    count = [0] * alphabet_len

    for i in range(N):
        count[ord(l2fs[i].rank) - ord("A") + norm] += 1

    for i in range(1, alphabet_len):
        count[i] += count[i - 1]

    # Build the output and l2f array
    l2f = [0]*N
    i = N - 1
    while i >= 0:
        index = ord(l2fs[i].rank) - ord("A") + norm
        output[count[index] - 1] = l2fs[i].rank
        l2f[l2fs[i].index] = count[index] - 1
        count[index] -= 1
        i -= 1
    return output, l2f

def FIND(SA:list, BWT:list, LTF:list, C:list,
         PATTERN:str, M:int,) -> int:
    """
    Returns int representing first index of first character of
    instance of sequence in reference. In the output value here is 1 indexed
    to align with SAM file consistencies Returns -1 if no matches were found.

    Returns boolean. True indicates multiple alignments found. 
    False indicates only one alignment was found

    Parameters
    ----------
    SA : list
        Suffix Array representation of REF
    BWT : list
        Burrows Wheeler transformation representation of REF
    LTF : list
        Last to first vector of Burrows Wheeler transformation
    C : list
        List indicating starting indicies for each character in order
        A, C, G, T
    PATTERN : str
        The pattern to match onto REF
    M : int
        Length of PATTERN, len(PATTERN)
    """

    dic = {"A":0, "C":1, "G":2, "N":3, "T":4}
    min = C[dic[PATTERN[M-1]]]
    max = C[dic[PATTERN[M-1]]+1]

    #multiple_matches = False
    matches = range(min, max)
    for i in range(1, M):
        new_matches = []
        for j in matches:
            char = BWT[j]
            if char == PATTERN[-i-1]:
                new_matches.append(LTF[j])
        # If no matches found return -1
        if len(new_matches) == 0:
            return -1 #, multiple_matches
        matches = new_matches

    # Denote if multiple matches were found
    #if matches > 1:
    #    multiple_matches = True
    # Return only first match for simplicity
    return SA[matches[0]] #, multiple_matches

def GET_ALIGNMENT(QNAME:str, TEMPLATE:str, QUAL:str, POS:int,
                  RNAME:str) -> str: #TODO
    """
    Returns tabular data of string alignment information according to 
    section 1.4 of the SAM file format specs:
    https://samtools.github.io/hts-specs/SAMv1.pdf

    Parameters
    ----------
    QNAME : str
        Query Name that corresponds to the read id in the fastq file input
    TEMPLATE : str
        This string stores the read that we aligned to the reference genome
    QUAL : str
        Quality scores taken directly from the fastq file.
        e.g. if given illumina reads, quality score decoding found from Illumina
    POS : int
        Mapped alignment location in the corresponding RNAME
    RNAME : int
        Reference Name. Often of the format Chr1 or 1.
    """
    output = ""
    # 1 QNAME String 
    output += QNAME + "\t"
    # 2 FLAG Int TODO flags need to be determined
    output += str(0) + "\t"
    # 3 RNAME String
    output += RNAME + "\t"
    # 4 POS Int
    output += str(POS) + "\t"
    # 5 MAPQ Int TODO 
    #   It looks like there is no standard way to implement these so I'll hold
    #   off until my tool is a little more developed.
    output += str(255) + "\t"
    # 6 CIGAR String
    output += str(len(TEMPLATE)) + "M" + "\t"
    # 7 RNEXT String TODO This is the location of the "next sequence"
    output += str(0) + "\t"
    # 8 PNEXT Int TODO This is the location of the "next sequence"
    output += str(0) + "\t"
    # 9 TLEN Int
    output += str(len(TEMPLATE)) + "\t"
    # 10 SEQ String
    output += TEMPLATE + "\t"
    # 11 QUAL
    output += QUAL + "\t"

    return output + "\n"

def GET_METRICS(TOT_RDS_LEN:int, TOT_REF_LEN:int, 
                NUM_RDS:int, NUM_NA:float, #NUM_UA:float, 
                A_TIME:float, AF_TIME:float, REF_TIME:float, REFB_TIME:float) -> str:
    """
    Takes inputs of several recorded data figures to calculate important metrics and
    generate an output string for the metrics file.

    Parameters
    ----------
    NUM_RDS : int
        Number of reads that attempted to align
    TOT_RDS_LEN : int
        Total number of bp's of all reads
    TOT_REF_LEN : int
        Total number of bp's of reference genome.
    NUM_NA : float
        Num of reads Not Aligned
    NUM_UA : float
        Num of reads Uniquely Aligned
    A_TIME : float
        Time spent aligning reads and outputing to sam file
    AF_TIME : float
        Time spent actually finding alignments
    REF_TIME : float
        Time spent reading from reference file and building auxilury data structures
    REFB_TIME : float
        Time actually spent building auxilury data structures
    """

    metrics = ""

    # Header Line
    cur_time = time.strftime("%d:%H:%M")
    metrics += "Metrics file wf-align run at " + str(cur_time) + "\n\n"

    # Run time metrics

    metrics += "Runtime Metrics:\n"

    metrics += "Total Time:\t" + str(A_TIME + REF_TIME) + "\n"
    metrics += "IO Reference Time:\t" + str(REF_TIME - REFB_TIME) + "\n"
    metrics += "Building Auxiliary Structures Time:\t" + str(REFB_TIME) + "\n"

    metrics += "IO Reads Time:\t" + str(A_TIME - AF_TIME) + "\n"
    metrics += "Search Algorithm Time:\t" + str(AF_TIME) + "\n"

    metrics += "\n"

    # Accuracy metrics
    metrics += "Accuracy Metrics:\n"

    metrics += "Total Num Reads:\t" + str(NUM_RDS) + "\n"
    metrics += "Total Length of Reads:\t" + str(TOT_RDS_LEN) + "\n"
    metrics += "Total Length of Reference:\t" + str(TOT_REF_LEN) + "\n"

    #metrics += "Percent reads uniquely aligned:\t" + (NUM_UA)/NUM_RDS + "\n"
    #metrics += "Percent reads non-uniquely aligned:\t" + (NUM_RDS-NUM_NA-NUM_UA)/NUM_RDS + "\n"

    metrics += "Percent reads aligned:\t" + str(((NUM_RDS-NUM_NA)/NUM_RDS) if NUM_RDS > 0 else 0 )  + "\n"
    metrics += "Percent reads unaligned:\t" + str(((NUM_NA)/NUM_RDS) if NUM_RDS > 0 else 0) + "\n"

    return metrics
