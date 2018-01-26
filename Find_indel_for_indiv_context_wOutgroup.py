#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio import SeqIO
import numpy as np

def find_indels(
    window, clade1L, 
    alignedFasta,
    step, lenClade1, 
    cladeOutL,
    lenCladeOut
    ):
    """
    Scans for small indels using a sliding window
    Collapses adjacent windows at the end and prints
    out results.
    
    Parameters
    ----------
    argv: window
        window size to go through fasta file
    argv: clade1L
        list with taxa names from clade 1
        taxa names should be in a single column
    argv: alignedFasta
        file should be an aligned fasta file
    argv: lenClade1
        number of indivs in clade 1
    argv: cladeOutL
        list with taxa names from clade outgroup
        taxa names should be in a single column
    argv: lenCladeOut
        number of indivs in clade outgroupls
    """

    # initialize clade dictionaries
    clade1D   = {}
    cladeOutD = {}
    # initialize length variable
    length    = ''
    
    # read in aligned fasta file
    format    = "fasta"
    handle    = open(alignedFasta)
    alignment = list(SeqIO.parse(handle, format))

    # fill clade1D dictionary with individual id (key) and sequence (value)
    for indiv in alignment:
        if indiv.id in clade1L:
            clade1D[indiv.id]   = indiv.seq
        # fill cladeOutD dictionary with individual id (key) and sequence (value)
        elif indiv.id in cladeOutL:
            cladeOutD[indiv.id] = indiv.seq

    # determine sequence length
    length    = len(alignment[0].seq)

    # DNA acid single letters and gap set
    nts       = set('ATCGNatcgn')
    gap       = set('-')

    # intialize list for np array of continuous identified indels
    indel_arr = []

    # loop through aligned sequence using step size
    for i in range(0, (int(length)+1) - int(step), int(step)):
        # initialize variables to hold sequence of interest being analyzed
        clade1seg   = ''
        cladeOutseg = ''
        # loop through all individuals from clade of interest and save 
        # window being analyzed into a concatenated string from all indivs
        for k, v in clade1D.items():
            clade1seg   += (v[i:i+int(window)])
        for k, v in cladeOutD.items():
            cladeOutseg += (v[i:i+int(window)])
       
        # test if clade1seg only contains nts and cladeOutseg only contains gaps
        if (set(str(clade1seg)) <= nts) and (set(str(cladeOutseg)) <= gap):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "gap2", str(clade1seg), str(cladeOutseg))

            # add pertinent information to indel_arr
            temp_list = []
            temp_list.append(i)
            temp_list.append(i+int(window))
            temp_list.append("outGap")
            temp_list.append("ins")
            indel_arr.append(temp_list)

        # test if clade1seg only contains gaps and clade2seg only contains aas
        elif (set(str(clade1seg)) <= gap) and (set(str(cladeOutseg)) <= nts):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "gap1", str(clade1seg), str(clade2seg), str(cladeOutseg))

            # add pertinent information to indel_arr
            temp_list = []
            temp_list.append(i)
            temp_list.append(i+int(window))
            temp_list.append("outNts")
            temp_list.append("del")
            indel_arr.append(temp_list)


    # convert indel_arr to np array
    indel_arr = np.asarray(indel_arr)
    # if no indels are found, exit
    if len(indel_arr) == 0:
        sys.exit()

    ## Looping over lines to collapse adjacently identified indels
    merged_data = []
    idx_ii      = 0
    idx_jj      = idx_ii + 1
    idx_opt     = True
    while_opt   = True
    # Looping over first line
    while(idx_opt):
        try:
            line_ii = indel_arr[idx_ii]
            line_jj = indel_arr[idx_jj]
            # Testing condition
            line_ii_opt = True
            # loop through indel_arr
            while(line_ii_opt):
                new_line = line_ii
                # if end line 1 (ii) is greater or equal to line 2 (jj) and
                # line 1 clade with gap is equal to line 2 clade with gap and
                # segment identity of outgroup is equal in line 1 and line 2
                if (line_ii[1] >= line_jj[0]) and (line_ii[3]==line_jj[3]) and (line_ii[2]==line_jj[2]):
                    idx_jj  += int(1)
                    # new line is line 1 start, line 2 stop, line 2 outgroup ID, line 2 gap ID
                    new_line = [line_ii[0], line_jj[1], line_jj[2], line_jj[3]]
                    line_ii  = new_line
                    line_jj  = indel_arr[idx_jj]
                elif (line_ii[1] >= line_jj[0]) and (line_ii[3]==line_jj[3]) and (line_ii[2]!=line_jj[2]):
                    idx_jj  += int(1)
                    new_line = [line_ii[0], line_jj[1], "outCombo", line_jj[3]]
                    line_ii  = new_line
                    line_jj  = indel_arr[idx_jj]
                else:
                    line_ii_opt = False
                    merged_data.append(new_line)
                    idx_ii = idx_jj
                    idx_jj = idx_ii + 1
        # except for end of indel_arr
        except IndexError:
            merged_data.append(line_ii)
            idx_opt = False
    
    # Convert to np array
    merged_data = np.asarray(merged_data)

    # print np array line by line tab delimited
    for c0, c1, c2, c3 in merged_data:
        size = int(c1)-int(c0)  
        print("{}\t{}\t{}\t{}\t{}\t".format(c0, c1, size, c2, c3))


def read_clades_into_list(
    window, clade1, 
    alignedFasta,
    step, cladeOutgroup
    ):
    """
    Scans for small indels 
    
    Parameters
    ----------
    argv: window
        window size to go through fasta file
    argv: clade1
        file with taxa names from clade 1
        taxa names should be in a single column
    argv: alignedFasta
        file should be an aligned fasta file
    argv: step
        step of sliding window
    argv: cladeOutgroup
        file with taxa anmes from clade outgroup
        taxa names should be in a single column
    """
    
    # initialize clade lists
    clade1L   = []
    cladeOutL = []

    # read clade files into lists
    clade1L   = clade1
    cladeOutL = [line.rstrip('\n') for line in open(cladeOutgroup)]

    # determine number of indivs in each clade
    lenClade1 = len(clade1L)
    lenCladeOut = len(cladeOutL)

    # pass to find_indels function
    find_indels(
        window, clade1L, 
        alignedFasta, 
        step, lenClade1, 
        cladeOutL,
        lenCladeOut
        )

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    window        = ''
    clade1        = ''
    alignedFasta  = ''
    step          = ''
    cladeOutgroup = ''

    try:
        opts, args = getopt.getopt(argv, "hw:o:t:i:s:g:")
    except getopt.GetoptError:
        # error message
        print("Error\nFor help use -h argument\n")
        sys.exit(2)
    # if no arguments are used print a help message
    if len(opts) == 0:
        # error message
        print("\nNo arguments provided...")
        print("For help use -h argument\n")
        sys.exit(2)
    # test for arguments
    for opt, arg in opts:
        if opt == '-h':
            # window size explanation
            print("\n-w\twindow size:")
            print("\tspecifies the number of nucleotides to be analyzing at once")
            # clade 1 file explanation
            print("\n-o\tindividual of interest:")
            print("\ta single column file with the fasta header names from clade 1")
            print("\tof interest. If the names do not match the fasta header, they")
            print("\twill not be recognized and accounted for.")
            # clade outgroup file explanation
            print("\n-g\tclade outgroup of interest:")
            print("\ta single columns file with the fasta header names from clade outgroup")
            print("\tof interest. This group should be the one used to contextualize the")
            print("\tindels identified between clade 1 and clade 2. For example, if gaps are")
            print("\tobserved in clade 1 and nts are observed in clade 2 and clade outgroup")
            print("\tthis script will report a deletion in clade 1. In the same instance, the")
            print("\toutgroup has both gaps and nts, it will be reported as ambiguous.")
            # aligned fasta file explanation
            print("\n-i\taligned fasta file:")
            print("\taligned nucleotide fasta file containing (but not limited to)")
            print("\ttaxa specified in clade 1 (-o) and clade 2 (-t). If taxa are in")
            print("\tthe fasta file not specified in the clade 1 and clade 2 files,")
            print("\tthe taxa will not be considered in the analysis.")
            # step size of analysis
            print("\n-s\tstep size:")
            print("\tsize of step as window moves through aligned fasta file\n")
            sys.exit()
        elif opt == '-w':
            if float(arg).is_integer():
                window = arg
                #print("\nWindow size is {}".format(window))
            else:
                # error message
                print("\n\nProvided window size is not an integer.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-o':
            if arg:
                clade1 = arg
                #print("\nFile with taxa from clade 1: {}".format(clade1))
            else:
                # error message
                print("\n\nInvalid input for argument -o.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                alignedFasta = arg
                #print("\nAlignment input fasta file: {}".format(alignedFasta))
            else:
                # error message
                print("\n\nThe specified file does not exist.\n")
                print("Check -i argument\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-g':
            if os.path.isfile(arg):
                cladeOutgroup = arg
                #print("\cladeOutgroup input file: {}".format(cladeOutgroup))
            else:
                # error message
                print("\n\nThe specified file does not exist.\n")
                print("Check -g argument\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-s':
            if float(arg).is_integer():
                step = arg
                #print("Step: {}".format(step))
            else:
                # error message
                print("\n\nProvided window size is not an integer.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to read_clades_into_list function
    read_clades_into_list(
        window, clade1, 
        alignedFasta,
        step, cladeOutgroup
        )

if __name__ == '__main__':
    main(sys.argv[1:])
