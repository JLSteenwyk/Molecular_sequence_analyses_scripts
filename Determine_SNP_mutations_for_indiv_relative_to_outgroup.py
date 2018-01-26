#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio import SeqIO
import numpy as np

def find_snps(
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

    # DNA single letters and gap set
    A    = set('Aa')
    T    = set('Tt')
    C    = set('Cc')
    G    = set('Gg')
    N    = set('Nn')
    gap  = set('-')

    butA = set('TtCcGgNn-')
    butT = set('AaCcGgNn-')
    butC = set('AaTtGgNn-')
    butG = set('AaTtCcNn-')

    # intialize list for np array of continuous identified snps
    snp_arr = [['start','stop','C1','OG','C1type']]

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
       
        ## outgroup As
        if set(str(cladeOutseg)) <= A:
            # clade1 Ts
            if set(str(clade1seg)) <= T:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # outgroup
                temp_list.append("A")
                # C1type
                temp_list.append("Tv")
                snp_arr.append(temp_list)
            # clade1 Cs
            elif set(str(clade1seg)) <= C:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # outgroup
                temp_list.append("A")
                # C1type
                temp_list.append("Tv")
                snp_arr.append(temp_list)
            # clade1 Gs
            elif set(str(clade1seg)) <= G:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # outgroup
                temp_list.append("A")
                # C1type
                temp_list.append("Ts")
                snp_arr.append(temp_list)
        ## outgroup Ts
        elif set(str(cladeOutseg)) <= T:
            # clade1 Ts
            if set(str(clade1seg)) <= A:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # outgroup
                temp_list.append("T")
                # C1type
                temp_list.append("Tv")
                snp_arr.append(temp_list)
            # clade1 Cs
            elif set(str(clade1seg)) <= C:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # outgroup
                temp_list.append("T")
                # C1type
                temp_list.append("Ts")
                snp_arr.append(temp_list)
            # clade1 Gs
            elif set(str(clade1seg)) <= G:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # outgroup
                temp_list.append("T")
                # C1type
                temp_list.append("Tv")
                snp_arr.append(temp_list)
        ## outgroup Cs
        elif set(str(cladeOutseg)) <= C:
            # clade1 Ts
            if set(str(clade1seg)) <= A:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # outgroup
                temp_list.append("C")
                # C1type
                temp_list.append("Tv")
                snp_arr.append(temp_list)
            # clade1 Cs
            elif set(str(clade1seg)) <= T:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # outgroup
                temp_list.append("C")
                # C1type
                temp_list.append("Ts")
                snp_arr.append(temp_list)
            # clade1 Gs
            elif set(str(clade1seg)) <= G:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # outgroup
                temp_list.append("C")
                # C1type
                temp_list.append("Tv")
                snp_arr.append(temp_list)
        ## outgroup Gs
        elif set(str(cladeOutseg)) <= G:
            # clade1 Ts
            if set(str(clade1seg)) <= A:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # outgroup
                temp_list.append("G")
                # C1type
                temp_list.append("Ts")
                snp_arr.append(temp_list)
            # clade1 Cs
            elif set(str(clade1seg)) <= T:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # outgroup
                temp_list.append("G")
                # C1type
                temp_list.append("Tv")
                snp_arr.append(temp_list)
            # clade1 Gs
            elif set(str(clade1seg)) <= C:
                # append pertinent information to snp_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # outgroup
                temp_list.append("G")
                # C1type
                temp_list.append("Tv")
                snp_arr.append(temp_list)


    #print(snp_arr)
    # convert snp_arr to np array
    snp_arr = np.asarray(snp_arr)

    # if no snps are found, exit
    if len(snp_arr) == 0:
        sys.exit()

    # print np array line by line tab delimited
    for c0, c1, c2, c3, c4 in snp_arr:
        print("{}\t{}\t{}\t{}\t{}".format(c0, c1, c2, c3, c4))

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
    find_snps(
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
    window        = 1
    clade1        = ''
    clade2        = ''
    alignedFasta  = ''
    step          = 1
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
            # clade 1 file explanation
            print("\n-o\tindividual of interest:")
            print("\tstring of individual of interest name in fasta file for which")
            print("\tsnps will be determined in.")
            # clade outgroup file explanation
            print("\n-g\tclade outgroup of interest:")
            print("\ta single columns file with the fasta header names from clade outgroup")
            print("\tof interest. This group should be the one used to contextualize the")
            print("\tSNP identified between clade 1 and clade 2. For example, if As are")
            print("\tobserved in clade 1 and Ts are observed in clade 2 and clade outgroup")
            print("\tthis script will report a T->A mutation in clade 1. In the same instance,")
            print("\tif the outgroup has Ns or -s, it will be reported as ambiguous.")
            # aligned fasta file explanation
            print("\n-i\taligned fasta file:")
            print("\taligned nucleotide fasta file containing (but not limited to)")
            print("\ttaxa specified in clade 1 (-o) and clade 2 (-t). If taxa are in")
            print("\tthe fasta file not specified in the clade 1 and clade 2 files,")
            print("\tthe taxa will not be considered in the analysis.")
            # output explanation
            print("\noutput explanation:")
            print("\tThe output will have 6 columns titled:")
            print("\tstart, stop, C1, OG, C1type")
            print("\t  - start refers to the starting position of the SNP and stop is the end")
            print("\t  - C1 refers to the nucleotide found in clade 1 specified using -o")
            print("\t  - OG refers to the nucleotide found in outgroup clade specified using -g")
            print("\t  - C1type refers to if the SNP found at C1 is a Transition (Ts) or")
            print("\t\tor Transversion (Tv). If no SNP is found, the column will contain NA")
            sys.exit()
        elif opt == '-o':
            if arg:
                clade1 = arg
                #print("\nIndividual of interest: {}".format(clade1))
            else:
                # error message
                print("\n\nInvalid input provided.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
                sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                alignedFasta = arg
                #print("\nAlignment input fasta file: {}".format(alignedFasta))
            else:
                # error message
                print("\n\nThe specified file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-g':
            if os.path.isfile(arg):
                cladeOutgroup = arg
                #print("\cladeOutgroup input file: {}".format(cladeOutgroup))
            else:
                # error message
                print("\n\nThe specified file does not exist.\n")
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