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
    Scans for small snps using a sliding window
    
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
    cladeD   = {}
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
            cladeD[indiv.id]   = indiv.seq
        # fill cladeOutD dictionary with individual id (key) and sequence (value)
        elif indiv.id in cladeOutL:
            cladeOutD[indiv.id] = indiv.seq

    # determine sequence length
    length    = len(alignment[0].seq)

    # intialize list for np array of continuous identified indels
    snp_arr = [['pos','C1','OG']]

    # loop through aligned sequence using step size
    for i in range(0, (int(length)+1) - int(step), int(step)):
        # initialize variables to hold sequence of interest being analyzed
        cladeSeg   = ''
        cladeOutseg = ''
        # loop through all individuals from clade of interest and save 
        # window being analyzed into a concatenated string from all indivs
        for k, v in cladeD.items():
            cladeSeg   += (v[i:i+int(window)])
        for k, v in cladeOutD.items():
            cladeOutseg += (v[i:i+int(window)])

        # If cladeSeg only contains 1 unique character
        if len(set(cladeSeg)) == 1:
        	# if the unique character is not found in cladeOutseg then print position
            if cladeSeg[0] not in set(cladeOutseg):
                # append pertinent information to snp_arr
                temp_list = []
                # position
                temp_list.append(i+1)
                # clade 
                temp_list.append(cladeSeg[0])
                # out
                temp_list.append(set(cladeOutseg))
                # append
                snp_arr.append(temp_list)

    # convert snp_arr to np array
    snp_arr = np.asarray(snp_arr)

    # if no snps are found, exit
    if len(snp_arr) == 0:
        sys.exit()

    # print np array line by line tab delimited
    for c0, c1, c2 in snp_arr:
        print("{}\t{}\t{}".format(c0, c1, ''.join(c2)))


def read_clades_into_list(
    window, clade, 
    alignedFasta,
    step, cladeOutgroup
    ):
    """
    Scans for small indels 
    
    Parameters
    ----------
    argv: window
        window size to go through fasta file
    argv: clade
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
    cladeL   = []
    cladeOutL = []

    # read clade files into lists
    cladeL   = [line.rstrip('\n') for line in open(clade)]
    cladeOutL = [line.rstrip('\n') for line in open(cladeOutgroup)]

    # determine number of indivs in each clade
    lenClade = len(cladeL)
    lenCladeOut = len(cladeOutL)

    # pass to find_indels function
    find_snps(
        window, cladeL, 
        alignedFasta, 
        step, lenClade, 
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
    clade         = ''
    alignedFasta  = ''
    step          = 1
    cladeOutgroup = ''

    try:
        opts, args = getopt.getopt(argv, "hw:o:i:s:g:")
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
            print("\n-o\tclade 1 of interest:")
            print("\ta single column file with the fasta header names from clade 1")
            print("\tof interest. If the names do not match the fasta header, they")
            print("\twill not be recognized and accounted for.")
            # clade outgroup file explanation
            print("\n-g\tclade outgroup of interest:")
            print("\ta single columns file with the fasta header names from clade outgroup")
            print("\tof interest. This group should be the one used to contextualize the")
            print("\tSNP identified between clade 1. For example, if As are")
            print("\tobserved in clade 1 and Ts and clade outgroup this script will report a")
            print("\tT->A mutation in clade 1. In the same instance,")
            print("\tif the outgroup has Ns or -s, it will be reported as ambiguous.")
            # aligned fasta file explanation
            print("\n-i\taligned fasta file:")
            print("\taligned fasta file containing (but not limited to)")
            print("\ttaxa specified in clade 1 (-o) and the outgroup (-g). If taxa are in")
            print("\tthe fasta file not specified in the clade 1 and outgroup files,")
            print("\tthe taxa will not be considered in the analysis.")
            # output explanation
            print("\noutput explanation:")
            print("\tThe output will have 6 columns titled:")
            print("\tstart, stop, C1, OG, class, C1type, C2type")
            print("\t  - start refers to the starting position of the SNP and stop is the end")
            print("\t  - C1 refers to the nucleotide found in clade 1 specified using -o")
            print("\t  - OG refers to the nucleotide found in outgroup clade specified using -g")
            print("\t  - class will be populated with either a C1, C2, or C1C2. This refers to")
            print("\t\twhich clade has the SNP. C1 refers to clade1, C2")
            print("\t\trefers to clade2 and C1C2 refers to both clades.")
            sys.exit()
        elif opt == '-o':
            if os.path.isfile(arg):
                clade = arg
                #print("\nFile with taxa from clade 1: {}".format(clade1))
            else:
                # error message
                print("\n\nThe specified file does not exist.\n")
                print("For detailed explanation use -h argument\n")
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
        window, clade, 
        alignedFasta,
        step, cladeOutgroup
        )

if __name__ == '__main__':
    main(sys.argv[1:])