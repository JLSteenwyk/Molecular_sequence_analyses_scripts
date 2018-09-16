#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio import SeqIO
import numpy as np

def find_conserved_positions(
    taxaL, alignedFasta
    ):
    """
    Scans for positions with conserved characters
    
    Parameters
    ----------
    argv: taxa
        file with taxa names from taxa of interest
        taxa names should be in a single column
    argv: alignedFasta
        input fasta file
    argv: character
        character argument that specifies nucl or prot
    """

    # initialize step and window
    step   = 1
    window = 1

    # initialize taxa dictionaries
    taxaD  = {}

    # initialize length variable
    length = ''
    
    # read in aligned fasta file
    format    = "fasta"
    handle    = open(alignedFasta)
    alignment = list(SeqIO.parse(handle, format))

    # fill taxaD dictionary with individual id (key) and sequence (value)
    for indiv in alignment:
        if indiv.id in taxaL:
            taxaD[indiv.id] = indiv.seq

    # determine sequence length
    length    = len(alignment[0].seq)

    # intialize list for np array of continuous identified indels
    snp_arr = [['start','stop','char','pos']]

    # loop through aligned sequence using step size
    for i in range(0, (int(length)+1) - int(step), int(step)):
        # initialize variables to hold sequence of interest being analyzed
        taxaseg = ''

        # loop through all individuals from taxa of interest and save 
        # window being analyzed into a concatenated string from all indivs
        for k, v in taxaD.items():
            taxaseg   += (v[i:i+int(window)])
            # make all characters the same case
            taxaseg=taxaseg.upper()
        taxaseg=set(taxaseg)
       
        # test if taxaseg only contains 1 character or more
        taxasegChar=','.join(str(s) for s in taxaseg)
        if (len(taxaseg) < 2) and (taxasegChar!="-"):
            # append pertinent information to indel_arr
            temp_list = []
            # start
            temp_list.append(i)
            # stop
            temp_list.append(i+int(window))
            # char
            temp_list.append(taxasegChar)
            # determine position of char in codon
            if (i+int(window))%3==0:
                temp_list.append("3rd")
            elif (i+int(window))%3==1:
                temp_list.append("1st")
            elif (i+int(window))%3==2:
                temp_list.append("2nd")
            snp_arr.append(temp_list)

    # convert snp_arr to np array
    snp_arr = np.asarray(snp_arr)

    # if no snps are found, exit
    if len(snp_arr) == 0:
        sys.exit()

    # print np array line by line tab delimited
    for c0, c1, c2, c3 in snp_arr:
        print("{}\t{}\t{}\t{}".format(c0, c1, c2, c3))

def read_clade_into_list(
    taxa, alignedFasta
    ):
    """
    Reads taxa of interest in to a list 
    
    Parameters
    ----------
    argv: taxa
        file with taxa names from taxa of interest
        taxa names should be in a single column
    argv: alignedFasta
        input fasta file
    """
    
    # initialize taxa list
    taxaL = []


    # read clade files into lists
    taxaL = [line.rstrip('\n') for line in open(taxa)]

    # determine number of indivs in each clade
    lenTaxaL = len(taxaL)

    # pass to find_conserved_positions function
    find_conserved_positions(
        taxaL, alignedFasta
        )

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    taxa         = ''
    alignedFasta = ''
    character    = ''

    try:
        opts, args = getopt.getopt(argv, "ht:i:")
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
            # script explanation
            print("\nThis script identifies positions in a fasta alignment file (-i parameter)")
            print("where for the taxa of interest (-t parameter), the script will determine the")
            print("positions that have the same character. For example, if position 1 in all taxa")
            print("have an 'A', then the script will report that position 1 has the same nucleotide")
            print("in all taxa.")
            # taxa file explanation
            print("\n-t\ttaxa of interest:")
            print("\ta single column file with the fasta header names for taxa")
            print("\tof interest. If the names do not match the fasta header, they")
            print("\twill not be recognized and accounted for.")
            # aligned fasta file explanation
            print("\n-i\taligned fasta file:")
            print("\taligned nucleotide fasta file containing (but not limited to)")
            print("\ttaxa specified in clade 1 (-o) and clade 2 (-t). If taxa are in")
            print("\tthe fasta file not specified in the clade 1 and clade 2 files,")
            print("\tthe taxa will not be considered in the analysis.")
            # output explanation
            print("\noutput explanation:")
            print("\tThe output will have 4 columns titled:")
            print("\tstart, stop, char, pos")
            print("\tcol1, start, represents where the conserved character of interest begins")
            print("\tcol2, stop, represents where the conserved character of interest ends")
            print("\tcol3, char, will be the nucleotide present in the character of interest")
            print("\tcol4, pos, will be the codon position of the character of interest\n")
            sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                taxa = arg
                #print("\nFile with taxa from clade 1: {}".format(clade1))
            else:
                # error message
                print("\n\nThe specified taxa file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                alignedFasta = arg
                #print("\nAlignment input fasta file: {}".format(alignedFasta))
            else:
                # error message
                print("\n\nThe specified fasta file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to read_clade_into_list function
    read_clade_into_list(
        taxa, alignedFasta
        )

if __name__ == '__main__':
    main(sys.argv[1:])