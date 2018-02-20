#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import numpy as np
from datetime import datetime

def shuffle(
    fasta, sequence_type
    ):
    """
    Reads the alignment list and taxa list
    
    Parameters
    ----------
    argv: fasta
        input fasta file
    argv: sequence_type
        specifies sequence type
    """

    # read in fasta file to alignment variable
    format    = "fasta"
    handle    = open(fasta)
    alignment = list(SeqIO.parse(handle, format))

    # initialize step size
    step = 1

    # create a variable with all taxa names
    taxa_list = []
    for record in alignment:
        taxa_list.append(record.id)

    # initialize dictionary to hold sequences for each taxa
    seqDict = {}
    # populate seqDict with individual id (key) and sequence (value)
    for indiv in alignment:
        if indiv.id in taxa_list:
            seqDict[indiv.id] = indiv.seq

    # determine length of sequence
    length    = len(alignment[0].seq)

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

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # start time of script
    startTime = datetime.now()

    # initialize argument variables
    alignment_list  = ''
    taxa_list       = ''
    prefix          = ''

    try:
        opts, args = getopt.getopt(argv, "hi:t:")
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
            ## explanation
            print("\nShuffle the nucleotides or amino acids at each site in a fasta file.")
            print("Because this script shuffles, the relative number of each nucleotide or")
            print("amino acid is maintained after shuffling. This script is specifically")
            print("designed to facilitate randomization tests. For example, the tree length")
            print("test for clonality -- see figure 4 in http://jcm.asm.org/content/41/2/703.full")
            ## options
            # fasta file
            print("\n-i <input fasta file>")
            print("\tSingle column file of the alignment files that will be concatenated.")
            # specify if nucleotide or protein files list
            print("\n-t <fasta files are either nucleotide or protein fastas>")
            print("\tArgument can be either 'prot' or 'nucl' to specify if the")
            print("\talignments contain protein or nucleotide sequences.\n")
            sys.exit()

        elif opt == '-i':
            if os.path.isfile(arg):
                fasta = arg
            else:
                # error message
                print("\n\nThe specified fasta file (-i) does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-t':
            if arg in ('prot', 'nucl'):
                sequence_type = arg 
            else:
                # error message
                print("\n\nThe specified input for (-t) should be 'prot' or 'nucl'.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()

    # pass to shuffle function
    shuffle(
        fasta, sequence_type
        )

if __name__ == '__main__':
    main(sys.argv[1:])