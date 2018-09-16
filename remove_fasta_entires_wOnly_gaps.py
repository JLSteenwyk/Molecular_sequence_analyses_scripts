#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import SeqIO

def parse(
    fastaFile
    ):
    """
    Scans for positions with conserved characters
    
    Parameters
    ----------
    argv: fastaFile
        fasta file with n number of sequences with only gaps
    """

    # read in fasta file
    format    = "fasta"
    handle    = open(fastaFile)
    fasta     = list(SeqIO.parse(handle, format))

    # fill taxaD dictionary with individual id (key) and sequence (value)
    for indiv in fasta:
        if set(indiv.seq) == set(["-"]):
            continue
        elif set(indiv.seq) != set(["-"]):
            print("{}{}\n{}".format(">", indiv.id, indiv.seq))

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    fastaFile = ''

    try:
        opts, args = getopt.getopt(argv, "hi:")
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
            print("\nThis script removes entries in mfasta files that are only filled with gaps.")
            # fasta file explanation
            print("\n-i\tfasta file")
            print("\tfasta file that contains n number of sequences that have only gaps")
            # output explanation
            print("\noutput explanation:")
            print("\tThe output is a fasta file with n number of sequences that have only gaps removed.\n")
            sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                fastaFile = arg
                #print("\nFile with taxa from clade 1: {}".format(clade1))
            else:
                # error message
                print("\n\nThe specified fasta file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to parse function
    parse(
        fastaFile
        )

if __name__ == '__main__':
    main(sys.argv[1:])