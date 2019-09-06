#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import AlignIO


def modFasta(
    fasta, formatOut
    ):
    """
    prints to stdout the modified fasta format file

    Parameters
    ----------
    argv: fasta
        alignment fasta file
    argv: formatOut
        output format
    """
    
    # convert
    alignment=AlignIO.read(open(fasta), "fasta")
    print(alignment.format(formatOut))


def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    fasta     = ''
    formatOut = ''

    try:
        opts, args = getopt.getopt(argv, "hi:f:")
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
            # general script explanation
            print("\n This script takes an input fasta file and outputs a file in a different format to stdout.")
            print("Choices of output include:")
            print("clustal, emboss, fasta, fasta-m10, ig, maf, mauve, nexus, phylip, phylip-sequential, phylip-relaxed, and stockholm.\n")
            print("Required arguments include:")
            print("-i: input fasta file")
            print("-f: output format")
            print("\n")
            sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                fasta = arg
            else:
                # error message
                print("\n\nThe specified fasta file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()   
        elif opt == '-f':
            if arg:
                formatOut = arg
            else:
                # error message
                print("\n\nFormat of output file\n")
                sys.exit()
    
    # pass to modFasta function
    modFasta(
        fasta, formatOut
        )

if __name__ == '__main__':
    main(sys.argv[1:])