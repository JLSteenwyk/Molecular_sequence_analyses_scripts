#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import AlignIO
from pprint import pprint


def modFasta(
    fasta
    ):
    """
    prints to stdout the modified fasta format file

    Parameters
    ----------
    argv: fasta
        alignment fasta file
    """
    
    # convert
    alignment=AlignIO.read(open(fasta), "fasta")
    num    = len([1 for line in open(fasta) if line.startswith(">")])
    length = len(alignment[0]._seq)
    
    print(num,"",length)
    for i in range(0,len(alignment)):
        print(alignment[i].name,"",alignment[i]._seq)
    #pprint(vars(alignment[0]))


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
            # general script explanation
            print("\n This script takes an input fasta file and outputs a paml fasta file.")
            print("")
            print("Required arguments include:")
            print("-i: input fasta file")
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
    
    # pass to modFasta function
    modFasta(
        fasta
        )

if __name__ == '__main__':
    main(sys.argv[1:])