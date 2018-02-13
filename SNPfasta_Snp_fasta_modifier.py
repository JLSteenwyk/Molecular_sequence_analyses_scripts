#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import SeqIO
import numpy as np

def replace_nucl(
    FASTA, ENTRY,
    COORD, OLD, NEW
    ):
    """
    replaces nucleotide with SNP
    
    Parameters
    ----------
    argv: FASTA
        fasta file
    argv: ENTRY
        fasta entry ID
    argv: COORD
        coordinate to modify
    argv: OLD
        old nucleotide
    argv: NEW
        new nucleotide
    """

    # read in fasta file
    FAdict = {}
    format = "fasta"
    handle = open(FASTA)
    FAdict = SeqIO.to_dict(SeqIO.parse(handle, format))

    # intialize ID and seq holder
    ID     = ''
    Seq    = ''
    NewSeq = ''

    # populate ID and Seq
    ID  = FAdict[ENTRY].id
    Seq = FAdict[ENTRY].seq

    # check if coordinate matches the old nucleotide and print the new fasta entry if it does
    COORD=int(COORD)-1
    if Seq[COORD] == OLD:
        # create string with new sequence
        NewSeq=Seq[:COORD] + NEW + Seq[COORD + 1:]
        print(">"+ID+"_NToriginal")
        print(Seq)
        print(">"+ID+"_NTmod")
        print(NewSeq)
        print(">"+ID+"_NToriginal")
        print(Seq.translate())        
        print(">"+ID+"_PROTmod")
        print(NewSeq.translate())
    else:
        print("Position", COORD+1, "is a", Seq[COORD], "and not a", OLD)
        print("Check that the correct position has been specified")
        sys.exit()

    # reassign seqs with the amino acid sequence
    Seq     = Seq.translate()
    NewSeq  = NewSeq.translate()
    
    # result lists
    res_arr   = []
    res_arr.append(["Pos", "Old", "New"])
    # summarize differences between the two AA strings
    for AA in range(0, (int(len(Seq))), 1):
        SeqAA    = Seq[AA:AA+1]
        NewSeqAA = NewSeq[AA:AA+1]
        #if strings are not equal
        if SeqAA != NewSeqAA:
            temp_arr = []
            temp_arr.append(AA+1)
            temp_arr.append(SeqAA[0])
            temp_arr.append(NewSeqAA[0])
            res_arr.append(temp_arr)


    # check if res_arr has no new additions and exit if so
    if res_arr==[["Pos", "Old", "New"]]:
        sys.exit()

    # convert to numpy array
    res_arr = np.asarray(res_arr)

    # print np array line by line tab delimited
    for c0, c1, c2 in res_arr:
        print("{}\t{}\t{}".format(c0, c1, c2))


def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    FASTA  = ''
    ENTRY  = ''
    COORD  = ''
    OLD    = ''
    NEW    = ''

    try:
        opts, args = getopt.getopt(argv, "hf:e:c:o:n:")
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
            print("\nModify a fasta file according to a user specified SNP.")
            print("The resulting fasta entry will contain the fasta entry modified to have a mutation.")
            # fasta file explanation
            print("\n-f\tprotein alignment:")
            print("\tSingle or multiple fasta file.")
            # fasta entry
            print("\n-e\tentry from fasta file:")
            print("\tTo accommodate multi-fasta files, the specific entry in the fasta file must be")
            print("\tspecified. This could be a specific gene or scaffold.")
            # starting coordinate
            print("\n-c\tmutation coordinate:")
            print("\tThe coordinate in the fasta file that will be mutagenized. For example, if the")
            print("\tuser puts a 5, the 5th position in the fasta entry will be mutagenized.")
            # the old nucleotide/aa
            print("\n-o\told nucl/aa:")
            print("\tAs a quality control check, the user must specify what the old nucl/aa was.")
            print("\tThe script will check if the coordinate you specified is the nucl/aa the user")
            print("\tthinks it is.")
            # the new nucleotide/aa
            print("\n-n\tnew nucl/aa:")
            print("\tSpecify what the old nucl/aa now is. For example, if the mutation was an A->T")
            print("\ttransversion, specify T here.\n")
            sys.exit()
        # read in fasta file name
        elif opt == '-f':
            if os.path.isfile(arg):
                FASTA = arg
            else:
                # error message
                print("\n\nThe specified fasta file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        # read in entry ID
        elif opt == '-e':
            ENTRY = arg
        # read in coordinate
        elif opt == '-c':
            if int(arg):
                COORD = arg
            else:
                # error message
                print("\n\nThe specified coordinate is not an integer.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        # read in old nucl/aa
        elif opt == '-o':
            if arg:
                OLD = arg
            else:
                # error message
                print("\n\nPlease specify a nucl/aa for -o argument\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        # read in old nucl/aa
        elif opt == '-n':
            if arg:
                NEW = arg
            else:
                # error message
                print("\n\nPlease specify a nucl/aa for -n argument\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to thread_dna function
    replace_nucl(
        FASTA, ENTRY,
        COORD, OLD, NEW
        )

if __name__ == '__main__':
    main(sys.argv[1:])