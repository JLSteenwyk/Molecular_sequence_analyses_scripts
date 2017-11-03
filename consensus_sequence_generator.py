#!/usr/bin/env python

import sys
import getopt
import os.path
import os
from Bio import AlignIO
from Bio.Align import AlignInfo

def create_consensus(
	input_file, threshold_in
    ):
    """
    Reads in alignment file and creates a consensus sequences from
    sequences in the alignment file using BioPython

    Parameters
    ----------
    argv: input_file
        MSA file
    argv: threshold_in
        residue threshold
    """

    # read in MSA file
    with open(input_file, "r") as f:
        # open fasta file
        format_MSA = "fasta"
        handle = open(input_file)
        alignment = AlignIO.read(handle, format_MSA)
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus(threshold = float(threshold_in), ambiguous='X')
    header = ">" + str(input_file) + ".consensus"
    print("{}\n{}".format(header,consensus))




def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    input_file    = ''
    threshold_in  = ''

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
            # description
            print("\n\nThis script will take a multiple sequence alignment (MSA) and create")
            print("a consensus sequence from the MSA. This is done using BioPython's")
            print("dumb_consensus utility. Output will print a sequence with the header")
            print("'>' + 'input_file' + '.consensus'")
            print("Written by Jacob Steenwyk, November 2 2017")
            # options
            print("\n-i <input alignment file>")
            print("\tfile should be a multiple sequence alignment (MSA) file")
            print("\tredirect output to a new file using '>'\n")
            print("\n-t <threshold>")
            print("\tthis threshold specifies how common a particular residue has to be at a position")
            print("\tbefore it is added. Value should be between 0 and 1 reflecting presence in")
            print("\t0-100 percent in the alignment.\n\n")
            sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                input_file = arg
            else:
                # error message
                print("\n\nThe specified input alignment file does not exist.\n")
                print("For detailed explanation of usage supply the -h argument\n")
                sys.exit()
        elif opt == '-t':
            if 0 <= float(arg) <= 1:
                threshold_in = arg
            else:
                # error message
                print("\n\nThe specified value is not between 0 and 1.\n")
                print("For detailed explanation of usage supply the -h argument\n")
                sys.exit()

    # pass to read_config parameter
    create_consensus(
        input_file, threshold_in
        )

if __name__ == '__main__':
    main(sys.argv[1:])