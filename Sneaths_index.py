#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio import SeqIO
from Bio.Seq import Seq
import warnings
from Bio import BiopythonWarning
import numpy as np
from pprint import pprint


def calculate(
    indexD,
    taxonD,
    length,
    outgroupD
    ):
    """
    determine sneath's index for conserved positions
    
    Parameters
    ----------
    argv: indexD
        index of sneath's values for various AA combinations
    argv: taxonD
        taxon of interest and its sequence
    argv: length
        length of sequence
    argv: outgroupD
        outgroup taxa and their sequence
    """

    # initialize step and window size
    step   = 1
    window = 1

    # intialize list for np array of continuous identified snps
    summary_arr = [['start','stop','taxon','out','Sneath']]

    # suppress unnecessary warning
    warnings.simplefilter('ignore', BiopythonWarning)

    # loop through aligned sequence using step size
    for i in range(0, (int(length)+1) - int(step), int(step)):
        # initialize variables to hold sequence of interest being analyzed
        taxonSeg = ''
        outSeg   = ''
        # loop through all individuals from clade of interest and save 
        # window being analyzed into a concatenated string from all indivs
        for k, v in taxonD.items():
            taxonSeg += (v[i:i+int(window)])
        for k, v in outgroupD.items():
            outSeg   += (v[i:i+int(window)])

        # if outSeg AAs are all the same, get the outSeg character
        if len(set(outSeg)) == 1:
            outChar = outSeg[0:1]
            # find match bewteen outChar and taxonSeg in indexD
            for k, v in indexD.items():
                taxon_out_chars = []
                taxon_out_chars = [taxonSeg._data, outChar._data]
                # append to summary_arr if there is a match
                if set(k) == set(taxon_out_chars):
                    temp_list = []
                    # start
                    temp_list.append(i)
                    # stop
                    temp_list.append(i+int(window))
                    # taxonSeg
                    temp_list.append(taxonSeg._data)
                    # outgroup character
                    temp_list.append(outChar._data)
                    # Sneath's index value
                    temp_list.append(v)
                    # append to summary_arr
                    summary_arr.append(temp_list)

    # convert summary_arr to np array
    summary_arr = np.asarray(summary_arr)

    # if no snps are found, exit
    if len(summary_arr) == 0:
        sys.exit()

    # print np array line by line tab delimited
    for c0, c1, c2, c3, c4 in summary_arr:
        print("{}\t{}\t{}\t{}\t{}".format(c0, c1, c2, c3, c4))

def sneath(
    outgroupD,
    taxonD,
    length
    ):
    """
    loads sneath's index
    
    Parameters
    ----------
    argv: taxonD
        taxon of interest and its sequence
    argv: outgroupD
        outgroup taxa and their sequence
    argv: length
        length of sequence
    """

    # Create Sneath's index key and value pairs
    # https://en.wikipedia.org/wiki/Amino_acid_replacement#Sneath's_index
    indexD = {}
    # Isoleucine and leucine
    indexD[('I', 'L')] = 5
    # Valine and leucine
    indexD[('V', 'L')] = 9
    # Valine and isoleucine
    indexD[('V', 'I')] = 7
    # Glycine and leucine
    indexD[('G', 'L')] = 24
    # Glycine and isoleucine
    indexD[('G', 'I')] = 25
    # Glycine and valine
    indexD[('G', 'V')] = 19
    # Alanine and leucine
    indexD[('A', 'L')] = 15
    # Alanine and isoleucine
    indexD[('A', 'I')] = 17
    # Alanine and valine
    indexD[('A', 'V')] = 12
    # Alanine and glycine
    indexD[('A', 'G')] = 9
    # Proline and leucine
    indexD[('P', 'L')] = 23
    # Proline and isoleucine
    indexD[('P', 'I')] = 24
    # Proline and valine
    indexD[('P', 'V')] = 20
    # Proline and glycine
    indexD[('P', 'G')] = 17
    # Proline and alanine
    indexD[('P', 'A')] = 16
    # Glutamine and leucine
    indexD[('Q', 'L')] = 22
    # Glutamine and isoleucine
    indexD[('Q', 'I')] = 24
    # Glutamine and valine
    indexD[('Q', 'V')] = 25
    # Glutamine and glycine
    indexD[('Q', 'G')] = 32
    # Glutamine and alanine
    indexD[('Q', 'A')] = 26
    # Glutamine and proline
    indexD[('Q', 'P')] = 33
    # Asparagine and leucine
    indexD[('N', 'L')] = 20
    # Asparagine and isoleucine
    indexD[('N', 'I')] = 23
    # Asparagine and valine
    indexD[('N', 'V')] = 23
    # Asparagine and glycine
    indexD[('N', 'G')] = 26
    # Asparagine and alanine
    indexD[('N', 'A')] = 25
    # Asparagine and proline
    indexD[('N', 'P')] = 31
    # Asparagine and glutamine
    indexD[('N', 'Q')] = 10
    # Methionine and leucine
    indexD[('M', 'L')] = 20
    # Methionine and isoleucine
    indexD[('M', 'I')] = 22
    # Methionine and valine
    indexD[('M', 'V')] = 23
    # Methionine and glycine
    indexD[('M', 'G')] = 34
    # Methionine and alanine
    indexD[('M', 'A')] = 25
    # Methionine and proline
    indexD[('M', 'P')] = 31
    # Methionine and glutamine
    indexD[('M', 'Q')] = 13
    # Methionine and Asparagine
    indexD[('M', 'N')] = 21
    # Threonine and leucine
    indexD[('T', 'L')] = 23
    # Threonine and isoleucine
    indexD[('T', 'I')] = 21
    # Threonine and valine
    indexD[('T', 'V')] = 17
    # Threonine and glycine
    indexD[('T', 'G')] = 20
    # Threonine and alanine
    indexD[('T', 'A')] = 20
    # Threonine and proline
    indexD[('T', 'P')] = 25
    # Threonine and glutamine
    indexD[('T', 'Q')] = 24
    # Threonine and asparagine
    indexD[('T', 'N')] = 19
    # Threonine and methionine
    indexD[('T', 'M')] = 25
    # Serine and leucine
    indexD[('S', 'L')] = 23
    # Serine and isoleucine
    indexD[('S', 'I')] = 25
    # Serine and valine
    indexD[('S', 'V')] = 20
    # Serine and glycine
    indexD[('S', 'G')] = 19
    # Serine and alanine
    indexD[('S', 'A')] = 16
    # Serine and proline
    indexD[('S', 'P')] = 24
    # Serine and glutamine
    indexD[('S', 'Q')] = 21
    # Serine and asparagine
    indexD[('S', 'N')] = 15
    # Serine and methionine
    indexD[('S', 'M')] = 22
    # Serine and threonine
    indexD[('S', 'T')] = 12
    # Cysteine and leucine
    indexD[('C', 'L')] = 24
    # Cysteine and isoleucine
    indexD[('C', 'I')] = 26
    # Cysteine and valine
    indexD[('C', 'V')] = 21
    # Cysteine and glycine
    indexD[('C', 'G')] = 21
    # Cysteine and alanine
    indexD[('C', 'A')] = 13
    # Cysteine and proline
    indexD[('C', 'P')] = 25
    # Cysteine and glutamine
    indexD[('C', 'Q')] = 22
    # Cysteine and asparagine
    indexD[('C', 'N')] = 19
    # Cysteine and methionine
    indexD[('C', 'M')] = 17
    # Cysteine and threonine
    indexD[('C', 'T')] = 19
    # Cysteine and serine
    indexD[('C', 'S')] = 13
    # Glutamic acid and leucine
    indexD[('E', 'L')] = 30
    # Glutamic acid and isoleucine
    indexD[('E', 'I')] = 31
    # Glutamic acid and valine
    indexD[('E', 'V')] = 31
    # Glutamic acid and glycine
    indexD[('E', 'G')] = 37
    # Glutamic acid and alanine
    indexD[('E', 'A')] = 34
    # Glutamic acid and proline
    indexD[('E', 'P')] = 43
    # Glutamic acid and glutamine
    indexD[('E', 'Q')] = 14
    # Glutamic acid and asparagine
    indexD[('E', 'N')] = 19
    # Glutamic acid and methionine
    indexD[('E', 'M')] = 26
    # Glutamic acid and threonine
    indexD[('E', 'T')] = 34
    # Glutamic acid and serine
    indexD[('E', 'S')] = 29
    # Glutamic acid and cysteine
    indexD[('E', 'C')] = 33
    # Aspartic acid and 
    indexD[('D', 'L')] = 25
    # Aspartic acid and isoleucine
    indexD[('D', 'I')] = 28
    # Aspartic acid and valine
    indexD[('D', 'V')] = 28
    # Aspartic acid and glycine
    indexD[('D', 'G')] = 33
    # Aspartic acid and alanine 
    indexD[('D', 'A')] = 30
    # Aspartic acid and proline
    indexD[('D', 'P')] = 40
    # Aspartic acid and glutamine
    indexD[('D', 'Q')] = 22
    # Aspartic acid and asparagine
    indexD[('D', 'N')] = 14
    # Aspartic acid and methionine
    indexD[('D', 'M')] = 31
    # Aspartic acid and threonine
    indexD[('D', 'T')] = 29
    # Aspartic acid and serine
    indexD[('D', 'S')] = 25
    # Aspartic acid and cysteine
    indexD[('D', 'C')] = 28
    # Aspartic acid and glutamic acid
    indexD[('D', 'E')] = 7
    # Lysine and leucine
    indexD[('K', 'L')] = 23
    # Lysine and isoleucine
    indexD[('K', 'I')] = 24
    # Lysine and valine
    indexD[('K', 'V')] = 26
    # Lysine and glycine
    indexD[('K', 'G')] = 31
    # Lysine and alanine
    indexD[('K', 'A')] = 26
    # Lysine and proline
    indexD[('K', 'P')] = 31
    # Lysine and glutamine
    indexD[('K', 'Q')] = 21
    # Lysine and asparagine
    indexD[('K', 'N')] = 27
    # Lysine and methionine
    indexD[('K', 'M')] = 24
    # Lysine and threonine
    indexD[('K', 'T')] = 34
    # Lysine and serine
    indexD[('K', 'S')] = 31
    # Lysine and cysteine
    indexD[('K', 'C')] = 32
    # Lysine and glutamic acid
    indexD[('K', 'E')] = 26
    # Lysine and aspartic acid
    indexD[('K', 'D')] = 34
    # Arginine and leucine
    indexD[('R', 'L')] = 33
    # Arginine and isoleucine
    indexD[('R', 'I')] = 34
    # Arginine and valine
    indexD[('R', 'V')] = 36
    # Arginine and glycine
    indexD[('R', 'G')] = 43
    # Arginine and alanine
    indexD[('R', 'A')] = 37
    # Arginine and proline
    indexD[('R', 'P')] = 43
    # Arginine and glutamine
    indexD[('R', 'Q')] = 23
    # Arginine and asparagine
    indexD[('R', 'N')] = 31
    # Arginine and methionine
    indexD[('R', 'M')] = 28
    # Arginine and threonine
    indexD[('R', 'T')] = 38
    # Arginine and serine
    indexD[('R', 'S')] = 37
    # Arginine and cysteine
    indexD[('R', 'C')] = 36
    # Arginine and glutamic acid
    indexD[('R', 'E')] = 31
    # Arginine and aspartic acid
    indexD[('R', 'D')] = 39
    # Arginine and lysine
    indexD[('R', 'K')] = 14
    # Tyrosine and leucine
    indexD[('Y', 'L')] = 30
    # Tyrosine and isoleucine
    indexD[('Y', 'I')] = 34
    # Tyrosine and valine
    indexD[('Y', 'V')] = 36
    # Tyrosine and glycine
    indexD[('Y', 'G')] = 36
    # Tyrosine and alanine
    indexD[('Y', 'A')] = 34
    # Tyrosine and proline
    indexD[('Y', 'P')] = 37
    # Tyrosine and glutamine
    indexD[('Y', 'Q')] = 29
    # Tyrosine and asparagine
    indexD[('Y', 'N')] = 28
    # Tyrosine and methionine
    indexD[('Y', 'M')] = 32
    # Tyrosine and serine
    indexD[('Y', 'S')] = 29
    # Tyrosine and cysteine
    indexD[('Y', 'C')] = 34
    # Tyrosine and glutamic acid
    indexD[('Y', 'E')] = 34
    # Tyrosine and aspartic acid
    indexD[('Y', 'D')] = 34
    # Tyrosine and lysine
    indexD[('Y', 'K')] = 34
    # Tyrosine and arginine
    indexD[('Y', 'R')] = 36
    # Phenylalanine and leucine
    indexD[('F', 'L')] = 19
    # Phenylalanine and isoleucine
    indexD[('F', 'I')] = 22
    # Phenylalanine and valine
    indexD[('F', 'V')] = 26
    # Phenylalanine and glycine
    indexD[('F', 'G')] = 29
    # Phenylalanine and alanine
    indexD[('F', 'A')] = 26
    # Phenylalanine and proline
    indexD[('F', 'P')] = 27
    # Phenylalanine and glutamine
    indexD[('F', 'Q')] = 24
    # Phenylalanine and asparagine
    indexD[('F', 'N')] = 24
    # Phenylalanine and methionine
    indexD[('F', 'M')] = 24
    # Phenylalanine and threonine
    indexD[('F', 'T')] = 28
    # Phenylalanine and serine
    indexD[('F', 'S')] = 25
    # Phenylalanine and cysteine
    indexD[('F', 'C')] = 29
    # Phenylalanine and glutamic acid
    indexD[('F', 'E')] = 35
    # Phenylalanine and aspartic acid
    indexD[('F', 'D')] = 35
    # Phenylalanine and lysine
    indexD[('F', 'K')] = 28
    # Phenylalanine and arginine
    indexD[('F', 'R')] = 34
    # Phenylalanine and tyrosine
    indexD[('F', 'Y')] = 13 
    # Tryptophan and leucine
    indexD[('W', 'L')] = 30
    # Tryptophan and isoleucine
    indexD[('W', 'I')] = 34
    # Tryptophan and valine
    indexD[('W', 'V')] = 37
    # Tryptophan and glycine
    indexD[('W', 'G')] = 39
    # Tryptophan and alanine
    indexD[('W', 'A')] = 36
    # Tryptophan and proline
    indexD[('W', 'P')] = 37
    # Tryptophan and glutamine
    indexD[('W', 'Q')] = 31
    # Tryptophan and asparagine
    indexD[('W', 'N')] = 32
    # Tryptophan and methionine
    indexD[('W', 'M')] = 31
    # Tryptophan and threonine
    indexD[('W', 'T')] = 38
    # Tryptophan and serine
    indexD[('W', 'S')] = 35
    # Tryptophan and cysteine
    indexD[('W', 'C')] = 37
    # Tryptophan and glutamic acid
    indexD[('W', 'E')] = 43
    # Tryptophan and aspartic acid
    indexD[('W', 'E')] = 45
    # Tryptophan and lysine
    indexD[('W', 'K')] = 34
    # Tryptophan and arginine
    indexD[('W', 'R')] = 36
    # Tryptophan and tyrosine
    indexD[('W', 'Y')] = 21
    # Tryptophan and phenylalanine
    indexD[('W', 'F')] =  13
    # Histidine and leucine
    indexD[('H', 'L')] = 25
    # Histidine and isoleucine
    indexD[('H', 'I')] = 28
    # Histidine and valine
    indexD[('H', 'V')] = 31
    # Histidine and glycine
    indexD[('H', 'G')] = 34
    # Histidine and alanine
    indexD[('H', 'A')] = 29
    # Histidine and proline
    indexD[('H', 'P')] = 36
    # Histidine and glutamine
    indexD[('H', 'Q')] = 27
    # Histidine and asparagine
    indexD[('H', 'N')] = 24
    # Histidine and methionine
    indexD[('H', 'M')] = 30
    # Histidine and threonine
    indexD[('H', 'T')] = 34
    # Histidine and serine
    indexD[('H', 'S')] = 28
    # Histidine and cysteine
    indexD[('H', 'C')] = 31
    # Histidine and glutamic acid
    indexD[('H', 'E')] = 27
    # Histidine and aspartic acid
    indexD[('H', 'D')] = 35
    # Histidine and lysine
    indexD[('H', 'K')] = 27
    # Histidine and arginine
    indexD[('H', 'N')] = 31
    # Histidine and tyrosine
    indexD[('H', 'Y')] = 23
    # Histidine and phenylalanine
    indexD[('H', 'F')] = 18
    # Histidine and typtophan
    indexD[('H', 'W')] = 25

    calculate(
        indexD,
        taxonD,
        length,
        outgroupD
        )

def read_seq_in(
    taxonL, 
    fasta, 
    outgroupL,
    ):
    """
    Reads sequences in
    
    Parameters
    ----------
    argv: taxonL
        list with taxon of interest
    argv: fasta
        input fasta file
    argv: outgroupL
        list of outgroup taxa
    """

    # initialize step and window
    step   = 1
    window = 1

    # initialize taxa dictionaries
    outgroupD = {}
    taxonD     = {}

    # initialize length variable
    length = ''
    
    # read in aligned fasta file
    format    = "fasta"
    handle    = open(fasta)
    alignment = list(SeqIO.parse(handle, format))

    # fill outgroupD dictionary with individual id (key) and sequence (value)
    for indiv in alignment:
        if indiv.id in outgroupL:
            outgroupD[indiv.id] = indiv.seq
    # fill taxonD dictionary with individual id (key) and sequence (value)
    for indiv in alignment:
        if indiv.id in taxonL:
            taxonD[indiv.id] = indiv.seq

    # determine sequence length
    length    = len(alignment[0].seq)

    sneath(
        outgroupD,
        taxonD,
        length
        )
    
def read_clades_into_list(
    taxon, 
    fasta,
    outgroup
    ):
    """
    reads taxon into lists
    
    Parameters
    ----------
    argv: taxon
        string of taxon name of interest
    argv: fasta
        file should be an aligned fasta file
    argv: outgroup
        file with taxa names from the outgroup
        taxa names should be in a single column
    """
    
    # initialize clade lists
    taxonL    = []
    outgroupL = []

    # read clade files into lists
    taxonL    = [taxon]
    outgroupL = [line.rstrip('\n') for line in open(outgroup)]

    # pass to sneath function
    read_seq_in(
        taxonL, 
        fasta, 
        outgroupL,
        )

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    fasta = ''

    try:
        opts, args = getopt.getopt(argv, "ht:i:g:")
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
            print("\nScans the outgroup (-g parameter) for conserved amino acids and then determines")
            print("if there is an amino acid substitution in the taxon of interest (-t parameter).")
            print("If there is a substitution, the Sneath's index value of dissimilarity is reported")
            print("Values for Sneath's index values can be found here: ")
            print("https://en.wikipedia.org/wiki/Amino_acid_replacement#Sneath's_index\n")
            print("Output is 5 columns that specify the position of the substition (cols 1 and 2),")
            print("the taxon of interests sequence (col 3), the outgroup's sequence (col 4), and the")
            print("Sneath's index at that position (col 5).\n\n")
            ## options
            # taxon of interest explanation
            print("\n-t\t<taxon of interest>:")
            print("\ta single taxon name present in the fasta file")
            # fasta file
            print("\n-i\t<input fasta file>:")
            print("\tInput fasta file.")
            # clade outgroup file explanation
            print("\n-g\t<clade outgroup of interest>:")
            print("\ta single columns file with the fasta header names from clade outgroup")
            print("\tof interest. This group should be the one used to contextualize the")
            print("\tAA substitution identified between the taxon of interest and the outgroup\n\n")
            sys.exit()

        elif opt == '-t':
            if arg:
                taxon = arg
                #print("\nFile with taxa from clade 1: {}".format(clade1))
            else:
                # error message
                print("\n\nPlease enter a name for the taxon of interest.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                fasta = arg
                #print("\nAlignment input fasta file: {}".format(alignedFasta))
            else:
                # error message
                print("\n\nThe specified file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-g':
            if os.path.isfile(arg):
                outgroup = arg
                #print("\cladeOutgroup input file: {}".format(cladeOutgroup))
            else:
                # error message
                print("\n\nThe specified file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to sneath function
    read_clades_into_list(
        taxon,
        fasta,
        outgroup
        )

if __name__ == '__main__':
    main(sys.argv[1:])
