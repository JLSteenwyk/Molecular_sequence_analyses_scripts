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
    summary_arr = [['start','stop','taxon','out','Epstein']]

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
                taxon_out_chars = (taxonSeg._data, outChar._data)
                # append to summary_arr if there is a match
                if k == taxon_out_chars:
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

def epstein(
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

    # Create Epstein's coefficients of difference key and value pairs
    # https://en.wikipedia.org/wiki/Amino_acid_replacement#Epstein's_coefficient_of_difference
    indexD = {}
    # F to M
    indexD[('F', 'M')] = 0.1
    # F to L
    indexD[('F', 'L')] = 0.15
    # F to I
    indexD[('F', 'I')] = 0.15
    # F to V
    indexD[('F', 'V')] = 0.2
    # F to P
    indexD[('F', 'P')] = 0.2
    # F to Y
    indexD[('F', 'Y')] = 0.2
    # F to W
    indexD[('F', 'W')] = 0.21
    # F to C
    indexD[('F', 'C')] = 0.28
    # F to A
    indexD[('F', 'A')] = 0.5
    # F to G
    indexD[('F', 'G')] = 0.61
    # F to S
    indexD[('F', 'S')] = 0.81
    # F to T
    indexD[('F', 'T')] = 0.81
    # F to H
    indexD[('F', 'H')] = 0.8
    # F to E
    indexD[('F', 'E')] = 1
    # F to Q
    indexD[('F', 'Q')] = 1
    # F to D
    indexD[('F', 'D')] = 1
    # F to N
    indexD[('F', 'N')] = 1
    # F to K
    indexD[('F', 'K')] = 1
    # F to R
    indexD[('F', 'R')] = 1
    # M to F
    indexD[('M', 'F')] = 0.05
    # M to L
    indexD[('M', 'L')] = 0.05
    # M to I
    indexD[('M', 'I')] = 0.05
    # M to V
    indexD[('M', 'V')] = 0.1
    # M to P
    indexD[('M', 'P')] = 0.1
    # M to Y
    indexD[('M', 'Y')] = 0.22
    # M to W
    indexD[('M', 'W')] = 0.24
    # M to C
    indexD[('M', 'C')] = 0.22
    # M to A
    indexD[('M', 'A')] = 0.45
    # M to G
    indexD[('M', 'G')] = 0.56
    # M to S
    indexD[('M', 'S')] = 0.8
    # M to T
    indexD[('M', 'T')] = 0.8
    # M to H
    indexD[('M', 'H')] = 0.8
    # M to E
    indexD[('M', 'E')] = 1
    # M to Q
    indexD[('M', 'Q')] = 1
    # M to D
    indexD[('M', 'D')] = 1
    # M to N
    indexD[('M', 'N')] = 1
    # M to K
    indexD[('M', 'K')] = 1
    # M to R
    indexD[('M', 'R')] = 1
    # L to F
    indexD[('L', 'F')] = 0.08
    # L to M
    indexD[('L', 'M')] = 0.03
    # L to I
    indexD[('L', 'I')] = 0
    # L to V
    indexD[('L', 'V')] = 0.05
    # L to P
    indexD[('L', 'P')] = 0.05
    # L to Y
    indexD[('L', 'Y')] = 0.22
    # L to W
    indexD[('L', 'W')] = 0.25
    # L to C
    indexD[('L', 'C')] = 0.21
    # L to A
    indexD[('L', 'A')] = 0.43
    # L to G
    indexD[('L', 'G')] = 0.54
    # L to S
    indexD[('L', 'S')] = 0.8
    # L to T
    indexD[('L', 'T')] = 0.8
    # L to E
    indexD[('L', 'E')] = 1
    # L to Q
    indexD[('L', 'Q')] = 1
    # L to D
    indexD[('L', 'D')] = 1
    # L to N
    indexD[('L', 'N')] = 1
    # L to K
    indexD[('L', 'K')] = 1
    # L to R
    indexD[('L', 'R')] = 1
    # I to F
    indexD[('I', 'F')] = 0.08
    # I to M
    indexD[('I', 'M')] = 0.03
    # I to L
    indexD[('I', 'L')] = 0
    # I to V
    indexD[('I', 'V')] = 0.05  
    # I to P
    indexD[('I', 'P')] = 0.05
    # I to Y
    indexD[('I', 'Y')] = 0.22
    # I to W
    indexD[('I', 'W')] = 0.25
    # I to C
    indexD[('I', 'C')] = 0.21
    # I to A
    indexD[('I', 'A')] = 0.43
    # I to G
    indexD[('I', 'G')] = 0.54
    # I to S
    indexD[('I', 'S')] = 0.8
    # I to T
    indexD[('I', 'T')] = 0.8
    # I to H
    indexD[('I', 'H')] = 1
    # I to E
    indexD[('I', 'E')] = 1
    # I to Q
    indexD[('I', 'Q')] = 1
    # I to D
    indexD[('I', 'D')] = 1
    # I to N
    indexD[('I', 'N')] = 1
    # I to K
    indexD[('I', 'K')] = 1
    # I to R
    indexD[('I', 'R')] = 1
    # V to F
    indexD[('V', 'F')] = 0.1
    # V to M
    indexD[('V', 'M')] = 0.1
    # V to L
    indexD[('V', 'L')] = 0.03
    # V to I
    indexD[('V', 'I')] = 0.03
    # V to P
    indexD[('V', 'P')] = 0
    # V to Y
    indexD[('V', 'Y')] = 0.24
    # V to W
    indexD[('V', 'W')] = 0.27
    # V to C
    indexD[('V', 'C')] = 0.2
    # V to A
    indexD[('V', 'A')] = 0.41
    # V to G
    indexD[('V', 'G')] = 0.52
    # V to S
    indexD[('V', 'S')] = 0.8
    # V to T
    indexD[('V', 'T')] = 0.8
    # V to H
    indexD[('V', 'H')] = 0.8
    # V to E
    indexD[('V', 'E')] = 1
    # V to Q
    indexD[('V', 'Q')] = 1
    # V to D
    indexD[('V', 'D')] = 1
    # V to N
    indexD[('V', 'N')] = 1
    # V to K
    indexD[('V', 'K')] = 1
    # V to R
    indexD[('V', 'R')] = 1.01
    # P to F
    indexD[('P', 'F')] = 0.1
    # P to M
    indexD[('P', 'M')] = 0.1
    # P to L
    indexD[('P', 'L')] = 0.03
    # P to I 
    indexD[('P', 'I')] = 0.03
    # P to V
    indexD[('P', 'V')] = 0
    # P to Y
    indexD[('P', 'Y')] = 0.24
    # P to W
    indexD[('P', 'W')] = 0.27
    # P to C
    indexD[('P', 'C')] = 0.2
    # P to A
    indexD[('P', 'A')] = 0.41
    # P to G
    indexD[('P', 'G')] = 0.52
    # P to S
    indexD[('P', 'S')] = 0.8
    # P to T
    indexD[('P', 'T')] = 0.8
    # P to H
    indexD[('P', 'H')] = 0.8
    # P to E
    indexD[('P', 'E')] = 1
    # P to Q
    indexD[('P', 'Q')] = 1
    # P to D
    indexD[('P', 'D')] = 1
    # P to N
    indexD[('P', 'N')] = 1
    # P to K
    indexD[('P', 'K')] = 1
    # P to R
    indexD[('P', 'R')] = 1.01
    # Y to F
    indexD[('Y', 'F')] = 0.21
    # Y to M
    indexD[('Y', 'M')] = 0.25
    # Y to L
    indexD[('Y', 'L')] = 0.28
    # Y to I
    indexD[('Y', 'I')] = 0.28
    # Y to V
    indexD[('Y', 'V')] = 0.32
    # Y to P
    indexD[('Y', 'P')] = 0.32
    # Y to W
    indexD[('Y', 'W')] = 0.05
    # Y to C
    indexD[('Y', 'C')] = 0.25
    # Y to A
    indexD[('Y', 'A')] = 0.4
    # Y to G
    indexD[('Y', 'G')] = 0.5
    # Y to S
    indexD[('Y', 'S')] = 0.62
    # Y to T
    indexD[('Y', 'T')] = 0.61
    # Y to H
    indexD[('Y', 'H')] = 0.6
    # Y to E
    indexD[('Y', 'E')] = 0.8
    # Y to Q
    indexD[('Y', 'Q')] = 0.8
    # Y to D
    indexD[('Y', 'D')] = 0.81
    # Y to N
    indexD[('Y', 'N')] = 0.81
    # Y to K
    indexD[('Y', 'K')] = 0.8
    # Y to R
    indexD[('Y', 'R')] = 0.8
    # W to F
    indexD[('W', 'F')] = 0.25
    # W to M
    indexD[('W', 'M')] = 0.32
    # W to L
    indexD[('W', 'L')] = 0.36
    # W to I
    indexD[('W', 'I')] = 0.36
    # W to V
    indexD[('W', 'V')] = 0.4
    # W to P
    indexD[('W', 'P')] = 0.4
    # W to Y
    indexD[('W', 'Y')] = 0.1
    # W to C
    indexD[('W', 'C')] = 0.35
    # W to A
    indexD[('W', 'A')] = 0.49
    # W to G
    indexD[('W', 'G')] = 0.58
    # W to S
    indexD[('W', 'S')] = 0.63
    # W to T
    indexD[('W', 'T')] = 0.63
    # W to H
    indexD[('W', 'H')] = 0.61
    # W to E
    indexD[('W', 'E')] = 0.81
    # W to Q
    indexD[('W', 'Q')] = 0.81
    # W to D
    indexD[('W', 'D')] = 0.81
    # W to N
    indexD[('W', 'N')] = 0.81
    # W to K
    indexD[('W', 'K')] = 0.81
    # W to R
    indexD[('W', 'R')] = 0.8
    # C to F 
    indexD[('C', 'F')] = 0.22
    # C to M
    indexD[('C', 'M')] = 0.21
    # C to L
    indexD[('C', 'L')] = 0.2
    # C to I
    indexD[('C', 'I')] = 0.2
    # C to V
    indexD[('C', 'V')] = 0.2
    # C to P
    indexD[('C', 'P')] = 0.2
    # C to Y
    indexD[('C', 'Y')] = 0.13
    # C to W
    indexD[('C', 'W')] = 0.18
    # C to A
    indexD[('C', 'A')] = 0.22
    # C to G
    indexD[('C', 'G')] = 0.34
    # C to S
    indexD[('C', 'S')] = 0.6
    # C to T
    indexD[('C', 'T')] = 0.6
    # C to H
    indexD[('C', 'H')] = 0.61
    # C to E
    indexD[('C', 'E')] = 0.8
    # C to Q
    indexD[('C', 'Q')] = 0.8
    # C to D
    indexD[('C', 'D')] = 0.8
    # C to N
    indexD[('C', 'N')] = 0.8
    # C to K
    indexD[('C', 'K')] = 0.8
    # C to R
    indexD[('C', 'R')] = 0.81
    # A to F
    indexD[('A', 'F')] = 0.43
    # A to M
    indexD[('A', 'M')] = 0.41
    # A to L
    indexD[('A', 'L')] = 0.43
    # A to I
    indexD[('A', 'I')] = 0.43
    # A to V
    indexD[('A', 'V')] = 0.4
    # A to P
    indexD[('A', 'P')] = 0.4
    # A to Y
    indexD[('A', 'Y')] = 0.27
    # A to W
    indexD[('A', 'W')] = 0.3
    # A to C
    indexD[('A', 'C')] = 0.25
    # A to G
    indexD[('A', 'G')] = 0.1
    # A to S
    indexD[('A', 'S')] = 0.4
    # A to T
    indexD[('A', 'T')] = 0.4
    # A to H
    indexD[('A', 'H')] = 0.42
    # A to E
    indexD[('A', 'E')] = 0.61
    # A to Q
    indexD[('A', 'Q')] = 0.61
    # A to D
    indexD[('A', 'D')] = 0.61
    # A to N
    indexD[('A', 'N')] = 0.61
    # A to K
    indexD[('A', 'K')] = 0.61
    # A to R
    indexD[('A', 'R')] = 0.62
    # G to F
    indexD[('G', 'F')] = 0.53
    # G to M
    indexD[('G', 'M')] = 0.42
    # G to L
    indexD[('G', 'L')] = 0.51
    # G to I
    indexD[('G', 'I')] = 0.51
    # G to V
    indexD[('G', 'V')] = 0.5
    # G to P
    indexD[('G', 'P')] = 0.5
    # G to Y
    indexD[('G', 'Y')] = 0.36
    # G to W
    indexD[('G', 'W')] = 0.39
    # G to C
    indexD[('G', 'C')] = 0.31
    # G to A
    indexD[('G', 'A')] = 0.1
    # G to S
    indexD[('G', 'S')] = 0.3
    # G to T
    indexD[('G', 'T')] = 0.31
    # G to H
    indexD[('G', 'H')] = 0.34
    # G to E
    indexD[('G', 'E')] = 0.52
    # G to Q
    indexD[('G', 'Q')] = 0.52
    # G to D
    indexD[('G', 'D')] = 0.51
    # G to N
    indexD[('G', 'N')] = 0.51
    # G to K
    indexD[('G', 'K')] = 0.52
    # G to R
    indexD[('G', 'R')] = 0.53
    # S to F
    indexD[('S', 'F')] = 0.81
    # S to M
    indexD[('S', 'M')] = 0.8
    # S to L
    indexD[('S', 'L')] = 0.8
    # S to I
    indexD[('S', 'I')] = 0.8
    # S to V
    indexD[('S', 'V')] = 0.8
    # S to P
    indexD[('S', 'P')] = 0.8
    # S to Y
    indexD[('S', 'Y')] = 0.62
    # S to W
    indexD[('S', 'W')] = 0.63
    # S to C
    indexD[('S', 'C')] = 0.6
    # S to A
    indexD[('S', 'A')] = 0.4
    # S to G
    indexD[('S', 'G')] = 0.32
    # S to T
    indexD[('S', 'T')] = 0.03
    # S to H
    indexD[('S', 'H')] = 0.1
    # S to E
    indexD[('S', 'E')] = 0.22
    # S to Q
    indexD[('S', 'Q')] = 0.22
    # S to D
    indexD[('S', 'D')] = 0.21
    # S to N
    indexD[('S', 'N')] = 0.21
    # S to K
    indexD[('S', 'K')] = 0.22
    # S to R
    indexD[('S', 'R')] = 0.24
    # T to F
    indexD[('T', 'F')] = 0.81
    # T to M
    indexD[('T', 'M')] = 0.8
    # T to L
    indexD[('T', 'L')] = 0.8
    # T to I
    indexD[('T', 'I')] = 0.8
    # T to V
    indexD[('T', 'V')] = 0.8
    # T to P
    indexD[('T', 'P')] = 0.8
    # T to Y
    indexD[('T', 'Y')] = 0.61
    # T to W
    indexD[('T', 'W')] = 0.63
    # T to C
    indexD[('T', 'C')] = 0.6
    # T to A
    indexD[('T', 'A')] = 0.41
    # T to G
    indexD[('T', 'G')] = 0.34
    # T to S
    indexD[('T', 'S')] = 0.03
    # T to H
    indexD[('T', 'H')] = 0.08
    # T to E
    indexD[('T', 'E')] = 0.21
    # T to Q
    indexD[('T', 'Q')] = 0.21
    # T to D
    indexD[('T', 'D')] = 0.2
    # T to N
    indexD[('T', 'N')] = 0.2
    # T to K
    indexD[('T', 'K')] = 0.21
    # T to R
    indexD[('T', 'R')] = 0.22
    # H to F
    indexD[('H', 'F')] = 0.8
    # H to M
    indexD[('H', 'M')] = 0.8
    # H to L
    indexD[('H', 'L')] = 0.81
    # H to I
    indexD[('H', 'I')] = 0.81
    # H to V
    indexD[('H', 'V')] = 0.81
    # H to P
    indexD[('H', 'P')] = 0.81
    # H to Y
    indexD[('H', 'Y')] = 0.6
    # H to W
    indexD[('H', 'W')] = 0.61
    # H to C
    indexD[('H', 'C')] = 0.62
    # H to A
    indexD[('H', 'A')] = 0.47
    # H to G
    indexD[('H', 'G')] = 0.42
    # H to S
    indexD[('H', 'S')] = 0.1
    # H to T
    indexD[('H', 'T')] = 0.08
    # H to E
    indexD[('H', 'E')] = 0.2
    # H to Q
    indexD[('H', 'Q')] = 0.2
    # H to D
    indexD[('H', 'D')] = 0.21
    # H to N
    indexD[('H', 'N')] = 0.21
    # H to K
    indexD[('H', 'K')] = 0.2
    # H to R
    indexD[('H', 'R')] = 0.2
    # E to F
    indexD[('E', 'F')] = 1
    # E to M 
    indexD[('E', 'M')] = 1
    # E to L
    indexD[('E', 'L')] = 1
    # E to I
    indexD[('E', 'I')] = 1
    # E to V
    indexD[('E', 'V')] = 1
    # E to P
    indexD[('E', 'P')] = 1
    # E to Y
    indexD[('E', 'Y')] = 0.8
    # E to W
    indexD[('E', 'W')] = 0.81
    # E to C
    indexD[('E', 'C')] = 0.81
    # E to A
    indexD[('E', 'A')] = 0.63
    # E to G
    indexD[('E', 'G')] = 0.56
    # E to S
    indexD[('E', 'S')] = 0.21
    # E to T
    indexD[('E', 'T')] = 0.21
    # E to H
    indexD[('E', 'H')] = 0.2
    # E to Q
    indexD[('E', 'Q')] = 0
    # E to D
    indexD[('E', 'D')] = 0.03
    # E to N
    indexD[('E', 'N')] = 0.03
    # E to K
    indexD[('E', 'K')] = 0
    # E to R
    indexD[('E', 'R')] = 0.05
    # Q to F
    indexD[('Q', 'F')] = 1
    # Q to M
    indexD[('Q', 'M')] = 1
    # Q to L
    indexD[('Q', 'L')] = 1
    # Q to I
    indexD[('Q', 'I')] = 1
    # Q to V
    indexD[('Q', 'V')] = 1
    # Q to P
    indexD[('Q', 'P')] = 1
    # Q to Y
    indexD[('Q', 'Y')] = 0.8
    # Q to W
    indexD[('Q', 'W')] = 0.81
    # Q to C
    indexD[('Q', 'C')] = 0.81
    # Q to A
    indexD[('Q', 'A')] = 0.63
    # Q to G
    indexD[('Q', 'G')] = 0.56
    # Q to S
    indexD[('Q', 'S')] = 0.21
    # Q to T
    indexD[('Q', 'T')] = 0.21
    # Q to H
    indexD[('Q', 'H')] = 0.2
    # Q to E
    indexD[('Q', 'E')] = 0
    # Q to D
    indexD[('Q', 'D')] = 0.03
    # Q to N
    indexD[('Q', 'N')] = 0.03
    # Q to K
    indexD[('Q', 'K')] = 0
    # Q to R
    indexD[('Q', 'R')] = 0.05
    # D to F
    indexD[('D', 'F')] = 1
    # D to M
    indexD[('D', 'M')] = 1
    # D to L
    indexD[('D', 'L')] = 1
    # D to I
    indexD[('D', 'I')] = 1
    # D to V
    indexD[('D', 'V')] = 1
    # D to P
    indexD[('D', 'P')] = 1
    # D to Y
    indexD[('D', 'Y')] = 0.81
    # D to W
    indexD[('D', 'W')] = 0.81
    # D to C
    indexD[('D', 'C')] = 0.8
    # D to A
    indexD[('D', 'A')] = 0.62
    # D to G
    indexD[('D', 'G')] = 0.54
    # D to S
    indexD[('D', 'S')] = 0.2
    # D to T
    indexD[('D', 'T')] = 0.2
    # D to H
    indexD[('D', 'H')] = 0.21
    # D to E
    indexD[('D', 'E')] = 0.03
    # D to Q
    indexD[('D', 'Q')] = 0.03
    # D to N
    indexD[('D', 'N')] = 0
    # D to K
    indexD[('D', 'K')] = 0.03
    # D to R
    indexD[('D', 'R')] = 0.08
    # N to F
    indexD[('N', 'F')] = 1
    # N to M
    indexD[('N', 'M')] = 1
    # N to L
    indexD[('N', 'L')] = 1
    # N to I
    indexD[('N', 'I')] = 1
    # N to V
    indexD[('N', 'V')] = 1
    # N to P
    indexD[('N', 'P')] = 1
    # N to Y
    indexD[('N', 'Y')] = 0.81
    # N to W
    indexD[('N', 'W')] = 0.81
    # N to C
    indexD[('N', 'C')] = 0.8
    # N to A
    indexD[('N', 'A')] = 0.62
    # N to G
    indexD[('N', 'G')] = 0.54
    # N to S
    indexD[('N', 'S')] = 0.2
    # N to T
    indexD[('N', 'T')] = 0.2
    # N to H
    indexD[('N', 'H')] = 0.21
    # N to E
    indexD[('N', 'E')] = 0.03
    # N to Q
    indexD[('N', 'Q')] = 0.03
    # N to D
    indexD[('N', 'D')] = 0
    # N to K
    indexD[('N', 'K')] = 0.03
    # N to R
    indexD[('N', 'R')] = 0.08
    # K to F
    indexD[('K', 'F')] = 1
    # K to M
    indexD[('K', 'M')] = 1
    # K to L
    indexD[('K', 'L')] = 1
    # K to I
    indexD[('K', 'I')] = 1
    # K to V
    indexD[('K', 'V')] = 1
    # K to P
    indexD[('K', 'P')] = 1
    # K to Y
    indexD[('K', 'Y')] = 0.8
    # K to W
    indexD[('K', 'W')] = 0.81
    # K to C
    indexD[('K', 'C')] = 0.81
    # K to A
    indexD[('K', 'A')] = 0.63
    # K to G
    indexD[('K', 'G')] = 0.56
    # K to S
    indexD[('K', 'S')] = 0.21
    # K to T
    indexD[('K', 'T')] = 0.21
    # K to H
    indexD[('K', 'H')] = 0.2
    # K to E
    indexD[('K', 'E')] = 0
    # K to Q
    indexD[('K', 'Q')] = 0
    # K to D
    indexD[('K', 'D')] = 0.03
    # K to N
    indexD[('K', 'N')] = 0.03
    # K to R
    indexD[('K', 'R')] = 0.0
    # R to F
    indexD[('R', 'F')] = 1
    # R to M
    indexD[('R', 'M')] = 1
    # R to L
    indexD[('R', 'L')] = 1.01
    # R to I
    indexD[('R', 'I')] = 1.01
    # R to V
    indexD[('R', 'V')] = 1.02
    # R to P
    indexD[('R', 'P')] = 1.02
    # R to Y
    indexD[('R', 'Y')] = 0.8
    # R to W
    indexD[('R', 'W')] = 0.8
    # R to C
    indexD[('R', 'C')] = 0.82
    # R to A
    indexD[('R', 'A')] = 0.67
    # R to G
    indexD[('R', 'G')] = 0.61
    # R to S
    indexD[('R', 'S')] = 0.24
    # R to T
    indexD[('R', 'T')] = 0.22
    # R to H
    indexD[('R', 'H')] = 0.2
    # R to E
    indexD[('R', 'E')] = 0.05
    # R to Q
    indexD[('R', 'Q')] = 0.05
    # R to D
    indexD[('R', 'D')] = 0.08
    # R to N
    indexD[('R', 'N')] = 0.08
    # R to K
    indexD[('R', 'K')] = 0.05

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

    # initialize taxa dictionaries
    outgroupD  = {}
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

    epstein(
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
            print("If there is a substitution, the Epsteins's coefficient of difference is determined.")
            print("Values for Epstein's coefficient of difference values can be found here: ")
            print("https://en.wikipedia.org/wiki/Amino_acid_replacement#Epstein's_coefficient_of_difference\n")
            print("Output is 5 columns that specify the position of the substition (cols 1 and 2),")
            print("the taxon of interests sequence (col 3), the outgroup's sequence (col 4), and the")
            print("Epsteins's index at that position (col 5).\n\n")
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

    # pass to read clades function
    read_clades_into_list(
        taxon,
        fasta,
        outgroup
        )

if __name__ == '__main__':
    main(sys.argv[1:])
