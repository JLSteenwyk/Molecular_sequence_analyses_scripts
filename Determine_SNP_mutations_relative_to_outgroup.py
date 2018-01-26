#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio import SeqIO
import numpy as np

def find_snps(
    window, clade1L, 
    clade2L, alignedFasta,
    step, lenClade1, 
    lenClade2, cladeOutL,
    lenCladeOut
    ):
    """
    Scans for small indels using a sliding window
    Collapses adjacent windows at the end and prints
    out results.
    
    Parameters
    ----------
    argv: window
        window size to go through fasta file
    argv: clade1L
        list with taxa names from clade 1
        taxa names should be in a single column
    argv: clade2L
        list with taxa names from clade 2
        taxa names should be in a single column
    argv: alignedFasta
        file should be an aligned fasta file
    argv: lenClade1
        number of indivs in clade 1
    argv: lenClade2
        number of indivs in clade 2
    argv: cladeOutL
        list with taxa names from clade outgroup
        taxa names should be in a single column
    argv: lenCladeOut
        number of indivs in clade outgroupls
    """

    # initialize clade dictionaries
    clade1D   = {}
    clade2D   = {}
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
            clade1D[indiv.id]   = indiv.seq
        # fill clade2D dictionary with individual id (key) and sequence (value)
        elif indiv.id in clade2L:
            clade2D[indiv.id]   = indiv.seq
        # fill cladeOutD dictionary with individual id (key) and sequence (value)
        elif indiv.id in cladeOutL:
            cladeOutD[indiv.id] = indiv.seq

    # determine sequence length
    length    = len(alignment[0].seq)

    # DNA single letters and gap set
    A    = set('Aa')
    T    = set('Tt')
    C    = set('Cc')
    G    = set('Gg')
    N    = set('Nn')
    gap  = set('-')

    butA = set('TtCcGgNn-')
    butT = set('AaCcGgNn-')
    butC = set('AaTtGgNn-')
    butG = set('AaTtCcNn-')

    # intialize list for np array of continuous identified indels
    snp_arr = [['start','stop','C1','C2','OG','class','C1type','C2type','pos']]

    # loop through aligned sequence using step size
    for i in range(0, (int(length)+1) - int(step), int(step)):
        # initialize variables to hold sequence of interest being analyzed
        clade1seg   = ''
        clade2seg   = ''
        cladeOutseg = ''
        # loop through all individuals from clade of interest and save 
        # window being analyzed into a concatenated string from all indivs
        for k, v in clade2D.items():
            clade2seg   += (v[i:i+int(window)])
        for k, v in clade1D.items():
            clade1seg   += (v[i:i+int(window)])
        for k, v in cladeOutD.items():
            cladeOutseg += (v[i:i+int(window)])
       
        ## clade 1 As
        # test if clade1seg only contains As and clade2seg only contains As
        if (set(str(clade1seg)) <= A) and (set(str(clade2seg)) <= A):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                1
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)

        # test if clade1seg only contains As and clade2seg only contains Ts
        if (set(str(clade1seg)) <= A) and (set(str(clade2seg)) <= T):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)

        # test if clade1seg only contains As and clade2seg only contains Gs
        elif (set(str(clade1seg)) <= A) and (set(str(clade2seg)) <= G):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list) 

        # test if clade1seg only contains As and clade2seg only contains Cs
        elif (set(str(clade1seg)) <= A) and (set(str(clade2seg)) <= C):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("A")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)  

        ## clade 1 Ts
        # test if clade1seg only contains Ts and clade2seg only contains Ts
        elif (set(str(clade1seg)) <= T) and (set(str(clade2seg)) <= T):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                1
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)

        # test if clade1seg only contains Ts and clade2seg only contains As
        elif (set(str(clade1seg)) <= T) and (set(str(clade2seg)) <= A):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)

        # test if clade1seg only contains Ts and clade2seg only contains Gs
        elif (set(str(clade1seg)) <= T) and (set(str(clade2seg)) <= G):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list) 

        # test if clade1seg only contains Ts and clade2seg only contains Cs
        elif (set(str(clade1seg)) <= T) and (set(str(clade2seg)) <= C):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("T")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)

        ## clade 1 Gs
        # test if clade1seg only contains Gs and clade2seg only contains Gs
        elif (set(str(clade1seg)) <= G) and (set(str(clade2seg)) <= G):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                1
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)

        # test if clade1seg only contains Gs and clade2seg only contains As
        elif (set(str(clade1seg)) <= G) and (set(str(clade2seg)) <= A):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)

        # test if clade1seg only contains Gs and clade2seg only contains Ts
        elif (set(str(clade1seg)) <= G) and (set(str(clade2seg)) <= T):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list) 

        # test if clade1seg only contains Gs and clade2seg only contains Cs
        elif (set(str(clade1seg)) <= G) and (set(str(clade2seg)) <= C):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("G")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)  

        ## clade 1 Cs
        # test if clade1seg only contains Cs and clade2seg only contains Cs
        elif (set(str(clade1seg)) <= C) and (set(str(clade2seg)) <= C):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("C")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                1

        # test if clade1seg only contains Cs and clade2seg only contains As
        elif (set(str(clade1seg)) <= C) and (set(str(clade2seg)) <= A):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("A")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)

        # test if clade1seg only contains Cs and clade2seg only contains Ts
        elif (set(str(clade1seg)) <= C) and (set(str(clade2seg)) <= T):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("T")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list) 

        # test if clade1seg only contains Cs and clade2seg only contains Gs
        elif (set(str(clade1seg)) <= C) and (set(str(clade2seg)) <= G):
            ## print window that contains indel and concatenated string from both clades
            #print(i, i+int(window), "A->T", str(clade1seg), str(clade2seg), str(cladeOutseg))

            ## test if cladeOutseg contains A, T, C, G, N, gaps, or a combination
            # As in outgroup seg
            if any(x in cladeOutseg for x in set(A)) and not any(x in cladeOutseg for x in set(butA)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("A")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("Ts")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Ts in outgroup seg
            elif any(x in cladeOutseg for x in set(T)) and not any(x in cladeOutseg for x in set(butT)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("T")
                # class
                temp_list.append("C1C2")
                # C1type
                temp_list.append("Ts")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            # Gs in outgroup seg
            elif any(x in cladeOutseg for x in set(G)) and not any(x in cladeOutseg for x in set(butG)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("G")
                # class
                temp_list.append("C1")
                # C1type
                temp_list.append("Tv")
                # C2type
                temp_list.append("NA")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)
            elif any(x in cladeOutseg for x in set(C)) and not any(x in cladeOutseg for x in set(butC)):
                # append pertinent information to indel_arr
                temp_list = []
                # start
                temp_list.append(i)
                # stop
                temp_list.append(i+int(window))
                # clade 1
                temp_list.append("C")
                # clade 2
                temp_list.append("G")
                # outgroup
                temp_list.append("C")
                # class
                temp_list.append("C2")
                # C1type
                temp_list.append("NA")
                # C2type
                temp_list.append("Tv")
                # determine position of snp in codon
                if (i+int(window))%3==0:
                    temp_list.append("3rd")
                elif (i+int(window))%3==1:
                    temp_list.append("1st")
                elif (i+int(window))%3==2:
                    temp_list.append("2nd")
                snp_arr.append(temp_list)

    #print(snp_arr)
    # convert snp_arr to np array
    snp_arr = np.asarray(snp_arr)

    # if no snps are found, exit
    if len(snp_arr) == 0:
        sys.exit()

    # print np array line by line tab delimited
    for c0, c1, c2, c3, c4, c5, c6, c7, c8 in snp_arr:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(c0, c1, c2, c3, c4, c5, c6, c7, c8))

def read_clades_into_list(
    window, clade1, 
    clade2, alignedFasta,
    step, cladeOutgroup
    ):
    """
    Scans for small indels 
    
    Parameters
    ----------
    argv: window
        window size to go through fasta file
    argv: clade1
        file with taxa names from clade 1
        taxa names should be in a single column
    argv: clade2
        file with taxa names from clade 2
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
    clade1L   = []
    clade2L   = []
    cladeOutL = []

    # read clade files into lists
    clade1L   = [line.rstrip('\n') for line in open(clade1)]
    clade2L   = [line.rstrip('\n') for line in open(clade2)]
    cladeOutL = [line.rstrip('\n') for line in open(cladeOutgroup)]

    # determine number of indivs in each clade
    lenClade1 = len(clade1L)
    lenClade2 = len(clade2L)
    lenCladeOut = len(cladeOutL)

    # pass to find_indels function
    find_snps(
        window, clade1L, 
        clade2L, alignedFasta, 
        step, lenClade1, 
        lenClade2, cladeOutL,
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
    clade1        = ''
    clade2        = ''
    alignedFasta  = ''
    step          = 1
    cladeOutgroup = ''

    try:
        opts, args = getopt.getopt(argv, "hw:o:t:i:s:g:")
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
            # clade 2 file explanation
            print("\n-t\tclade 2 of interest:")
            print("\ta single column file with the fasta header names from clade 2")
            print("\tof interest. If the names do not match the fasta header, they")
            print("\twill not be recognized and accounted for.")
            # clade outgroup file explanation
            print("\n-g\tclade outgroup of interest:")
            print("\ta single columns file with the fasta header names from clade outgroup")
            print("\tof interest. This group should be the one used to contextualize the")
            print("\tSNP identified between clade 1 and clade 2. For example, if As are")
            print("\tobserved in clade 1 and Ts are observed in clade 2 and clade outgroup")
            print("\tthis script will report a T->A mutation in clade 1. In the same instance,")
            print("\tif the outgroup has Ns or -s, it will be reported as ambiguous.")
            # aligned fasta file explanation
            print("\n-i\taligned fasta file:")
            print("\taligned nucleotide fasta file containing (but not limited to)")
            print("\ttaxa specified in clade 1 (-o) and clade 2 (-t). If taxa are in")
            print("\tthe fasta file not specified in the clade 1 and clade 2 files,")
            print("\tthe taxa will not be considered in the analysis.")
            # output explanation
            print("\noutput explanation:")
            print("\tThe output will have 6 columns titled:")
            print("\tstart, stop, C1, C2, OG, class, C1type, C2type")
            print("\t  - start refers to the starting position of the SNP and stop is the end")
            print("\t  - C1 refers to the nucleotide found in clade 1 specified using -o")
            print("\t  - C2 refers to the nucleotide found in clade 2 specified using -t")
            print("\t  - OG refers to the nucleotide found in outgroup clade specified using -g")
            print("\t  - class will be populated with either a C1, C2, or C1C2. This refers to")
            print("\t\twhich clade has the SNP. C1 refers to clade1, C2")
            print("\t\trefers to clade2 and C1C2 refers to both clades.")
            print("\t  - C1type refers to if the SNP found at C1 is a Transition (Ts) or")
            print("\t\tor Transversion (Tv). If no SNP is found, the column will contain NA")
            print("\t  - C2type contains the same information but for clade 2.")
            print("\t  - pos refers to the position of the SNP in a codon. That is, is the SNP in")
            print("\t\tthe wobble position (3rd) or the 1st or 2nd position.\n")
            sys.exit()
        elif opt == '-o':
            if os.path.isfile(arg):
                clade1 = arg
                #print("\nFile with taxa from clade 1: {}".format(clade1))
            else:
                # error message
                print("\n\nThe specified file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                clade2 = arg
                #print("\nFile with taxa from clade 2: {}".format(clade2))
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
        window, clade1, 
        clade2, alignedFasta,
        step, cladeOutgroup
        )

if __name__ == '__main__':
    main(sys.argv[1:])