#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio import SeqIO
import numpy as np

def identify_mutations(
        collapsed_indels
        ):
    """
    identifies and characterizes mutations in homopolymers

    Parameters
    ----------
    argv: collapsed_indels
        collapsed homopolymer runs with taxon and outgroup sequences
    """

    # initialize summary array
    summary_arr = [['start', 'stop', 'len', 'hom', 'homlen', 'mut', 'subs', 'ins', 'del', 'subsTr', 'insTr', 'delTr', 'taxon', 'outgroup']]

    ## identify mutational landscape of taxon Seq compared to outgroup Seq
    # loop through each entry in collapsed_indels
    for line in collapsed_indels:
        # extract sequences
        taxonSeq    = line[3].upper()
        outgroupSeq = line[4].upper()
        # initialize substition, indel, and total muts counters
        subs = 0
        inse = 0
        dele = 0
        muts = 0

        # initialize lists to keep track of where
        # mutations occur in homopolymers
        subs_tr = []
        inse_tr = []
        dele_tr = []

        # loop through each character in the outgroup and taxon seqs
        for tChar, oChar in zip(taxonSeq, outgroupSeq):
            # if tChar and oChar are the same then add 0s to position 
            # tracker lists for mutations
            if tChar == oChar:
                subs_tr.append(0)
                inse_tr.append(0)
                dele_tr.append(0)
            # if tChar and oChar aren't the same, check if subs or indel
            else:
                # check for deletion in tChar
                if (tChar == '-') and (oChar != '-'):
                    dele += 1
                    subs_tr.append(0)
                    inse_tr.append(0)
                    dele_tr.append(1)
                # check for insertion in tChar
                elif (tChar != '-') and (oChar == '-'):
                    inse += 1
                    subs_tr.append(0)
                    inse_tr.append(1)
                    dele_tr.append(0)
                # check for subs in tChar
                elif (tChar != '-') and (oChar != '-'):
                    subs += 1
                    subs_tr.append(1)
                    inse_tr.append(0)
                    dele_tr.append(0)
        
        # initialize line for the tracker 
        subs_tr = ','.join(map(str, subs_tr)) 
        inse_tr = ','.join(map(str, inse_tr))
        dele_tr = ','.join(map(str, dele_tr))
            
        # # correct the number of indels counted by codon size
        # dele = dele/3
        # inse = inse/3
        
        # sum the number of mutations by adding indels and subs
        muts = subs + dele + inse
        
        # determine length of the homopolymer
        leng = 0
        leng = len(outgroupSeq)

        # determine length of homopolymer excluding any insertions in the homopolymer
        homlen = 0
        homlen = len(outgroupSeq.replace('-', ''))

        # create line for summary array and add it to summary array
        # the line, temp_list, is created to have the following order
        # start, stop, leng, hom, homlen, mut, subs, ins, dele, subsTr, inseTr, deleTr, taxonSeq, outgroupSeq
        temp_list = []
        temp_list.append(line[0])
        temp_list.append(line[1])
        temp_list.append(leng)
        temp_list.append(line[2])
        temp_list.append(homlen)
        temp_list.append(int(muts))
        temp_list.append(int(subs))
        temp_list.append(int(inse))
        temp_list.append(int(dele))
        temp_list.append(subs_tr)
        temp_list.append(inse_tr)
        temp_list.append(dele_tr)
        temp_list.append(taxonSeq)
        temp_list.append(outgroupSeq)
        summary_arr.append(temp_list)

        #print(temp_list, dele)

    ## remove homopolymers of length 1
    # initialize list of idx with length 1 and counter
    idx_len1 = []
    cnt      = 0
    for line in summary_arr:
        if line[4] == 1:
            idx_len1.append(cnt)
            cnt += 1    
        else:
            cnt += 1

    # Convert to np array
    summary_arr = np.asarray(summary_arr)
    # remove indices
    for index in sorted(idx_len1, reverse=True):
        del summary_arr[index]
        
    # print np array line by line tab delimited
    for start, stop, leng, hom, homlen, mut, subs, ins, dele, subsTr, inseTr, deleTr, taxonSeq, outgroupSeq in summary_arr:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(start, stop, leng, hom, homlen, mut, \
            subs, ins, dele, subsTr, inseTr, deleTr, taxonSeq, outgroupSeq))

def identify_homopolymer_seqs(
        collapsed_indels, alignment,
        taxonL, outgroupL
        ):
    """
    identifes the sequences in the outgroup and ingroup taxa
    for each homopolymer run

    Parameters
    ----------
    argv: collapsed_indels
        collapsed homopolymer runs accounting for indels
    argv: alignment
        alignment file
    """

    # obtain taxon Seq and representative outgroup Seq
    taxonSeq    = ''
    outgroupSeq = ''
    for indiv in alignment:
        if indiv.id in taxonL:
            taxonSeq = indiv.seq._data
        elif indiv.id in outgroupL[0]:
            outgroupSeq = indiv.seq._data
    
    for line in collapsed_indels:
        seqStart = line[0]-1
        seqEnd   = line[1]-1
        line.append(taxonSeq[seqStart:seqEnd])
        line.append(outgroupSeq[seqStart:seqEnd])

    # identify mutations in homopolymers
    identify_mutations(
        collapsed_indels
        )

def collapse_indels(
        collapsed_data, alignment, 
        taxonL, outgroupL
        ):
    """
    collapses homopolymers that are adjacent to gaps

    Parameters
    ----------
    argv: collapsed_data
        collapsed array of homopolymers
    argv: alignment
        alignment file
    """

    ## Looping over lines to collapsing homopolymers and
    ##  those with gaps between them and those with gaps before them
    collapsed_indels = []
    idx_ii           = 0
    idx_jj           = idx_ii + 1
    idx_opt          = True
    while_opt        = True
    # Looping over first line
    while(idx_opt):
        try:
            line_ii = collapsed_data[idx_ii]
            line_jj = collapsed_data[idx_jj]
            # Testing condition
            line_ii_opt = True
            # loop through collapsed_data
            while(line_ii_opt):
                new_line = line_ii
                ## collapse gaps that occur within or after homopolymers
                # if end line 1 (ii) is greater or equal to line 2 (jj) and
                # line 1 conserved sequence is equal to line 2 conserved sequence
                # then collapse the windows
                if (line_ii[1] >= line_jj[0]) and ((line_ii[2]==line_jj[2]) or (line_jj[2]=='-')):
                    idx_jj  += int(1)
                    # new line is line 1 start, line 2 stop, line 2 outgroup ID, line 2 gap ID
                    new_line = [line_ii[0], line_jj[1], line_ii[2]]
                    line_ii  = new_line
                    line_jj  = collapsed_data[idx_jj]
                ## collapse gaps that occur before homopolymers
                elif (line_ii[1] >= line_jj[0]) and (line_ii[2]=='-'):
                    idx_jj  += int(1)
                    # new line is line 1 start, line 2 stop, line 2 outgroup ID, line 2 gap ID
                    new_line = [line_ii[0], line_jj[1], line_jj[2]]
                    line_ii  = new_line
                    line_jj  = collapsed_data[idx_jj]                
                else:
                    line_ii_opt = False
                    collapsed_indels.append(new_line)
                    idx_ii = idx_jj
                    idx_jj = idx_ii + 1
        # except for end of indel_arr
        except IndexError:
            collapsed_indels.append(line_ii)
            idx_opt = False

    ## remove homopolymers of gaps
    # initialize list of idx with length 1 and counter
    idxes = []
    cnt   = 0
    for line in collapsed_indels:
        if set(line[2]) == set('-'):
            idxes.append(cnt)
            cnt += 1    
        else:
            cnt += 1
    # remove indices
    for index in sorted(idxes, reverse=True):
        del collapsed_indels[index]

    ## remove polymers of length 1
    # initialize list of idx with length 1 and counter
    idx_len1 = []
    cnt      = 0
    for line in collapsed_indels:
        if line[1]-line[0] == 1:
            idx_len1.append(cnt)
            cnt += 1    
        else:
            cnt += 1
    # remove indices
    for index in sorted(idx_len1, reverse=True):
        del collapsed_indels[index]

    identify_homopolymer_seqs(
        collapsed_indels, alignment,
        taxonL, outgroupL
        )

def collapse_homopolymers(
        conserved_arr, alignment, 
        taxonL, outgroupL
        ):
    """
    collapsed adjacent windows with the same
    character in conserved_arr

    Parameters
    ----------
    argv: conserved_arr
        array of conserved sites in the outgroup
    argv: alignment
        alignment file
    """

    ## Looping over lines to collapsing homopolymers
    collapsed_data = []
    idx_ii      = 0
    idx_jj      = idx_ii + 1
    idx_opt     = True
    while_opt   = True
    # Looping over first line
    while(idx_opt):
        try:
            line_ii = conserved_arr[idx_ii]
            line_jj = conserved_arr[idx_jj]
            # Testing condition
            line_ii_opt = True
            # loop through conserved_arr
            while(line_ii_opt):
                new_line = line_ii
                # if end line 1 (ii) is greater or equal to line 2 (jj) and
                # line 1 conserved sequence is equal to line 2 conserved sequence
                # then collapse the windows
                if (line_ii[1] >= line_jj[0]) and (line_ii[2]==line_jj[2]):
                    idx_jj  += int(1)
                    # new line is line 1 start, line 2 stop, line 2 outgroup ID, line 2 gap ID
                    new_line = [line_ii[0], line_jj[1], line_jj[2]]
                    line_ii  = new_line
                    line_jj  = conserved_arr[idx_jj]
                else:
                    line_ii_opt = False
                    collapsed_data.append(new_line)
                    idx_ii = idx_jj
                    idx_jj = idx_ii + 1
        # except for end of conserved_arr
        except IndexError:
            collapsed_data.append(line_ii)
            idx_opt = False

    # pass to collapse indels function
    collapse_indels(
        collapsed_data, alignment,
        taxonL, outgroupL
        )

def find_conserved_positions(
    taxonL, 
    fasta, 
    outgroupL,
    ):
    """
    Scans for positions with conserved characters
    
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
    outgroupD  = {}

    # initialize length variable
    length = ''
    
    # read in aligned fasta file
    format    = "fasta"
    handle    = open(fasta)
    alignment = list(SeqIO.parse(handle, format))

    # fill taxaD dictionary with individual id (key) and sequence (value)
    for indiv in alignment:
        if indiv.id in outgroupL:
            outgroupD[indiv.id] = indiv.seq

    # determine sequence length
    length    = len(alignment[0].seq)

    # intialize list for np array of continuous chars identified
    conserved_arr = []

    # loop through aligned sequence using step size
    for i in range(0, (int(length)+1) - int(step), int(step)):
        # initialize variables to hold sequence of interest being analyzed
        outgroupSeg = ''

        # loop through all individuals from taxa of interest and save 
        # window being analyzed into a concatenated string from all indivs
        for k, v in outgroupD.items():
            outgroupSeg += (v[i:i+int(window)])
            # make all characters the same case
            outgroupSeg=outgroupSeg.upper()
        outgroupSeg=set(outgroupSeg)
       
        # test if outgroupSeg only contains 1 character or more
        outgroupChar=','.join(str(s) for s in outgroupSeg)
        if (len(outgroupSeg) < 2):
            # append pertinent information to indel_arr
            temp_list = []
            # start
            temp_list.append(i+1)
            # stop
            temp_list.append(i+int(window)+1)
            # char
            temp_list.append(outgroupChar)
            conserved_arr.append(temp_list)

    # if no snps are found, exit
    if len(conserved_arr) == 0:
        sys.exit()

    # pass to collapse homopolymers function
    collapse_homopolymers(
        conserved_arr, alignment,
        taxonL, outgroupL
        )

def read_clades_into_list(
    taxon, 
    fasta,
    outgroup
    ):
    """
    reads taxa into lists
    
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

    # pass to find_conserved_positions function
    find_conserved_positions(
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
    taxon        = ''
    fasta        = ''
    outgroup     = ''

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
            #script explanation
            print("\nThis script will identify homopolymer runs in the outgroup taxa sequence (-g parameter).")
            print("The script will then report substitutions if there are any in a taxon of interest (-t parameter).")
            print("Additionally the script will report any insertions or deletions in homopolymer sequences.")
            print("")
            # clade 1 file explanation
            print("\n-t\ttaxon of interest:")
            print("\ta single column file with the fasta header names from clade 1")
            print("\tof interest. If the names do not match the fasta header, they")
            print("\twill not be recognized and accounted for.")
            # clade outgroup file explanation
            print("\n-g\tclade outgroup of interest:")
            print("\ta single columns file with the fasta header names from clade outgroup")
            print("\tof interest. This group should be the one used to contextualize the")
            print("\tSNP identified between clade 1. For example, if As are")
            print("\tobserved in clade 1 and Ts and clade outgroup this script will report a")
            print("\tT->A mutation in clade 1. In the same instance,")
            print("\tif the outgroup has Ns or -s, it will be reported as ambiguous.")
            # aligned fasta file explanation
            print("\n-i\taligned fasta file:")
            print("\taligned fasta file containing (but not limited to)")
            print("\ttaxon of interest (-t) and the outgroup (-g). If taxa are in")
            print("\tthe fasta file not specified in the clade 1 and outgroup files,")
            print("\tthe taxa will not be considered in the analysis.")
            # output explanation
            print("\noutput explanation:")
            print("\tThe output will have 7 columns titled:")
            print("\tstart, stop, hom, len, mut, subs, ins, del, subsTr, insTr, delTr, Cl, OG")
            print("\t  - start refers to the starting position of the homopolymer in the outgroup")
            print("\t  - stop refers to the ending position of the homopolymer outgroup")
            print("\t  - len refers to the length of the homopolymer including gaps")
            print("\t  - hom refers to what the homopolymer character is")
            print("\t  - homlen refers to how long the homopolymer is excluding gaps")
            print("\t  - mut refers to whether or not a mutation was identified in the taxon of interest")
            print("\t  - subs refers to whether or not there is a substition found in the taxon")
            print("\t\tof interest. If there is no substition, the script will report '0'")
            print("\t  - ins refers to whether or not there is a insertion found in the taxon")
            print("\t\tof interest. If there is no insertion, the script will report '0'. If not 0, it is the")
            print("\t\tnumber of inserted nucleotides.")
            print("\t  - del refers to whether or not there is a deletion found in the taxon")
            print("\t\tof interest. If there is no deletion, the script will report '0'. If not 0, it is the")
            print("\t\tnumber of deleted nucleotides")
            print("\t  - subsTr refers to a tracker for where substituions are found. This is formatted")
            print("\t\tin a boolean way such that a homopolymer of length three with a substituion in the")
            print("\t\tfirst position will be reported as 1,0,0")
            print("\t  - insTr refers to a tracker for where insertions are found. This is formatted")
            print("\t\tin the same way as subsTr")
            print("\t  - delTr refers to a tracker for where insertions are found. This is formatted")
            print("\t\tin the same way as subsTr")
            print("\t  - Cl refers to the sequence found in ingroup taxon specified using -t")
            print("\t  - OG refers to the sequence found in outgroup clade specified using -g")
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

    # pass to read_clades_into_list function
    read_clades_into_list(
        taxon, 
        fasta,
        outgroup
        )

if __name__ == '__main__':
    main(sys.argv[1:])