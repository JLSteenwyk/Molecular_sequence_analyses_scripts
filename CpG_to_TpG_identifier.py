#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import SeqIO
import numpy as np


def print_res(
    summary_arr
    ):
    """
    Prints out results of GC dinucleotides with and without mutations

    Parameters
    ----------
    argv: summary_arr
        array with results of GC dinucleotides with and without mutations
        summary_arr is organized in the following way:
        1: start (counting from 0)
        2: end 
        3: outgroup sequence
        4: taxon of interest sequence
        5: boolean for mutation (1=mutation ; 0=no mutation)
        6: the taxon of interest nucleotide at the C position
    """

    # for start and stop add 1 to account for python counting from 0
    for line in summary_arr:
    	line[0] += 1
    	line[1] += 1

    # add a header to summary_arr
    summary_arr.insert(0, ['start', 'end', 'og', 'tax', 'mut', 'CpGmut'])

    ## remove lists that aren't 6 in length
    # initialize list of idx with length 1 and counter
    idxes = []
    cnt   = 0
    for line in summary_arr:
        if len(line) != 6:
            idxes.append(cnt)
            cnt += 1    
        else:
            cnt += 1
    # remove indices
    for index in sorted(idxes, reverse=True):
        del summary_arr[index]

    # Convert to np array
    summary_arr = np.asarray(summary_arr)

    # print np array line by line tab delimited
    for start, end, og, tax, mut, CpGmut in summary_arr:
        print("{}\t{}\t{}\t{}\t{}\t{}".format(start, end, mut, CpGmut, tax, og))


def filter_overlapping(
    GCdinucleotides_arr
    ):
    """
    Loops through GCdinucleotides_arr and filters
    out any substitutions that are being double counted
    because they are adjacent GC sequences

    Parameters
    ----------
    argv: GCdinucleotides_arr
        array of GC dinucleotides and the taxon of
        interest dinucleotide sequence
    """

    ## Looping over lines to filter out double counted substitutions
    GCdinucl_filtered = []
    idx_ii            = 0
    idx_jj            = idx_ii + 1
    idx_opt           = True
    while_opt         = True
    # Looping over first line
    while(idx_opt):
        try:
            line_ii = GCdinucleotides_arr[idx_ii]
            line_jj = GCdinucleotides_arr[idx_jj]
            # Testing condition
            line_ii_opt = True
            # loop through collapsed_data
            while(line_ii_opt):
                new_line = line_ii
                if (line_ii[1] > line_jj[0]) and ((line_ii[4]==1) and (line_jj[4]==1)):
                    idx_jj  += int(1)
                    new_line = line_ii
                    line_jj  = GCdinucleotides_arr[idx_jj]
                else:
                    line_ii_opt = False
                    GCdinucl_filtered.append(line_ii)
                    idx_ii = idx_jj
                    idx_jj = idx_ii + 1
        # except for end of indel_arr
        except IndexError:
            GCdinucl_filtered.append(new_line)
            idx_opt = False


    # initialize summary array
    summary_arr = [] 
    # remove redundant elements in GCdinucl_filtered
    for ele in GCdinucl_filtered: 
        if ele not in summary_arr: 
            summary_arr.append(ele) 
    
    print_res(
    	summary_arr
    	)

def find_muts(
    GCdinucleotides_arr
    ):
    """
    Loops through GCdinucleotides_arr and finds
    those with any type of mutation

    Parameters
    ----------
    argv: GCdinucleotides_arr
        array of GC dinucleotides and the taxon of
        interest dinucleotide sequence
    """

    # loop through GC dinucleotides
    for line in GCdinucleotides_arr:
        # if the outgroup sequence is GC
        if line[2] == 'GC':
            # if the taxon sequence is AC
            if line[3][0] == 'A':
                # boolean for a mutation
                line.append(1)
                # what the CpG mutation is
                line.append('T')
            # if the taxon sequence is TC
            elif line[3][0] == 'T':
                # boolean for a mutation
                line.append(1)
                # what the CpG mutation is
                line.append('A')
            # if the taxon sequence is GG
            elif line[3][0] == 'C':
                # boolean for a mutation
                line.append(1)
                # what the mutation is
                line.append('G')
            # if there is no mutation
            elif line[3][0] == 'G':
                # boolean for a mutation
                line.append(0)
                # what the mutation is
                line.append('C')
            # if there is a gap
            elif line[3][0] == '-':
                # boolean for a mutation
                line.append(1)
                # what the mutation is
                line.append('-')
        # if the outgroup sequence is CG
        if line[2] == 'CG':
            # if the taxon sequence is GT
            if line[3][0] == 'T':
                # boolean for a mutation
                line.append(1)
                # what the CpG mutation is
                line.append('T')
            # if the taxon sequence is GA
            elif line[3][0] == 'A':
                # boolean for a mutation
                line.append(1)
                # what the CpG mutation is
                line.append('A')
            # if the taxon sequence is GG
            elif line[3][0] == 'G':
                # boolean for a mutation
                line.append(1)
                # what the CpG mutation is
                line.append('G')
            # if there is no mutation
            elif line[2][0] == 'C':
                # boolean for a mutation
                line.append(0)
                # what the CpG mutation is
                line.append('C')
            # if there is a gap
            elif line[3][0] == '-':
                # boolean for a mutation
                line.append(1)
                # what the mutation is
                line.append('-')

    # new order of information for GCdinucleotides_arr is the following
    # 1: start (counting from 0)
    # 2: end
    # 3: outgroup sequence
    # 4: taxon of interest sequence
    # 5: boolean for mutation (1=mutation ; 0=no mutation)
    # 6: the taxon of interest nucleotide at the C position

    filter_overlapping(
    	GCdinucleotides_arr
    	)

def read_GCdinucleotides(
    taxonL, 
    fasta, 
    outgroupL,
    ):
    """
    Reads the fasta file into a list of lists of dinucleotides
    for the outgroup sequences
    
    Parameters
    ----------
    argv: taxonL
        list with taxon of interest
    argv: fasta
        input fasta file
    argv: outgroupL
        list of outgroup taxa
    """

    ## read in alignment
    # initialize taxa dictionaries
    outgroupD = {}
    taxonSeq  = ''

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
    # fill taxaD dictionary with individual id (key) and sequence (value)
    for indiv in alignment:
        if indiv.id in taxonL:
            taxonSeq = (indiv.seq).upper()

    ## loop through alignment
    # determine sequence length
    length    = len(alignment[0].seq)

    # intialize list for np array of continuous chars identified
    GCdinucleotides_arr = []

    # initialize step and window
    step   = 1
    window = 2

    # loop through aligned sequence using step size
    for i in range(0, (int(length)+1) - int(step), int(step)):
        # initialize variable to hold outgroup dinucleotides of interest being analyzed
        outgroupDict   = {}
        outgroupDiNucl = ''

        # loop through all individuals from taxa of interest and save 
        # window being analyzed into a concatenated string from all indivs
        for k, v in outgroupD.items():
            # skip the last site because it won't be a dinucleotide
            if len(v[i:i+int(window)]) < 2:
                continue
            else:
            	# create dictionary of each line
                outgroupDict[k] = (v[i:i+int(window)]._data).upper()
                outgroupDiNucl  = (v[i:i+int(window)]._data).upper()
        # if all outgroups have the same dinucleotide sequences then
        # keep the dinucleotide  
        if len(dict(zip(outgroupDict.values(),outgroupDict.keys())))==1 and \
            (outgroupDiNucl=='GC' or outgroupDiNucl=='CG'):
            # append pertinent information to GCdinucleotides_arr
            temp_list = []
            # start
            temp_list.append(i)
            # stop
            temp_list.append(i+int(window))
            # outgroup dinucleotide
            temp_list.append(outgroupDiNucl)
            # taxon dinucleotide
            temp_list.append(taxonSeq[i:i+int(window)]._data)
            # add to GCdinucleotides_arr
            GCdinucleotides_arr.append(temp_list)

    find_muts(
    	GCdinucleotides_arr
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

    # pass to read_GCdinucleotides function
    read_GCdinucleotides(
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
            print("\nThis script will parse through a fasta file (-i parameter) and identify CpG -> TpG")
            print("mutations suggestive of 5-mC mutations. To do so, it will first identify CpGs among")
            print("outgroup taxa (-g parameter) and then compare the outgroup sequence to a taxon of ")
            print("interest (-t parameter).")
            print("\nOutput is in the following order: start, end, mut, CpGmut, tax, og")
            print("start: CpG start")
            print("end: CpG end")
            print("mut: boolean for if there is a mutation")
            print("CpGmut: in the context of CpG, what is the C mutated to")
            print("tax: taxon sequence")   
            print("og: outgroup sequence")
            # taxon of interest explanation
            print("\n-t\ttaxon of interest:")
            print("\ta string argument that matches the fasta entry identifier for the taxon of interest")
            # clade outgroup file explanation
            print("\n-g\tclade outgroup of interest:")
            print("\ta single columns file with the fasta header names from clade outgroup")
            print("\tof interest. This group should be the one used to contextualize the")
            print("\tSNP identified in the taxon of interest.")
            # aligned fasta file explanation
            print("\n-i\taligned fasta file:")
            print("\taligned fasta file containing taxon of interest (-t) and the outgroup (-g).")
            # output explanation
            print("\noutput explanation:")
            print("\tThe output will have x columns titled:")
            print("\txx")
            print("\t")
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