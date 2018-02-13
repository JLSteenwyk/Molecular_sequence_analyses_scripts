#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio import SeqIO
import numpy as np


def print_results(
    ID, IndelSNPSeq, OGSeq
    ):
    """
    summarizes differences between the two 
    nucleotide sequences at the protein level
    
    Parameters
    ----------
    argv: ID
        ID of sequence
    argv: IndelSNPSeq
        new sequence
    argv: OGSeq
        original sequence
    """
   
    # make sure sequences are divisible by three
    if len(OGSeq) % 3 != 0:
        # determine the number of Ns (remainder) to add to create a codon length gene
        remainder=3-(len(OGSeq)%3)
        # add Ns to end of Seq
        OGSeq=OGSeq+(remainder*'N')
    if len(IndelSNPSeq) % 3 != 0:
        # determine the number of Ns (remainder) to add to create a codon length gene
        remainder=3-(len(IndelSNPSeq)%3)
        # add Ns to end of Seq
        IndelSNPSeq=IndelSNPSeq+(remainder*'N')

    # make OGSeq and IndelSNPSeq the same length
    OGlen          = len(OGSeq)
    IndelSNPSeqLen = len(IndelSNPSeq)
    if OGlen > IndelSNPSeqLen:
        difference=OGlen-IndelSNPSeqLen
        IndelSNPSeq=IndelSNPSeq+(difference*'N')
    elif OGlen < IndelSNPSeqLen:
        difference=IndelSNPSeqLen-OGlen
        OGSeq=OGSeq+(difference*'N')

    # reassign seqs with the amino acid sequence
    OGSeq        = OGSeq.translate()
    IndelSNPSeq  = IndelSNPSeq.translate()
    
    # result lists
    res_arr   = []
    res_arr.append(["Pos", "Old", "New"])
    # summarize differences between the two AA strings
    for AA in range(0, (int(len(OGSeq))), 1):
        SeqAA         = OGSeq[AA:AA+1]
        IndelSNPSeqAA = IndelSNPSeq[AA:AA+1]
        #if strings are not equal
        if SeqAA != IndelSNPSeqAA:
            temp_arr = []
            temp_arr.append(AA+1)
            temp_arr.append(SeqAA[0])
            temp_arr.append(IndelSNPSeqAA[0])
            res_arr.append(temp_arr)

    # check if res_arr has no new additions and exit if so
    if res_arr==[["Pos", "Old", "New"]]:
        sys.exit()

    # convert to numpy array
    res_arr = np.asarray(res_arr)

    # print old and new sequence
    print(">"+ID+"|old"+"\n"+OGSeq)
    print(">"+ID+"|new"+"\n"+IndelSNPSeq+"\n")

    # print np array line by line tab delimited and show differences
    # between the old and new sequence
    for c0, c1, c2 in res_arr:
        print("{}\t{}\t{}".format(c0, c1, c2))

    # determine percentage of sequence that differences from original sequence
    percentage=float((len(OGSeq)/(len(res_arr)-1))*100)
    print("\n"+str(len(res_arr)-1)+"/"+str(len(OGSeq)), \
        "\nnucleotide differences between old and new sequence")


def indel_of_nucl(
    SNPSeq, ID,
    INDEL_arr, OGSeq
    ):
    """
    modifies fasta sequence to have the indel
    
    Parameters
    ----------
    argv: SNPSeq
        sequence that has the snps included in the fasta file
    argv: ID
        fasta entry ID
    argv: INDEL_arr
        array of indel info 
    argv: Seq
        original nucleotide sequence
    """

    # intialize holders for coordinate (COORD), type (TYPE), modified nucleotide (MOD), and IndelSNPSeq
    COORD       = ''
    TYPE        = ''
    MOD         = ''
    IndelSNPSeq = ''

    # initialize counter for how many total nucleotides have been removed or deleted
    indelCNT = 0

    # if there are indels, modify the sequence accordingly
    if len(INDEL_arr) > 0:
        # loop through indel mutations and modify the sequence accordingly
        for indel in INDEL_arr:

            # assign coordinate (COORD), old nucleotide (OLD), and new/snp nucleotide (NEW)
            COORD = indel[1]
            TYPE  = indel[2]
            MOD   = indel[3]
            COORD = int(COORD)-1+indelCNT

            # if a deletion check that the nucleotide matches
            if TYPE == '-':
                # check if coordinate matches the old nucleotide
                if SNPSeq[COORD] == SNPSeq[COORD+len(MOD)]: 
                    # populate indelSNPSeq with SNPSeq without the MOD nucleotide 
                    IndelSNPSeq = SNPSeq[:COORD] + SNPSeq[COORD + 1:]
                    # reset SNPSeq to be the new nucleotide sequence
                    SNPSeq = IndelSNPSeq
                    # subtract length of MOD sequence
                    indelCNT=indelCNT-len(MOD)
                else:
                    print(SNPSeq[(COORD):(COORD+2)])
                    print("\nPosition", COORD+1, "is a", SNPSeq[COORD], "and not a", MOD)
                    print("Check that the correct position has been specified\n")
                    sys.exit()

            # if a deletion check that the nucleotide matches
            elif TYPE == '+':
                # populate indelSNPSeq with SNPSeq without the MOD nucleotide 
                IndelSNPSeq = SNPSeq[:COORD+1] + MOD + SNPSeq[COORD + 1:]
                # reset SNPSeq to be the new nucleotide sequence
                SNPSeq = IndelSNPSeq
                # add length of MOD sequence
                indelCNT=indelCNT+len(MOD)

        # pass to print_results function
        print_results(
            ID, IndelSNPSeq, OGSeq
            )

    # if there are no indels, pass to print_results function
    elif len(INDEL_arr) == 0:
        IndelSNPSeq = SNPSeq
        # pass to print_results
        print_results(
            ID, IndelSNPSeq, OGSeq
            )

def replace_nucl_with_SNP(
    FASTA, ENTRY,
    SNP_arr, INDEL_arr
    ):
    """
    replaces nucleotide with SNP
    
    Parameters
    ----------
    argv: FASTA
        fasta file
    argv: ENTRY
        fasta entry ID
    argv: SNP_arr
        array of snp info 
    argv: OLD
        array of indel info
    """

    # read in fasta file
    FAdict = {}
    format = "fasta"
    handle = open(FASTA)
    FAdict = SeqIO.to_dict(SeqIO.parse(handle, format))

    # intialize ID and seq holder
    ID     = ''
    Seq    = ''
    SNPSeq = ''

    # intialize holders for coordinate (COORD), old nucleotide (OLD), and new/snp nucleotide (NEW)
    COORD  = ''
    OLD    = ''
    NEW    = ''

    # populate ID and Seq
    ID     = FAdict[ENTRY].id
    Seq    = FAdict[ENTRY].seq
    OGSeq  = FAdict[ENTRY].seq
    
    # if length of SNP_arr is greater than 0 then replace snps
    if len(SNP_arr) > 0:
        # loop through SNP mutations and replace coordinate with sequence
        for snp in SNP_arr:

            # assign coordinate (COORD), old nucleotide (OLD), and new/snp nucleotide (NEW)
            COORD = snp[1]
            OLD   = snp[2]
            NEW   = snp[3]

            # check if coordinate matches the old nucleotide
            COORD=int(COORD)-1
            if Seq[COORD] == OLD:
               # create string with new sequence
                SNPSeq=Seq[:COORD] + NEW + Seq[COORD + 1:]
                # reset Seq to have the SNPs
                Seq=SNPSeq
            else:
                print("Position", COORD+1, "is a", Seq[COORD], "and not a", OLD)
                print("Check that the correct position has been specified")
                sys.exit()

        # pass to indel_of_nucl
        indel_of_nucl(
            SNPSeq, ID,
            INDEL_arr, OGSeq
            )

    # if SNP_arr is of length 0, then move to next function
    elif len(SNP_arr) == 0:
        SNPSeq = OGSeq
        # pass to indel_of_nucl
        indel_of_nucl(
            SNPSeq, ID,
            INDEL_arr, OGSeq
            )


    

def read_mutat_file(
    FASTA, ENTRY,
    MUTAT
    ):
    """
    read mutation file and parse out indels and snps
    
    Parameters
    ----------
    argv: FASTA
        fasta file
    argv: ENTRY
        fasta entry ID
    argv: MUTAT
        mutation summary file
    """

    # initialize lists for each mutation
    SNP_arr   = []
    INDEL_arr = []
    # initialize line counter
    lineCNT   = 1

    with open(MUTAT) as mutSummary:
        for line in mutSummary:
            line = re.split(r'\t+', line.rstrip('\n'))
            if line[0] == "snp":
                SNP_arr.append(line)
            elif line[0] == "indel":
                INDEL_arr.append(line)
            elif line[0] not in ("snp", "indel"):
                print("Line", lineCNT, "does not contain 'snp' or 'indel'")
                print("Please check the mutation summary (-m) file")
                print("Exiting now...")
                sys.exit()
            lineCNT+=1

    # pass to replace_nucl_with_SNP function
    replace_nucl_with_SNP(
        FASTA, ENTRY,
        SNP_arr, INDEL_arr
        )

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    FASTA  = ''
    ENTRY  = ''
    MUTAT  = ''

    try:
        opts, args = getopt.getopt(argv, "hf:e:m:")
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
            # mutation summary file
            print("\n-m\tmutation summary file:")
            print("\tThe mutation summary file should be tab delimited and in the following format:")
            print("\tcol 1: snp/indel")
            print("\tcol 2: coordinate")
            print("\tcol 3: for snps, the original nucleotide. for indels, +/-")
            print("\tcol 4: for snps, new nucleotide. for indels, inserted or deleted nucleotide.")
            print("\tFor example, two lines of a mutation summary file should be the following:")
            print("\tsnp\t2\tT\tG")
            print("\tindel\t3\t-\tG\n")
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
        # read in mutation summary file name
        elif opt == '-m':
            if os.path.isfile(arg):
                MUTAT = arg
            else:
                # error message
                print("\n\nThe specified mutation summary file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()


    # pass to read_mutat_file function
    read_mutat_file(
        FASTA, ENTRY,
        MUTAT
        )

if __name__ == '__main__':
    main(sys.argv[1:])