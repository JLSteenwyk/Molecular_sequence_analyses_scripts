#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio import SeqIO
from Bio.Seq import Seq

def translate(
    inFasta, TTable
    ):
    """
    Translates nucleotide sequences into protein sequences.

    Parameters
    ----------
    argv: inFasta
        input fasta file
    argv: TTable
    	translation table NCBI code
    """
    # read in fasta file
    format   = "fasta"
    handle   = open(inFasta)
    fastaSeq = list(SeqIO.parse(handle, format))

    # loop through fasta entries and translate each 
    for indiv in fastaSeq:
        header=''
        header=">"+indiv.description
        print(header)
        pSeq = str(indiv.seq.translate(table=TTable))
        pSeq = pSeq.replace("*", "")
        print(pSeq)

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    window        = ''
    clade1        = ''
    clade2        = ''
    alignedFasta  = ''
    step          = ''
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
            # general description
            print("\nThis script will take as an input file a multi-fasta file and translate the nucleotides")
            print("into their corresponding protein sequences. The input file is specified with the -i argument")
            print("and the -t argument specifies the translation table.\n")
            # input fasta file explanation
            print("\n-i\tfasta file:")
            print("\tinput fasta file of nucleotide sequences")
            # explanation of boolean for whether or not to include or exclude stop codons
            print("\n-s\tboolean for stop codon")
            print("\tT or F for keeping or removing a stop codon in a sequence, respectively.")
            # translation table explanattion
            print("\n-t\ttranslation table:")
            # translation code explanations
            print("\tTranslation tables are available through NCBI. Argument should be a number and")
            print("\tnot a string that describes the table.\n")
            print("\tThe following genetic codes are described here:")
            print("\t\t1. The Standard Code")
            print("\t\t2. The Vertebrate Mitochondrial Code")
            print("\t\t3. The Yeast Mitochondrial Code")
            print("\t\t4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code")
            print("\t\t5. The Invertebrate Mitochondrial Code")
            print("\t\t6. The Ciliate, Dasycladacean and Hexamita Nuclear Code")
            print("\t\t9. The Echinoderm and Flatworm Mitochondrial Code")
            print("\t\t10. The Euplotid Nuclear Code")
            print("\t\t11. The Bacterial, Archaeal and Plant Plastid Code")
            print("\t\t12. The Alternative Yeast Nuclear Code")
            print("\t\t13. The Ascidian Mitochondrial Code")
            print("\t\t14. The Alternative Flatworm Mitochondrial Code")
            print("\t\t16. Chlorophycean Mitochondrial Code")
            print("\t\t21. Trematode Mitochondrial Code")
            print("\t\t22. Scenedesmus obliquus Mitochondrial Code")
            print("\t\t23. Thraustochytrium Mitochondrial Code")
            print("\t\t24. Pterobranchia Mitochondrial Code")
            print("\t\t25. Candidate Division SR1 and Gracilibacteria Code")
            print("\t\t26. Pachysolen tannophilus Nuclear Code")
            print("\t\t27. Karyorelict Nuclear Code")
            print("\t\t28. Condylostoma Nuclear Code")
            print("\t\t29. Mesodinium Nuclear Code")
            print("\t\t30. Peritrich Nuclear Code")
            print("\t\t31. Blastocrithidia Nuclear Code")
            print("\t\t33. Cephalodiscidae Mitochondrial UAA-Tyr Code\n")
            print("Additional details of translation tables can be found here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi\n\n")
            sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                inFasta = arg
            else:
                # error message
                print("\n\nThe specified input file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-t':
            if arg:
                TTable = arg
                #print("\nFile with taxa from clade 1: {}".format(clade1))
            else:
                # error message
                print("\n\nThe specified file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-s':
            if arg in ('T', 'F'):
                Stop = arg
                #print("\nFile with taxa from clade 1: {}".format(clade1))
            else:
                # error message
                print("\n\nBoolean argument to keep stop codons or not was not provided.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()        


    # pass to read_clades_into_list function
    translate(
        inFasta, TTable
        )

if __name__ == '__main__':
    main(sys.argv[1:])