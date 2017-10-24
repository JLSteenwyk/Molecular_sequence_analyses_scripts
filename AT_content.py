#!/usr/bin/python

##############################################################
#							     #
# Author: Jacob Steenwyk				     #
# Date: April 24 2016					     #
# For description of script, append -h, --help or --usage    #
# for help or usage description				     #
#							     #
##############################################################

import sys, getopt, os.path, re, math

def calc_ATs(fastadata):
    """ calcs the occurence of ATs in fastadata """
    #initialize lists
    G_occurence = []
    A_occurence = []
    T_occurence = []
    C_occurence = []

    #count G|A|T|C
    G_occurence = fastadata.count('G')
    A_occurence = fastadata.count('A')
    T_occurence = fastadata.count('T')
    C_occurence = fastadata.count('C')
        
    GC = ''
    AT = ''
    total = ''

    GC = G_occurence + C_occurence
    AT = A_occurence + T_occurence
    total = GC + AT
    print(AT/total)    

def remove_star(fastadata):
    """ removes stop codons from fasta data list """
    fastadata = [x for x in fastadata if '*' not in x]
    calc_ATs(fastadata)
    
def AAs_list(ProteinFasta):
    """ Reads protein fasta file and creates list """
    fastadata = []
    with open(ProteinFasta, 'r') as file:
        file.readline()
        for line in file:
            line = line.rstrip('\n')
            if not '>' in line:
                for char in line:
                    fastadata.append(char)
    remove_star(fastadata)
    
""" Reads arguments and passes them to compare method """
ProteinFasta = '' # used for reading open sys.argv[1] file

try:	
	if sys.argv[1] == '-h':
		print('Calc AT content\nFirst and only argument should be an input nucleotide fasta file')
		sys.exit()
	elif sys.argv[1] == '--help':
		print('Calc AT content\nFirst and only argument should be an input nucleotide fasta file')
		sys.exit()
	else:
		if os.path.isfile(sys.argv[1]):
			ProteinFasta = sys.argv[1] #save argv 1
			AAs_list(ProteinFasta)
		else:
			print("\nThe specified file does not exist\nPlease check file pathway or file name\n\nFor help or usage message please append argument -h, --help or --usage\n\nThank you!\n")
except IndexError:
	print("\nNo file provided as argument\n\nFor help or usage message please append argument -h, --help or --usage\n\nThank you!\n")
