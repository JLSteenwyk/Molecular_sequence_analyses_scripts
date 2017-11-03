#!/usr/bin/env python

import sys
import getopt
import os.path
import os
import subprocess
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import configparser
from collections import Counter

def create_fas_per_USCO(
    orthoMCLout, fastaFiles,
    taxon_occupancy, mafft_path,
    trimAl_path, concat_fa,
    SCOG_partition, geneIds,
    fasta_dir, SCOGs
    ):
    """
    creates fasta files per SCOG that passes the 
    taxon occupancy criterion

    Parameters
    ----------
    argv: orthoMCLout
        orthoMCL output
    argv: fastaFiles
        list of fasta files
    argv: taxon_occupancy
        value between 0 and 1 for taxon occupancy per USCO
    argv: mafft_path
        pathway to mafft program
    argv: trimAl_path
        pathway to trimAl program
    argv: concat_fa
        output file name for concat fa
    argv: SCOG_partition
        output file name for the partition file
    argv: geneIds
        list with gene ids
    argv: fasta_dir
        directory with fasta files
    """

    # initialize variable for fastaFiles len
    fastaFilesLength = ''
    # create empty files to populate fasta seqs to
    for SCOGid in SCOGs:
        SCOGfasta = SCOGid+".fa"
        open(SCOGfasta, "w")

    # loop through list of USCOs that pass taxon occupancy
    fastaFilesLength = len(fastaFiles)
    for fasta in fastaFiles:
       
        # open busco output table
        fastaFileName = fasta_dir+fasta
        # check if tableFileName exists and exit if not
        if os.path.isfile(fastaFileName):
            1
        else:
            print(fastaFileName)
            print("fasta file does not exist. fasta_directory should end in \"/\" character")
            sys.exit()


def determine_SCOGs(
    orthoMCLout, fastaFiles,
    taxon_occupancy, mafft_path,
    trimAl_path, concat_fa,
    SCOG_partition, geneIds,
    fasta_dir
    ):
    """
    determines single copy orthologous groups
    
    Parameters
    ----------
    argv: orthoMCLout
        orthoMCL output
    argv: fastaFiles
        list of fasta files
    argv: taxon_occupancy
        value between 0 and 1 for taxon occupancy per USCO
    argv: mafft_path
        pathway to mafft program
    argv: trimAl_path
        pathway to trimAl program
    argv: concat_fa
        output file name for concat fa
    argv: SCOG_partition
        output file name for the partition file
    argv: geneIds
        list with gene ids
    argv: fasta_dir
        directory with fasta files
    """

    # initialize counts dictionary to count number of genes per OG
    counts = {}
    SCOGsD = {}
    # and OG or SCOG to hold OG ID or SCOG ids
    OG     = []
    SCOGs  = []

    # set variables for determining if OG is SC
    numTaxaInAnalysis = len(geneIds)

    # open orthoMCL output
    with open(orthoMCLout, "r") as f:
        # loop through lines in file
        for line in f:
            # remove text inside parantheses ie. '(.*):'
            line = re.sub(r'\([^)]*\):', '', line)
            # extract OG identity
            OG = (line.split()[0])
            
            temp = []
            tempD = {}
            for indiv in geneIds:
                if line.count(indiv) == 1:
                    temp.append(indiv) 

            if len(temp)/numTaxaInAnalysis >= float(taxon_occupancy):
                SCOGs.append(OG)
                line = re.split('\(|\)| ',line.rstrip('\n'))
                if len([e for e in line if e.startswith(indiv)]) == 1:
                    SCOGsD[OG] = [e for e in line if e.startswith(indiv)]
                    print(SCOGsD)



            # # modify line to extract everything within parantheses
            # regex = re.compile(".*?\((.*?)\)")
            # line = re.findall(regex, line)
            # # counts of each string occurence in list
            # counts = dict((x,line.count(x)) for x in set(line))
            
            # # determine number of total genes in OG
            # total_genes = sum(counts.values())
            # # determine number of indivs in OG
            # number_of_indivs = len(counts)
            
            # # if SCOG, save SCOG ID. otherwise, skip
            # num_taxa_single_copy = []
            # for k, v in counts.items():
            #     if v == 1:
            #         num_taxa_single_copy.append(k)
            # if len(num_taxa_single_copy) >= (numTaxaInAnalysis*float(taxon_occupancy)):
            #     SCOGs.append(OG)


    # create_fas_per_USCO(
    #     orthoMCLout, fastaFiles,
    #     taxon_occupancy, mafft_path,
    #     trimAl_path, concat_fa,
    #     SCOG_partition, geneIds,
    #     fasta_dir, SCOGs
    #     )

    ## SCOGs are single copy OGs that need to be aligned, trimmed, and concated

def read_lists(
    orthoMCLout, fasta_list,
    taxon_occupancy, mafft_path,
    trimAl_path, concat_fa,
    SCOG_partition, gene_ids,
    fasta_dir
	):
    """
    Reads fasta and gene ids list into a python list 

    Parameters
    ----------
    argv: orthoMCLout
        orthoMCL output
    argv: fasta_list
        single column file of fasta files
    argv: taxon_occupancy
        value between 0 and 1 for taxon occupancy per USCO
    argv: mafft_path
        pathway to mafft program
    argv: trimAl_path
        pathway to trimAl program
    argv: concat_fa
        output file name for concat fa
    argv: SCOG_partition
        output file name for the partition file
    argv: gene_ids
        file that contains the gene ids in a single col
    argv: fasta_dir
        directory with fasta files
    """

    # initialize lists
    fastaFiles = []
    geneIds    = []

    # read lists into python lists
    fastaFiles = [line.rstrip('\n') for line in open(fasta_list)]
    geneIds = [line.rstrip('\n') for line in open(gene_ids)]

    determine_SCOGs(
        orthoMCLout, fastaFiles,
        taxon_occupancy, mafft_path,
        trimAl_path, concat_fa,
        SCOG_partition, geneIds,
        fasta_dir
        )

def read_config(
    config_file
    ):
    """
    Reads config file

    Parameters
    ----------
    argv: configuration file
        provides all necessary parameters to run
    """
    orthoMCLout     = ''
    fasta_list      = ''
    fasta_dir       = ''
    taxon_occupancy = ''
    mafft_path      = ''
    trimAl_path     = ''
    concat_fa       = ''
    SCOG_partition  = ''
    gene_ids        = ''

    settings = configparser.ConfigParser()
    settings.read(config_file)

    # read and check busco_out path
    orthoMCLout = settings.get('input_files', 'orthoMCL_out_file')
    if os.path.isfile(orthoMCLout):
        1
    else:
        print("No or incorrect pathway to orthoMCL_out_file provided in configuration file!\nExiting now...")    
        sys.exit()

    # read and check fasta_files path
    fasta_list = settings.get('input_files', 'fasta_files_list')
    if os.path.isfile(fasta_list):
        1
    else:
        print("No or incorrect pathway to fasta_files provided in configuration file!\nExiting now...")    
        sys.exit()

    # read and checks fasta dir path
    fasta_dir = settings.get('input_files', 'fasta_directory')
    if os.path.exists(fasta_dir):
        1
    else:
        print("No or incorrect pathway to fasta_files provided in configuration file!\nExiting now...")
        sys.exit()

    # read and check taxon_occupancy is between 0 and 1
    taxon_occupancy = settings.get('parameters', 'occupancy')
    if 0 <= float(taxon_occupancy) <= 1:
        1
    else:
        print("No or incorrect value for occupancy parameter!\nExiting now...")    
        sys.exit()

    # read and check mafft path
    mafft_path = settings.get('programs', 'mafft')
    if os.path.isfile(mafft_path):
        1
    else:
        print("No or incorrect pathway to mafft provided in configuration file!\nExiting now...")    
        sys.exit()

    # read and check trimal path
    trimAl_path = settings.get('programs', 'trimAl')
    if os.path.isfile(trimAl_path):
        1
    else:
        print("No or incorrect pathway to trimal provided in configuration file!\nExiting now...")    
        sys.exit()

    # read output concatenation file path
    concat_fa = settings.get('output_files', 'concat_fasta')
    if concat_fa:
        1
    else:
        print("No concat_fasta specified in configuration file!\nExiting now...")    
        sys.exit()

    # read and check partition_file path
    SCOG_partition = settings.get('output_files', 'partition_file')
    if SCOG_partition:
        1
    else:
        print("No partition_file specified in configuration file!\nExiting now...")    
        sys.exit()

    # read gene name ids path
    gene_ids = settings.get('input_files', 'gene_ids_file')
    if gene_ids:
        1
    else:
        print("No gene_ids_file specified in configuration file!\nExiting now...")    
        sys.exit()


    # pass to read_lists function
    read_lists(
        orthoMCLout, fasta_list,
        taxon_occupancy, mafft_path,
        trimAl_path, concat_fa,
        SCOG_partition, gene_ids,
        fasta_dir
        )


def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    config_file = ''

    try:
        opts, args = getopt.getopt(argv, "hc:t")
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
            # options
            print("\n\n-c configuration file")
            print("\tfile that contains all necessary parameters and file names")
            print("\n-t print configuration template")
            print("\tprints a template configuration file for the user to fill out")
            print("\n")
            print("Contents of configuration file")
            print("------------------------------")
            ## programs
            # mafft
            print("mafft: <pathway to mafft program>")
            print("\tpathway to mafft program. Version(s) tested: 7.294b")
            # trimal
            print("\ntrimAl: <pathway to trimAl program>")
            print("\tpathway to trimal program. Version(s) tested: 1.4.rev11")
            ## input_files explanation
            # busco_out_list
            print("\torthoMCL_out_file: <orthomcl output file>")
            print("\tsingle column file of busco output dirs")
            print("\torder should match the list of fasta file")
            print("\tspecified in fasta_files parameter")
            # fasta_files_list explanation
            print("\nfasta_files_list: <list of fasta files>")
            print("\tsingle column file of fasta files")
            print("\torder should match the list of busco dirs")
            print("\tspecified in busco_out parameter")
            print("\tNOTE: fasta file headers per fa file")
            print("\tshould be named in the following convention...")
            print("\t>ID@1")
            print("\t>ID@2")
            print("\t>ID@3")
            print("\t> ...")
            print("\t>ID@N")
            # fasta_directory explained
            print("\nfasta_directory: <pathway to fasta file directory>")
            print("\tpathway to the directory that contains the fasta files")
            # gene_ids_file explained
            print("\ngene_ids_file: <format of gene IDs")
            print("\tcontains the formatting for how genes are named. The")
            print("\torder of gene names needs to appear in the same file at the")
            print("\tfasta_files_list entries.")
            ## parameters 
            # occupancy explanation
            print("\noccupancy: <taxon occupancy>")
            print("\ttaxon occupancy for orthologous groups")
            print("\tvalue should be between 0 and 1\n")
            sys.exit()
        elif opt == '-c':
            if os.path.isfile(arg):
                config_file = arg
            else:
                # error message
                print("\n\nThe specified configuration file does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-t':
            print("\n[programs]")
            print("mafft: pathway to mafft")
            print("trimAl: pathway to trimAl")
            print("\n[input_files]")
            print("orthoMCL_out_file: orthomcl output file")
            print("fasta_files_list: single column file with names of fasta files")
            print("fasta_directory: directory with fasta files, #should end in \"/\"")
            print("\n[output_files]")
            print("concat_fasta: output concatenation fasta file name")
            print("partition_file: RAxML partition file")
            print("\n[parameters]")
            print("occupancy: taxon occupancy per orthologous group, a value between 0 and 1\n")
            sys.exit()

    # pass to read_config parameter
    read_config(
        config_file
        )

if __name__ == '__main__':
    main(sys.argv[1:])