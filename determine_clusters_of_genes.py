#!/usr/bin/env python

import sys
import getopt
import os.path
import re
import numpy as np
from collections import OrderedDict
import itertools

def create_and_print_clusters(
    clusters,
    bedGOIL
    ):
    """
    prints out final report of gene clusters

    Parameters
    ----------
    argv: clusters
        clusters for genes of interest and genes inbetween
    argv: bedGOIL
        list of lists of bed file for genes of interest
    """

    clusterArr = []

    # collapse clusters that share genes
    ## Looping over lines to collapse clusters with the same genes in them
    collapsed_clusters = []
    idx_ii             = 0
    idx_jj             = idx_ii + 1
    idx_opt            = True
    while_opt          = True
    # Looping over first line
    while(idx_opt):
        try:
            line_ii = clusters[idx_ii]
            line_jj = clusters[idx_jj]
            # Testing condition
            line_ii_opt = True
            # loop through clusters
            while(line_ii_opt):
                new_line = line_ii
                iiGenes = []
                jjGenes = []
                for iiEntries in line_ii:
                    iiGenes.append(iiEntries[3])
                for jjEntries in line_jj:
                    jjGenes.append(jjEntries[3])
                # if genes are shared between two cluster entries, collapse them into one entry
                if len(list(set(iiGenes).intersection(jjGenes))) >= 2:
                    idx_jj  += int(1)
                    # new line is line 1 start, line 2 stop, line 2 outgroup ID, line 2 gap ID
                    new_line = line_ii + line_jj
                    new_line.sort()
                    new_line = list(new_line for new_line,_ in itertools.groupby(new_line))
                    line_ii  = new_line
                    line_jj  = clusters[idx_jj]
                # if genes are not shared between two cluster entries do not collapse them into one entry
                else:
                    line_ii_opt = False
                    collapsed_clusters.append(new_line)
                    idx_ii = idx_jj
                    idx_jj = idx_ii + 1
        # except for end of clusters
        except IndexError:
            collapsed_clusters.append(line_ii)
            idx_opt = False

    for cluster in collapsed_clusters:
        # determine chr, start, end, and number of genes in cluster
        clusterChr     = cluster[0][0]
        clusterStart   = cluster[0][1]
        clusterGeneNum = len(cluster)
        clusterEnd     = cluster[clusterGeneNum-1][2]
        # determine the genes in the cluster
        clusterGenes     = []
        clusterGenesDict = {}
        for gene in cluster:
            clusterGenes.append(gene[3])
        for gene in clusterGenes:
            for line in bedGOIL:
                if gene == line[3]:
                    clusterGenesDict[gene] = line[5]
        # add non homolog cluster genes to clusterGenesDict
        for gene in clusterGenes:
            if gene in clusterGenesDict:
                1
            else:
                clusterGenesDict[gene] = 'NA'
        # get homolog information in same order as cluster genes
        clusterGeneHomID = []
        for gene in clusterGenes:
            clusterGeneHomID.append((clusterGenesDict[gene]))

        ## get number of distinct homologs -- do not count NAs
        clusterUniqHoms = []
        numberOfUniqHoms = 0
        # remove 'NA's
        clusterUniqHoms = [x for x in clusterGeneHomID if x != 'NA']
        # determine length of set
        numberOfUniqHoms = len(set(clusterUniqHoms))


        # append all pertinent cluster info to clusterArr
        clusterArrEntry = []
        clusterArrEntry.append(clusterChr)
        clusterArrEntry.append(clusterStart)
        clusterArrEntry.append(clusterEnd)
        clusterArrEntry.append(numberOfUniqHoms)
        clusterArrEntry.append(clusterGeneNum)
        genesStr = ';'.join(clusterGenes)
        clusterArrEntry.append(genesStr)
        geneHomIDStr = ';'.join(clusterGeneHomID)
        clusterArrEntry.append(geneHomIDStr)
        clusterArr.append(clusterArrEntry)
    for chro, start, stop, UniqHoms, geneNum, genes, homologInf in clusterArr:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chro, start, stop, UniqHoms, geneNum, genes, homologInf))


def obtain_genes_in_the_middle_of_cluster(
    near_entries_reduced,
    featuresL,
    bedGOIL
    ):
    """
    extracts genes in the middle of the cluster

    Parameters
    ----------
    argv: near_entries_reduced
        clustered genes of interest
    argv: featuresL
        NCBI features file
    argv: bedGOIL
        list of lists of bed file for genes of interest
    """

    # initialize variable to hold clusters
    clusters        = []
    featuresCluster = []

    #print(near_entries_reduced)
    # loop through clusters and find the genes inbetween multi-gene clusters
    for cluster in near_entries_reduced:
        if len(cluster) > 2:
            featureIndexStart = cluster[0][6]
            lastGeneIndex     = len(cluster)-1
            featureIndexEnd   = cluster[lastGeneIndex][6]+1
            # populate temp with all genes the cluster from the features table
            tempCluster = []
            tempCluster.append(featuresL[featureIndexStart:featureIndexEnd])
            for entry in tempCluster:
                clusterEntry = []
                for ent in entry:
                    temp = []
                    temp.append(ent[6])
                    temp.append(ent[7])
                    temp.append(ent[8])
                    temp.append(ent[10])
                    clusterEntry.append(temp) 
            clusters.append(clusterEntry)                  
        if len(cluster) == 2:
            featureIndexStart = cluster[0][6]
            lastGeneIndex     = len(cluster)-1
            featureIndexEnd   = cluster[lastGeneIndex][6]+1
            # populate temp with all genes the cluster from the features table
            tempCluster = []
            tempCluster.append(featuresL[featureIndexStart:featureIndexEnd])
            for entry in tempCluster:
                clusterEntry = []
                for ent in entry:
                    temp = []
                    temp.append(ent[6])
                    temp.append(ent[7])
                    temp.append(ent[8])
                    temp.append(ent[10])
                    clusterEntry.append(temp)   
            clusters.append(clusterEntry)
        if len(cluster) == 1:
            featureIndexStart = cluster[0][6]
            lastGeneIndex     = len(cluster)-1
            featureIndexEnd   = cluster[lastGeneIndex][6]+1
            # populate temp with all genes the cluster from the features table
            tempCluster = []
            tempCluster.append(featuresL[featureIndexStart:featureIndexEnd])
            for entry in tempCluster:
                clusterEntry = []
                for ent in entry:
                    temp = []
                    temp.append(ent[6])
                    temp.append(ent[7])
                    temp.append(ent[8])
                    temp.append(ent[10])
                    clusterEntry.append(temp)   
            clusters.append(clusterEntry)

    create_and_print_clusters(
    	clusters,
    	bedGOIL)

def create_cluster_profiles(
    featuresL,
    bedGOIL,
    geneDistBound
    ):
    """
    Read in input files into lists of lists

    Parameters
    ----------
    argv: featuresL
        list of lists of NCBI feature file
    argv: bedGOIL
        list of lists of bed file for genes of interest
    """

    near_entries = []

    # loop through bedGOIL
    for GOIentry in bedGOIL:
        # initialize variables
        gene = ''
        chro = ''
        homologID = ''
        
        # fill variables 
        gene = GOIentry[3]
        chro = GOIentry[0]
        homologID = GOIentry[5]

        # loop through features file
        for feature in featuresL:
            # if gene id match
            if feature[10] == gene:
                # create list entry ID
                Lentry = ''
                Lentry = feature[20]
                # obtain boundaries
                geneDistBound=int(geneDistBound)
                lowerLimit = Lentry-geneDistBound
                upperLimit = Lentry+geneDistBound+1
                nearFeatures = featuresL[lowerLimit:upperLimit]
                # remove nearFeatures entries on different scaffolds
                nearFeatures_singScaf = []
                for nearFeature in nearFeatures:
                	if nearFeature[6] == chro:
                	    nearFeatures_singScaf.append(nearFeature)

                # determine if any nearFeatures entries are near any GOIentries
                single_entry = []
                for nearFeature in nearFeatures_singScaf:
                    for GIentry in bedGOIL:
                        if nearFeature[10] == GIentry[3]:
                            tempEntry=[]
                            GIentry.append(nearFeature[20])
                            single_entry.append(GIentry)
                near_entries.append(single_entry)
    

    #print(near_entries)      
    # remove duplicate entries
    near_entries_reduced = []
    for sublist in near_entries:
        if sublist not in near_entries_reduced:
            near_entries_reduced.append(sublist)

    
    obtain_genes_in_the_middle_of_cluster(
        near_entries_reduced,
        featuresL,
        bedGOIL
        )


def read_inputs(
    feature,
    bedGOI,
    geneDistBound
    ):
    """
    Read in input files into lists of lists

    Parameters
    ----------
    argv: feature
        NCBI feature file
    argv: bedGOI
        bed file for genes of interest
    """

    # initialize lists
    featureL  = []
    featuresL = []
    bedGOIL   = []

    # read files into lists
    featureL = [line.rstrip('\n').split("\t") for line in open(feature)]
    bedGOIL  = [line.rstrip('\n').split("\t") for line in open(bedGOI)]

    cnt=0

    for line in featureL:
        if line[0] == 'CDS':
            line.append(cnt)
            featuresL.append(line)
            cnt+=1

    # pass to create_cluster_profiles function
    create_cluster_profiles(
        featuresL, 
        bedGOIL,
        geneDistBound
        )

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    feature       = ''
    geneDistBound = ''
    bedGOI        = ''


    try:
        opts, args = getopt.getopt(argv, "hf:b:g:")
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
            print("\nTakes as input an NCBI features table and a bed file of the genes of interest.")
            print("The script will then determine if the genes of interest are in a cluster or not.")
            print("Additionally, this script will find genes inbetween clustered genes.")
            print("Output is a print out of gene clustering (or lack thereof). More specifically,")
            print("col1: scaffolds genes are on")
            print("col2: cluster start")
            print("col3: cluster stop")
            print("col4: number of unique homolog identifiers")
            print("col5: number of total genes in the cluster")
            print("col6: clustered gene identifiers in order of genomic appearance")
            print("col7: homolog identifiers that correspond to the genes in col6\n")
            print("The -b parameter file should contain: col1: scaffold, col2: gene start, col3: gene stop")
            print("col4: gene ID, col5: variable column, and col6: the homolog identier. For example, for the MAL locus,")
            print("the homolog identifier could be MALx3 for MALx3 homologs.")
            print("For example, using the gliotoxin gene cluster of Aspergillus fumigatus Af293...")
            print("\nCM000169.1\t1434575\t1436118\tEAL88212.1\tAspergillus_fumigatus_Af293.GCA_000002655.1_ASM265v1\tGliA")
            print("CM000169.1\t2271529\t2272650\tEAL90209.1\tAspergillus_fumigatus_Af293.GCA_000002655.1_ASM265v1\tGliC")
            print("CM000169.1\t2678281\t2693904\tEAL90366.1\tAspergillus_fumigatus_Af293.GCA_000002655.1_ASM265v1\tGliP")
            print("...\t\t...\t...\t...\t\t...\t\t\t\t\t\t\t...\n")
            ## options
            # feature table
            print("\n-f\t<feature table>:")
            print("\tfeatures table from NCBI")
            # bed file
            print("\n-b\t<bed file of genes of interest>:")
            print("\tInput bed file.")
            # gene number boundary
            print("\n-g\t<gene distance boundary>:")
            print("\tAn integer that represents how many genes away to look for linked clusters.\n\n")
            sys.exit()

        elif opt == '-f':
            if arg:
                feature = arg
            else:
                # error message
                print("\n\nThe specified file of arg -f does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-g':
            if int(arg):
                geneDistBound = arg
            else:
                # error message
                print("\n\nThe specified input of arg -g must be an integer.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-b':
            if os.path.isfile(arg):
                bedGOI  = arg
            else:
                # error message
                print("\n\nThe specified file of arg -b does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to read inputs function
    read_inputs(
        feature,
        bedGOI,
        geneDistBound
        )

if __name__ == '__main__':
    main(sys.argv[1:])