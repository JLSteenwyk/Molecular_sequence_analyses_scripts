#!/usr/bin/env python

import sys
import getopt
import os.path
from datetime import datetime


def output(
    fastaD,
    fastaDfiltered,
    alleleCountD,
    alleleCountDFiltered,
    prefix,
    startTime
    ):
    """
    uses fastaD, fastaDfiltered, and allelecountD to create summary files
    
    Parameters
    ----------
    argv: fastaD
        dictionary with taxa keys and fasta sequences for each individual taxon
    argv: fastaDfiltered
        dictionary with taxa key and fasta seuqneces for each individual taxon that pass maf threshold
    argv: alleleCountD
        counts of each allele per site
    argv: alleleCountDFiltered
        counts of each allele per maf filtered site 
    argv: startTime
        startTime
    """

    print("Writing output files...")

    # create the names of the output files
    fastaFull         = prefix + ".all.fa"
    fastaFiltered     = prefix + ".mafFiltered.fa"
    partitionFull     = prefix + ".all.partition"
    partitionFiltered = prefix + ".mafFiltered.partition"
    
    # write out full fasta matrix to fastaFull
    with open(fastaFull, "a+") as output_handle:
        for k, v in fastaD.items():
            entry  = ''
            entry += ">" + k + "\n" + "".join(v) + "\n"
            output_handle.write(entry)
    # write out filtered fasta matrix to fastaFiltered
    with open(fastaFiltered, "a+") as output_handle:
        for k, v in fastaDfiltered.items():
            entry  = ''
            entry += ">" + k + "\n" + "".join(v) + "\n"
            output_handle.write(entry)

    # write out partition file with allele counts to partitionFull
    with open(partitionFull, "a+") as output_handle:
        cnt = 0
        for k, v in alleleCountD.items():
            # count ATCGs and determine frequency
            try:
                Acount = v['A']
            except KeyError:
                Acount = 0
            try:
                Tcount = v['T']
            except KeyError:
                Tcount = 0
            try:
                Ccount = v['C']
            except KeyError:
                Ccount = 0
            try:
                Gcount = v['G']
            except KeyError:
                Gcount = 0
            total  = Acount + Tcount + Ccount + Gcount
            Afreq  = Acount/total
            Tfreq  = Tcount/total
            Cfreq  = Ccount/total
            Gfreq  = Gcount/total
	        
            # create reports
            countReport = ','.join([str(Acount),str(Tcount),str(Ccount),str(Gcount)])
            freqReport  = ','.join([str(Afreq),str(Tfreq),str(Cfreq),str(Gfreq)])
            entry  = ''
            entry += str(cnt) + " " + k + " " + "A,T,C,G" + " " + countReport + " " + freqReport + "\n"
            output_handle.write(entry)
            cnt+=1

    # write out partition file with allele counts to partitionFiltered
    with open(partitionFiltered, "a+") as output_handle:
        cnt = 0
        for k, v in alleleCountDFiltered.items():
            # count ATCGs and determine frequency
            try:
                Acount = v['A']
            except KeyError:
                Acount = 0
            try:
                Tcount = v['T']
            except KeyError:
                Tcount = 0
            try:
                Ccount = v['C']
            except KeyError:
                Ccount = 0
            try:
                Gcount = v['G']
            except KeyError:
                Gcount = 0
            total  = Acount + Tcount + Ccount + Gcount
            Afreq  = Acount/total
            Tfreq  = Tcount/total
            Cfreq  = Ccount/total
            Gfreq  = Gcount/total
            
            # create reports
            countReport = ','.join([str(Acount),str(Tcount),str(Ccount),str(Gcount)])
            freqReport  = ','.join([str(Afreq),str(Tfreq),str(Cfreq),str(Gfreq)])
            entry  = ''
            entry += str(cnt) + " " + k + " " + "A,T,C,G" + " " + countReport + " " + freqReport + "\n"
            output_handle.write(entry)
            cnt+=1

    # print log message
    print("Complete!\n")
    print("Total time:", datetime.now() - startTime)
    print("")


def summarize(
    inFileD,
    referenceL,
    maf,
    prefix,
    fastaD,
    ref,
    fastaDfiltered,
    startTime
    ):
    """
    uses inFileD and referenceL to create summary fasta and partition files
    
    Parameters
    ----------
    argv: inFileD
        dictionary with taxa keys and values are lists of lists of variable
        sites in a given taxon. In the values, the information is stored
        in the following format:
        - ele 0 is scaffold in the reference
        - ele 1 is the position in the reference
        - ele 2 is query character (i.e., the snp)
    argv: referenceL
        across all unique sites in inFileD, the scaffold, position, and reference
        character is provided in this list of lists 
    argv: maf
        minor allele frequency
    argv: prefix
        output prefix names
    argv: fastaD
        dictionary with taxa keys and empty values. This will eventually be populated
        with the fasta sequences for each individual taxon
    argv: ref
        reference string
    argv: fastaDfiltered
        dictionary with taxa key and empty values. This will eventually be populated 
        with the fasta seuqneces for each individual taxon that pass maf threshold
    argv: startTime
        startTime
    """

    # create a new dictionary to hold the bare bones of the fasta file
    # i.e., the key will be the taxon and the value will be the fasta entry

    # sort list of lists by first and then second 
    referenceL.sort(key=lambda k: (k[0], k[int(1)]), reverse = False)

    # create an empty dictionary with sites as keys and all alleles as values for 
    # use to calcualte minor allele frequency and filtering using maf argument
    allelesPerSite = {}
    for entry in referenceL:
        # create a string that will be used for keeping track of a site
        entryStr = ''
        entryStr = entry[0] + ' ' + entry[1]
        # make the site into a key with an empty value in allelesPerSite dictionary
        allelesPerSite[entryStr] = []
  
    # loop through entries in the reference list
    for entry in referenceL:
        # create a string that will be used to look for matching keys in inFileD items 
        entryStr = ''
        entryStr = entry[0] + ' ' + entry[1]
        # loop through taxon and snps in inFileD
        for taxon, snps in inFileD.items():
        	# loop for any snps that have a matching position as in the reference list
            if entryStr in snps.keys():
                # if there is a snp, append the identity of the snp to the value of the taxon key 
                fastaD[taxon].append(snps[entryStr])
                # additionally, append the identity of the snp to the value of the variable site
                allelesPerSite[entryStr].append(snps[entryStr])
            else:
                # if there is no snp, append the identity of the refence character to the value of the taxon key 
                fastaD[taxon].append(entry[2])
                # additionally, append the identity of the snp to the value of the nonvariable site
                allelesPerSite[entryStr].append(entry[2])
        # append the reference character to the reference
        fastaD[ref].append(entry[2])  
        # additionally, append the identity of the character in the reference
        allelesPerSite[entryStr].append(entry[2])  

    ## determine what sites do not pass maf threshold
    # create an empty dictionary for alleles counts
    alleleCountD={}
    # loop through allele counts and count the occurence of each allele
    for k, v in allelesPerSite.items():
    	alleleCountD[k]={x:v.count(x) for x in v}
    # loop through alleleCountD and remove sites that do not pass maf threshold. 
    # store sites that pass this threshold in initialized list mafFiltered
    mafFiltered = []
    values = []
    alleleCountDFiltered = {}
    # loop through alleleCountD
    for k, v in alleleCountD.items():
        # determine the observed maf
        totalTaxa = len(inFileD)
        values=list(v.values())
        length = len(values)
        values.sort()
        mafSite = values[length-2]
        mafSite = mafSite/totalTaxa
        # if site does not pass maf, then filter out the site
        if mafSite > maf:
            mafFiltered.append(k)
            alleleCountDFiltered[k] = v

    ## create a maf filtered reference list
    # initialize maf filtered reference list
    mafFilteredRef = []
    # loop through maf filtered sites and get reference list entry
    for site in mafFiltered:
        site = site.split(' ')
        for entry in referenceL:
            tempEle = []
            tempEle.append(entry[0])
            tempEle.append(entry[1])
            if site == tempEle:
                mafFilteredRef.append(entry)

    # loop through entries in the reference list
    for entry in mafFilteredRef:
        # create a string that will be used to look for matching keys in inFileD items 
        entryStr = ''
        entryStr = entry[0] + ' ' + entry[1]
        # loop through taxon and snps in inFileD
        for taxon, snps in inFileD.items():
        	# loop for any snps that have a matching position as in the reference list
            if entryStr in snps.keys():
                # if there is a snp, append the identity of the snp to the value of the taxon key 
                fastaDfiltered[taxon].append(snps[entryStr])
                # additionally, append the identity of the snp to the value of the variable site
                allelesPerSite[entryStr].append(snps[entryStr])
            else:
                # if there is no snp, append the identity of the refence character to the value of the taxon key 
                fastaDfiltered[taxon].append(entry[2])
                # additionally, append the identity of the snp to the value of the nonvariable site
                allelesPerSite[entryStr].append(entry[2])
        # append the reference character to the reference
        fastaDfiltered[ref].append(entry[2])  
        # additionally, append the identity of the character in the reference
        allelesPerSite[entryStr].append(entry[2])  

    # print an update to the user about the length of the unfiltered and filtered alignment lengths
    print()
    for k, v in fastaD.items():
        print("Unfiltered alignment length:", len(v))
        break
    for k, v in fastaDfiltered.items():
        print("MAF filtered alignment length:", len(v))
        break
    print()

    output(
        fastaD,
        fastaDfiltered,
        alleleCountD,
        alleleCountDFiltered,
        prefix, 
        startTime
        )



def read(
	maf, taxa, 
    prefix, ref,
    inFiles, startTime
    ):
    """
    Reads in input files and saves them as lists of lists.
    Then determine snp sites in query sequences and reference
    
    Parameters
    ----------
    argv: maf
        minor allele frequency
    argv: taxa
        list of taxa
    argv: prefix
        output prefix names
    argv: ref
        reference taxon name
    argv: inFile
        input snps files from mummer
    argv: startTime
        startTime
    """

    # initialize lists and dictionaries
    taxaL       = []
    inFileL     = []
    inFileD     = {}
    referenceL  = []
    fastaD      = {}
    fastaDfiltered = {} 


    # read lists into python lists
    taxaL   = taxa.split(',')
    inFileL = inFiles.split(',')

    # create keys in dictionary. Each key has a taxa. The values will be populated with the file contents
    for i in taxaL: 
        inFileD[i] = {}
    for i in taxaL: 
        fastaD[i] = []
    fastaD[ref] = []
    for i in taxaL: 
        fastaDfiltered[i] = []
    fastaDfiltered[ref] = []    
    
    ## create a dictionary with taxa keys and snp file values. Also keep track of reference positions
    ## and append reference positions and character identities to a separate list of lists 
    # loop through files and taxa lists
    for file, taxa in zip(inFileL, taxaL):
        with open(file) as file:
            # skip first four lines because they are mummer header information
            for _ in range(4):
                next(file)
            for line in file:
                # strip new line characters and split at tabs
                line = line.rstrip('\n')
                line = line.split('\t')

                ## populate tempList with important formation from the snps file
                # temp list is for each taxa
                tempStrKey   = []
                # ele 0 is scaffold in the reference
                tempStrKey = line[8] + ' ' + line[0]
                # append tempList to inFileD
                inFileD[taxa][tempStrKey] = line[2]

                ## populate refList with important information from the snps file 
                # refList is for the reference list of lists
                refList = []
                # ele 0 is the scaffold in the reference
                refList.append(line[8])
                # ele 1 is the position in the reference
                refList.append(line[0])
                # ele 2 is the reference character
                refList.append(line[1])
                # append refList to referenceL if scaffold, position, and character haven't been stored yet 
                if refList in referenceL:
                    1
                else:
                    referenceL.append(refList)

    # pass to the summarize function to create summary files
    summarize(
        inFileD,
        referenceL,
        maf,
        prefix,
        fastaD,
        ref,
        fastaDfiltered,
        startTime
        )

def main(
    argv
    ):
    """
    Reads arguments 
    """
    # start time of script
    startTime = datetime.now()

    # initialize argument variables
    maf     = 0
    taxa    = ''
    prefix  = ''
    ref     = ''
    inFiles = ''

    try:
        opts, args = getopt.getopt(argv, "hi:t:p:r:m:")
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
            print("")
            print("This script will summarize the information across multiple MUMMER snps files.")
            print("The summary will be presented as a mfasta file summarizing the SNPs, a partition")
            print("file showing where the snps are found, and the frequency of nucleotides at each position.")
            print("")
            print("The following input files are required:")
            print("-i: a list of input snp files separated by a comma (e.g., snp1.file,snp2.file,snp3.file...)")
            print("-t: a list of input snp files separated by a comma (e.g., taxa1,taxa2,taxa3...)")
            print("\tThese must be in the same as the order of files in the -i file. These names will be used in the mfasta.")
            print("-r: the name of the reference taxon name that will be used in the mfasta file")
            print("-p: prefix of the output file names that will be made")
            print("-m: minor allele frequency (MAF) to be used in the creation of the mfasta file (default: 0)")
            print("")
            print("This script will create four output files:")
            print("(1) .all.fa contains a fasta file for all taxa with all positions")
            print("(2) .all.partition contains a summary of where the snp is found and allele frequencies for every site\n")
            print("Input MUMMER script files were obtained using the following 2-step pipeline: ")
            print("(1) > ./nucmer --prefix=$prefix $reference $query")
            print("(2) > ./show-snps -CTIr $delta > $snps.")
            print("This has been tested using the output from MUMMER version 4.0.0beta2.")
            print("")
            sys.exit()
        # read in maf cut-off 
        elif opt == '-m':
            print(arg)
            if 0 <= float(arg) <= 1:
                maf = float(arg)
            else:
                maf = 0
        # read in file of taxa
        elif opt == '-t':
            if arg:
                taxa = arg
            else:
                # error message
                print("\n\nSpecify -t argument.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-p':
            if arg:
                prefix = str(arg)
            else:
                # error message
                print("\n\nSpecify -p argument.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-r':
            if arg:
                ref = str(arg)
            else:
                # error message
                print("\n\nSpecify -r argument.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-i':
            if arg:
                inFiles = arg
            else:
                # error message
                print("\n\nThe specified -i file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()        

    # print log output
    print("\n"+"-"*(len("- Argument parameters -")))
    print("| Argument parameters |")
    print("-"*(len("- Argument parameters -")))
    print("-m", maf)
    print("-t", taxa)
    print("-p", prefix)
    print("-r", ref)
    print("-i", inFiles)

    # pass to read function to read in input files
    read(
        maf, taxa, 
        prefix, ref,
        inFiles, startTime
        )

if __name__ == '__main__':
    main(sys.argv[1:])