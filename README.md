# Scripts to enable phylogenetic, molecular, and related analyses

This repository houses numerous scripts to facilitate the analyses of phylogenetic trees and fasta files

## Scripts

A short description of each script

### Find_indel_unqiue_to_clade.py
Identifies indels unique to specific clades within an alignment fasta file.
Indels are identified using user specified step and window sizes.
Separate files for each clade should be provided that contains one taxa name
per line. <br />
For detailed information use the -h argument <br />
```
python3
|- sys
|- getopt
|- os.path
|- re
|- Bio
   |- SeqIO
|- numpy
```
Basic usage: python Find_indel_unqiue_to_clade.py -w window -o clade1.file -t clade2.file -i alignment.fasta -s step <br />
Original author: [Jacob Steenwyk](https://jsteenwyk.github.io/)

### Find_indel_unqiue_to_clade_context_wOutgroup.py
Identifies indels unique to specific clades within an alignment fasta file.
Indels are identified using user specified step and window sizes.
Separate files for each clade should be provided that contains one taxa name
per line. Unlike 'Find_indel_unqiue_to_clade.py,' this script will contextualize
findings according to the outgroup status. That is to say that if clade 1 has gaps 
and clade 2 has nts and the outgroup clade has nts, clade 1 will be classified as 
a deletion. <br />
For detailed information use the -h argument <br />
```
python3
|- sys
|- getopt
|- os.path
|- re
|- Bio
   |- SeqIO
|- numpy
```
Basic usage: python Find_indel_unqiue_to_clade.py -w window -o clade1.file -t clade2.file -i alignment.fasta -s step -g outgroup_clade.file <br />
Original author: [Jacob Steenwyk](https://jsteenwyk.github.io/)

### Calculate_distance_between_two_taxa.py
Calculates phylogenetic distance between two taxa in a newick tree file.
Taxa names are provided as arguments. <br />
For detailed information use the -h argument <br />
```
python3
|- sys
|- getopt
|- os.path
|- re
|- Bio
   |- Phylo
      |- BaseTree
         |- TreeMixin
```
Basic usage: python Calculate_distance_between_two_taxa.py -i newick_tree_file -o taxa1 -t taxa2 <br />
Original author: [Jacob Steenwyk](https://jsteenwyk.github.io/)

### Calculate_pairwise_distances_among_taxa.py
Calculates all pairwise phylogenetic distances between two taxa in a newick tree file.
Taxa names are provided as a secondary file where taxa names are a single column. <br />
For detailed information use the -h argument <br />
```
python3
|- sys
|- getopt
|- os.path
|- re
|- Bio
   |- Phylo
      |- BaseTree
         |- TreeMixin
|- itertools
```
Basic usage: python Calculate_pairwise_distances_among_taxa.py -i newick_tree_file -l target_taxa_list <br />
Original author: [Jacob Steenwyk](https://jsteenwyk.github.io/)

### calculate_average_protein_sequence_identity.bash
Calculates the average protein sequence identity between two whole genome
protein fasta files representing two different strains or species. This
script implements a reciprocal best blast hit approach to determine average
protein sequence identity and implements scripts made available by [Harvard](http://archive.sysbio.harvard.edu/csb/resources/computational/scriptome/UNIX/Protocols/Sequences.html). <br />
```
perl v5.10.1
ncbi-blast-2.3.0+
|- makeblastdb
|- blastp
```
Basic usage: bash calculate_average_protein_sequence_identity.bash A.pep.fasta B.pep.fasta <br />
Original author: [Jacob Steenwyk](https://jsteenwyk.github.io/)

### calculate_clade1_clade2_branch_len.bash
Calculate the internode branch length of the branch that leads up to clade 1 and clade 2.
Input is a file of newick trees, which is the same as the input for ASTRAL coalescence based
phylogenetic inference. <br />
Variables Clade1, Clade2, and All_Clade are hardcoded and should be changed for each use.<br />
```
newick utilities v1.6
|- nw_clade
|- nw_labels
|- nw_distance
awk v3.1.7
```
Basic usage: bash calculate_clade1_clade2_branch_len.bash file_of_newick_trees <br />
Original author: [Jacob Steenwyk](https://jsteenwyk.github.io/)

### reorder_fasta.py
Reorder a multi-fasta file according to the order of header names in a third party file.
Input is the fasta file and a file with the header names in desired order of output
```
python3
|- sys
```
Basic usage: reorder_fasta.py file.fasta file.reference <br />
Original author: [dariober](http://seqanswers.com/forums/showthread.php?t=29558)

### seq_length.py
Prints fasta header and fasta length in a multi-fasta file.
```
python3
|- sys
|- Bio
   |- SeqIO
```
Basic usage: python seq_length.py mfasta.file <br />
Original author: [GummyBear](https://bioexpressblog.wordpress.com/2014/04/15/calculate-length-of-all-sequences-in-an-multi-fasta-file/)

### Split_Fasta_by_Header.sh
Renames headers in a fasta file to the following format >ID@xyz where 
xyz will be 0-n where n is the number of fasta headers minus one. <br /> 
First argument should be the ID and the second argument should be the multi-fasta file ending in ".fa"
```
awk v3.1.7
```
Basic usage: bash Split_Fasta_by_Header.sh Spp_ID fasta_file_name <br />
Original author: [Jacob Steenwyk](https://jsteenwyk.github.io/)

### busco2alignment.py
Takes the full table output of busco runs and creates a concatenated fasta file for phylogenetic inference using the concatenation method. Only argument is a configuration file specified with -c.
Configuration file should specify the path to [programs] mafft, trimAl, [input_files] a single column file with the names of the busco output directories, a single column file of the fasta files, [output_files] concatenation fasta file name, partition file output name, and [parameters] taxon occupancy with a value between 0 and 1. Example format is the following: <br /> 
[programs] <br />
mafft: pathway to mafft <br />
trimAl: pathway to trimAl <br /> 
<br />
[input_files] <br />
busco_out_list: busco_dirs.list <br />
fasta_files_list: fasta_files.list <br />
<br /> 
[output_files] <br />
concat_fasta: concat.fa <br />
partition_file: USCO_partition.txt <br /> 
<br />
[parameters] <br />
occupancy: value between 0 and 1 <br />
<br />
Exemplary template file can be printed out using the -t option
NOTE: busco_out and fasta_files should have files in the same order <br />
and the fasta file header names should be formatted in the following manner:<br /> 
\>indivID|1<br /> 
...<br /> 
\>indivID|2<br /> 
...<br /> 
where '|' is used in every gene name to denote specific genes but text prior to denotes a unique identifier for all genes from that fasta file.
For detailed information about usage, use the -h argument <br />
```
python3
|- sys
|- Bio
   |- SeqIO
   |- SeqRecord
      |- SeqRecord
   |- Alphabet
      |- IUPAC
   |- Seq
      |- Seq
|- getopt
|- os.path
|- os
|- re
|- subprocess
|- re
|- configparser
```
Basic usage: python busco2alignment.py -c config.busco2alignment<br />
Original author: [Jacob Steenwyk](https://jsteenwyk.github.io/)

This script is also broken up into two parts busco2uscofa.py and trimal2concat.py. <br />
Input and usage are the same but alignment and trimming must be done without these scripts. For trimal2concat.py, trimal files are assumed to have the following naming scheme buscoID.fa.mafft.trimal.

### busco_occupancy.py
Takes the full table output of busco runs and creates an occupancy matrix where columns are taxa and rows are busco IDs. The input file is a list of busco output dirs which is the same as "busco_out_list" in busco2alignment.py. Specifically, it is a single column file with busco output dir file names with one 
```
python3
|- sys
|- getopt
|- os.path
|- os
|- re
```
Basic usage: python busco_occupancy.py -b busco_output_dirs -o output_file<br />
Original author: [Jacob Steenwyk](https://jsteenwyk.github.io/)

### GC_content.py and AT_content.py
Takes a nucleotide fasta file as the first and only argument and determines the GC (or AT) content of the entirety of the file. This is used to determine genome GC (or AT) content.
```
python3
|- sys
|- getopt
|- os.path
|- math
|- re
```
Basic usage: python GC_content.py nucleotide.fasta or python AT_content.py nucleotide.fasta<br />
Original author: [Jacob Steenwyk](https://jsteenwyk.github.io/)

## Authors

* **Jacob Steenwyk** - [Github page](https://jsteenwyk.github.io/)
* The online community of bioinformaticians. Other/original authors are listed per script.

## Acknowledgments

* Rokas lab personnel
* ACCRE
* John Soghigian 
* Xingxing Shen

