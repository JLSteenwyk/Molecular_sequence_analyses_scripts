# Scripts to enable phylogenetic and related analyses

This repository houses numerous scripts to facilitate the analyses of phylogenetic trees and fasta files

## Scripts

A short description of each script

### Find_indel_unqiue_to_clade.py
Requires
```
python3
|- sys
|- getopt
|- os.path
|- re
|- Bio
|- numpy
```
Identifies indels unique to specific clades within an alignment fasta file.
Indels are identified using user specified step and window sizes.
Separate files for each clade should be provided that contains one taxa name
per line. <br />
For detailed information use the -h argument <br />
Basic usage: python Find_indel_unqiue_to_clade.py -w window -o clade1.file -t clade2.file -i alignment.fasta -s step

### Calculate_distance_between_two_taxa.py
Calculates phylogenetic distance between two taxa in a newick tree file.
Taxa names are provided as arguments. <br />
For detailed information use the -h argument <br />
Basic usage: python Calculate_distance_between_two_taxa.py -i newick_tree_file -o taxa1 -t taxa2

### calculate_average_protein_sequence_identity.bash
Calculates the average protein sequence identity between two whole genome
protein fasta files representing two different strains or species. This 
script implements a reciprocal best blast hit approach to determine average
protein sequence identity and implements scripts made available by [Harvard](http://archive.sysbio.harvard.edu/csb/resources/computational/scriptome/UNIX/Protocols/Sequences.html). <br />
Basic usage: bash calculate_average_protein_sequence_identity.bash A.pep.fasta B.pep.fasta

### calculate_clade1_clade2_branch_len.bash
Calculate the internode branch length of the branch that leads up to clade 1 and clade 2.
Input is a file of newick trees, which is the same as the input for ASTRAL coalescence based
phylogenetic inference. <br />
Variables Clade1, Clade2, and All_Clade are hardcoded and should be changed for each use.<br />
Basic usage: bash calculate_clade1_clade2_branch_len.bash file_of_newick_trees

## Authors

* **Jacob Steenwyk** - [Github page](https://jsteenwyk.github.io/)

## Acknowledgments

* Rokas lab personnel
* ACCRE

