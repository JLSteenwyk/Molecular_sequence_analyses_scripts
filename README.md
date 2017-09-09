# Scripts to enable phylogenetic and related analyses

This repository houses numerous scripts to facilitate the analyses of phylogenetic trees and fasta files

### Prerequisites

Required software

```
python3
```

## Scripts

A short description of each script

### Find_indel_unqiue_to_clade.py
Identifies indels unique to specific clades within an alignment fasta file.
Indels are identified using user specified step and window sizes.
Separate files for each clade should be provided that contains one taxa name
per line.__
For detailed information use the -h argument

### Calculate_distance_between_two_taxa.py
Calculates phylogenetic distance between two taxa is a newick tree file.
Taxa names are provided as arguments.__
For detailed information use the -h argument

### calculate_average_protein_sequence_identity.bash
Calculates the average protein sequence identity between two whole genome
protein fasta files representing two different strains or species. This 
script implements a reciprocal best blast hit approach to determine average
protein sequence identity and implements scripts made available by [Harvard](http://archive.sysbio.harvard.edu/csb/resources/computational/scriptome/UNIX/Protocols/Sequences.html).__
Basic usage is "bash calculate_average_protein_sequence_identity.bash A.pep.fasta B.pep.fasta"

## Authors

* **Jacob Steenwyk** - [Github page](https://jsteenwyk.github.io/)

## Acknowledgments

* Rokas lab personnel
* ACCRE

