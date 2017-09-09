# Scripts to enable phylogenetic and related analyses

This repository houses numerous scripts to facilitate the analyses of phylogenetic trees and fasta files

### Prerequisites

Required software
(assuming python3 is being used)

```
sys
getopt
os.path
re
Bio
numpy
```

how to install them

```
pip install "some package"
```

## Scripts

A short description of each script

**Find_indel_unqiue_to_clade.py**
Identifies indels unique to specific clades within an alignment fasta file.
Indels are identified using user specified step and window sizes.
Separate files for each clade should be provided that contains one taxa name
per line.
For detailed information use the -h argument

**Calculate_distance_between_two_taxa.py**
Calculates phylogenetic distance between two taxa is a newick tree file.
Taxa names are provided as arguments.
For detailed information use the -h argument

## Authors

* **Jacob Steenwyk** - [Github page](https://jsteenwyk.github.io/)

## Acknowledgments

* Rokas lab personnel
* ACCRE

