# Scripts to enable molecular sequence analyses

This repository houses numerous scripts to facilitate molecular sequence analyses

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
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

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
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

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
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

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
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

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
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### consensus_sequence_generator.py
Takes a multiple sequence alignment fasta file and a specified threshold value to create a consensus sequence among sequences. The threshold value should be between 0 and 1 representing 0-100% threshold for minimum AA present in a column to be considered the consensus.
```
python3
|- sys
|- getopt
|- os.path
|- os
|- Bio
   |- AlignIO
   |- Align
      |- AlignInfo
```
Basic usage: python consensus_sequence_generator.py -i protein.MSA.fasta -t threshold<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### remove_fasta_entries_wOnly_gaps.py
Takes a multi-fasta file and remove entries that contain only gaps ('-').
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- SeqIO
```
Basic usage: python remove_fasta_entries_wOnly_gaps.py -i fasta.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

## Authors

* **Jacob Steenwyk** - [Github page](https://jlsteenwyk.github.io/)
* The online community of bioinformaticians. Other/original authors are listed per script.

## Acknowledgments

* [Rokas lab personnel](https://as.vanderbilt.edu/rokaslab/people/)
* [ACCRE](http://www.accre.vanderbilt.edu/)
* John Soghigian - [Github page](http://www.vector-eco-evo.com/)
* Xingxing Shen - [Github page](https://xingxingshen.github.io/)

