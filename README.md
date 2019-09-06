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

### find_positions_with_same_character_in_fasta.py
This script identifies positions in a fasta alignment file (-i parameter)
where for the taxa of interest (-t parameter), the script will determine the
positions that have the same character. For example, if position 1 in all taxa
have an 'A', then the script will report that position 1 has the same nucleotide
in all taxa..
```import sys
import getopt
import os.path
import re
from Bio import SeqIO
import numpy as np
python3
|- sys
|- getopt
|- os.path
|- re
|- Bio
   |- SeqIO
|- numpy
```
Basic usage: python remove_fasta_entries_wOnly_gaps.py -t taxa.list -i aligned_fasta.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### Determine_fixed_SNP_in_ingroup_and_no_where_in_outgroup.py
This script will identifies snps in one clade (-o parameter), which
can be one or more taxa. The script will ensure that the SNP is identified
by comparing the nucleotide or amino acid to an outgroup set of one or
more taxa (-g parameter). The output will detail the position of the SNP
and nucleotide or amino acid status in the ingroup and outgroup.
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
Basic usage: python Determine_fixed_SNP_in_ingroup_and_no_where_in_outgroup.py -o ingroup_taxa.list -g outgroup_taxa.list -i aligned_fasta.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### homopolymer_mutation_identifier.py
This script will identify homopolymer runs in the outgroup taxa sequence (-g parameter).
The script will then report substitutions if there are any in a taxon of interest (-t parameter).
Additionally the script will report any insertions or deletions in homopolymer sequences.
Warning: this script is written for CODON BASED ALIGNMENTS ONLY!
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
Basic usage: python homopolymer_mutation_identifier.py.py -t taxon -g outgroup_taxa.list -i aligned_fasta.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### CpG_to_TpG_identifier.py
This script will parse through a fasta file (-i parameter) and identify CpG -> TpG mutations suggestive of 5-mC mutations. To do so, it will first identify CpGs among outgroup taxa (-g parameter) and then compare the outgroup sequence to a taxon of interest (-t parameter). Output is in the following order: start, end, mut, CpGmut, tax, og
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- SeqIO
|- numpy
```
Basic usage: python CpG_to_TpG_identifier.py.py -i fasta.file -g outgroup_taxa.list -t taxon_of_interest <br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### Epsteins_coefficient_of_difference.py
Scans the outgroup (-g parameter) for conserved amino acids and then determines if there is an amino acid substitution in the taxon of interest (-t parameter). If there is a substitution, the Epsteins's coefficient of difference is determined. Values for Epstein's coefficient of difference values can be found here: https://en.wikipedia.org/wiki/Amino_acid_replacement#Epstein's_coefficient_of_difference Output is 5 columns that specify the position of the substition (cols 1 and 2), the taxon of interests sequence (col 3), the outgroup's sequence (col 4), and the Epsteins's index at that position (col 5).
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- SeqIO
   |- BiopythonWarning
|- Bio.Seq
   |- Seq
|- warnings
|- numpy
|- pprint
   |- pprint
```
Basic usage: python Epsteins_coefficient_of_difference.py -i fasta.file -g outgroup_taxa.list -t taxon_of_interest<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### Sneaths_index.py
Scans the outgroup (-g parameter) for conserved amino acids and then determines if there is an amino acid substitution in the taxon of interest (-t parameter). If there is a substitution, the Sneath's index value of dissimilarity is determined. Values for Sneath's index values can be found here: https://en.wikipedia.org/wiki/Amino_acid_replacement#Sneath's_index Output is 5 columns that specify the position of the substition (cols 1 and 2), the taxon of interests sequence (col 3), the outgroup's sequence (col 4), and the Sneath's index at that position (col 5).
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- SeqIO
   |- BiopythonWarning
|- Bio.Seq
   |- Seq
|- warnings
|- numpy
|- pprint
   |- pprint
```
Basic usage: python Sneaths_index.py -i fasta.file -g outgroup_taxa.list -t taxon_of_interest<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### determine_clusters_of_genes.py
Takes as input an NCBI features table and a bed file of the genes of interest.<br />
The script will then determine if the genes of interest are in a cluster or not.<br />
Additionally, this script will find genes inbetween clustered genes<br />
Output is a print out of gene clustering (or lack thereof). More specifically,<br />
col1: scaffolds genes are on<br />
col2: cluster start<br />
col3: cluster stop<br />
col4: number of unique homolog identifiers<br />
col5: number of total genes in the cluster<br />
col6: clustered gene identifiers in order of genomic appearance<br />
col7: homolog identifiers that correspond to the genes in col6<br />
The -b parameter file should contain: col1: scaffold, col2: gene start, col3: gene stop<br />
col4: gene ID, col5: variable column, and col6: the homolog identier. For example, for the MAL locus,<br />
the homolog identifier could be MALx3 for MALx3 homologs.<br />
For example, using the gliotoxin gene cluster of Aspergillus fumigatus Af293...<br />
CM000169.1    1434575    1436118    EAL88212.1    Aspergillus_fumigatus_Af293.GCA_000002655.1_ASM265v1    GliA<br />
CM000169.1    2271529    2272650    EAL90209.1    Aspergillus_fumigatus_Af293.GCA_000002655.1_ASM265v1    GliC<br />
CM000169.1    2678281    2693904    EAL90366.1    Aspergillus_fumigatus_Af293.GCA_000002655.1_ASM265v1    GliP<br />
...    ...    ...    ...        ...                             ...<br />
arguments
<br />-f    <feature table> -- as downloaded from NCBI
<br />-b    <bed file of genes of interest>
<br />-g    <gene distance boundary> -- An integer that represents how many genes away to look for linked clusters.
```
python3
|- sys
|- getopt
|- os.path
|- re
|- collections
   |- OrderedDict
|- numpy
```
Basic usage: python determine_clusters_of_genes.py -f feature.file -g gene_distance_boundary -b bed.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### fasta2paml.py
Converts a fasta file to a paml file. This was written explicitly for conducting dN/dS analysis but paml files can be used as input into various softwares.<br />
```
python3
|- sys
|- getopt
|- os.path
|- Bio 
   |- AlignIO
|- pprint
   |- pprint
```
Basic usage: python fasta2paml.py -i fasta.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### alignment_format_converter.py
This script takes an input fasta file and outputs a file in a different format to stdout. Choices of output include: clustal, emboss, fasta, fasta-m10, ig, maf, mauve, nexus, phylip, phylip-sequential, phylip-relaxed, and stockholm.<br />
Required arguments include:<br />
-i: input fasta file<br />
-f: output format<br />

```
python3
|- sys
|- getopt
|- os.path
|- Bio 
   |- AlignIO
```

### translate_DNA2Protein.py
This script will take as an input file a multi-fasta file and translate the nucleotides into their corresponding protein sequences.<br />
The input file is specified with the -i argument and the -t argument specifies the translation table.<br />
-i    fasta file<br />
input fasta file of nucleotide sequences<br />
-s    boolean for stop codon<br />
T or F for keeping or removing a stop codon in a sequence, respectively.<br />
-t    translation table<br /><br />
Translation tables are available through NCBI. Argument should be a number and not a string that describes the table.<br />
The following genetic codes are described here:<br />
    1. The Standard Code<br />
    2. The Vertebrate Mitochondrial Code<br />
    3. The Yeast Mitochondrial Code<br />
    4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code<br />
    5. The Invertebrate Mitochondrial Code<br />
    6. The Ciliate, Dasycladacean and Hexamita Nuclear Code<br />
    9. The Echinoderm and Flatworm Mitochondrial Code<br />
    10. The Euplotid Nuclear Code<br />
    11. The Bacterial, Archaeal and Plant Plastid Code<br />
    12. The Alternative Yeast Nuclear Code<br />
    13. The Ascidian Mitochondrial Code<br />
    14. The Alternative Flatworm Mitochondrial Code<br />
    16. Chlorophycean Mitochondrial Code<br />
    21. Trematode Mitochondrial Code<br />
    22. Scenedesmus obliquus Mitochondrial Code<br />
    23. Thraustochytrium Mitochondrial Code<br />
    24. Pterobranchia Mitochondrial Code<br />
    25. Candidate Division SR1 and Gracilibacteria Code<br />
    26. Pachysolen tannophilus Nuclear Code<br />
    27. Karyorelict Nuclear Code<br />
    28. Condylostoma Nuclear Code<br />
    29. Mesodinium Nuclear Code<br />
    30. Peritrich Nuclear Code<br />
    31. Blastocrithidia Nuclear Code<br />
    33. Cephalodiscidae Mitochondrial UAA-Tyr Code<br />
    Additional details of translation tables can be found here [https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)<br /><br />

```
python3
|- sys
|- getopt
|- os.path
|- Bio 
   |- SeqIO
|- Bio.Seq
   |- Seq
```

## Authors

* **Jacob Steenwyk** - [Github page](https://jlsteenwyk.github.io/)
* The online community of bioinformaticians. Other/original authors are listed per script.

## Acknowledgments

* [Rokas lab personnel](https://as.vanderbilt.edu/rokaslab/people/)
* [ACCRE](http://www.accre.vanderbilt.edu/)
* John Soghigian - [Github page](http://www.vector-eco-evo.com/)
* Xingxing Shen - [Github page](https://xingxingshen.github.io/)

