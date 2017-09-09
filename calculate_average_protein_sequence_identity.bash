# usage:
# bash calculate_average_protein_sequence_identity.bash A.fasta B.fasta
#
# explanation:
# this script calculates average protein sequence similarity using reciprocal 
# best blast hits between two genomes using a reciprocal best blast hit approach

### Step 00
## store arguments as variables
genomeA=$1
genomeB=$2

### Step 01
## Blast genome A versus genome B and vice versa

# create blast dbs for each genome
makeblastdb -in $genomeA -dbtype prot
makeblastdb -in $genomeB -dbtype prot

# blast genome A to B and vice versa
blastp -query $genomeA -db $genomeB -evalue 3 -outfmt 6 -num_threads 12 -out $genomeA.dbB.fmt6.blast
blastp -query $genomeB -db $genomeA -evalue 3 -outfmt 6 -num_threads 12 -out $genomeB.dbA.fmt6.blast

### Step 02
## execute RBBH analyses using scripts from harvard resource
## http://archive.sysbio.harvard.edu/csb/resources/computational/scriptome/UNIX/Protocols/Sequences.html
# 0
perl -e '$name_col=0; $score_col=11; while(<>) {s/\r?\n//; @F=split /\t/, $_; ($n, $s) = @F[$name_col, $score_col]; if (! exists($max{$n})) {push @names, $n}; if (! exists($max{$n}) || $s > $max{$n}) {$max{$n} = $s; $best{$n} = ()}; if ($s == $max{$n}) {$best{$n} .= "$_\n"};} for $n (@names) {print $best{$n}}' $genomeA.dbB.fmt6.blast > A_B.best
# 1
perl -e '$name_col=0; $score_col=11; while(<>) {s/\r?\n//; @F=split /\t/, $_; ($n, $s) = @F[$name_col, $score_col]; if (! exists($max{$n})) {push @names, $n}; if (! exists($max{$n}) || $s > $max{$n}) {$max{$n} = $s; $best{$n} = ()}; if ($s == $max{$n}) {$best{$n} .= "$_\n"};} for $n (@names) {print $best{$n}}' $genomeB.dbA.fmt6.blast > B_A.best
# 2
perl -e '$col1=1; $col2=0;' -e '($f1,$f2)=@ARGV; open(F1,$f1); while (<F1>) {s/\r?\n//; @F=split /\t/, $_; $line1{$F[$col1]} .= "$_\n"} open(F2,$f2); while (<F2>) {s/\r?\n//;@F=split /\t/, $_; if ($x = $line1{$F[$col2]}) {$x =~ s/\n/\t$_\n/g; print $x}}' A_B.best B_A.best > A_B_A
# 3
perl -e '$colm=0; $coln=13; $count=0; while(<>) {s/\r?\n//; @F=split /\t/, $_; if ($F[$colm] eq $F[$coln]) {print "$_\n"; $count++}}' A_B_A > A_B_A.recip
# 4 - modified to determine average protein similarity
perl -e '@cols=(0, 1, 2); while(<>) {s/\r?\n//; @F=split /\t/, $_; print join("\t", @F[@cols]), "\n"}' A_B_A.recip | awk '{ total += $3 } END { print total/NR }'

### Step 03
## clean
rm $genomeA.dbB.fmt6.blast
rm $genomeB.dbA.fmt6.blast
rm A_B.best
rm B_A.best
rm A_B_A
rm A_B_A.recip
rm $genomeA.phr
rm $genomeA.pin
rm $genomeA.psq
rm $genomeB.phr
rm $genomeB.pin
rm $genomeB.psq
