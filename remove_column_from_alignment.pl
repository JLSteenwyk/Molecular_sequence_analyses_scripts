#!/usr/bin/perl
use strict;
use warnings;

$ARGV[1] or die "use deleteAlignColumn.pl COLUMNS FASTA\n";
my @cols = split (/,/, $ARGV[0]);
@cols = sort {$b <=> $a} @cols;
open FA, "$ARGV[1]" or die;
$/ = "\n>";
while (<FA>) {
    s/>//g;
    my ($id, @seq) = split (/\n/, $_);
    my $seq = join "", @seq;
    foreach my $pos (@cols) {
        substr ($seq, $pos - 1, 1) = '';
    }
    print ">$id\n";
    while ($seq) {
         print substr($seq, 0, 80);
         print "\n";
         substr($seq, 0, 80) = '';
    }
}
close FA;
