#!/usr/bin/perl
use strict;
my $i=0;
my $infile=shift;
my $flag=shift;
open IN,$infile || die $!;
while(<IN>){
chomp;
my @info=split(/\t/);
my $model=$info[2];
if ($model eq "mRNA"){
$i+=1;
}
if ($model eq "CDS" or $model eq "exon" ){
my $target=$1 if ($info[8]=~/Parent=(\S+);/);
print "$info[0]\t$flag\tnucleotide_to_protein_match\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\tID=match.$flag.$i;Target=$target\n";
}
}
close IN;
