#!/usr/bin/perl
use strict;
open IN,@ARGV[0] || die $!;
my $g=1;
my $e=1;
while (<IN>){
chomp;
my @info=split(/\t/);
my $chr=$info[0];
my $predict=$info[1];
my $model=$info[2];
if ($model eq "mRNA" ){
my $id=$1 if ($info[8]=~/ID=(\S+\d);+/);
my $gene="$chr".".g$g";
print "$chr\t$predict\tgene\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\tID=$gene;Name=$predict model $id\n";
print "$chr\t$predict\t$model\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\tID=$id;Parent=$gene\n";
$g++;
}
if ($model eq "CDS"){
my $Parent=$1 if ($info[8]=~/Parent=(\S+\d);+/);
my $i="$chr".".e$e";
print "$chr\t$predict\texon\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\tID=$i;Parent=$Parent\n";
print "$chr\t$predict\t$model\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\tID=cds.of.$Parent;Parent=$Parent\n";
$e++;
}
}
close IN;
