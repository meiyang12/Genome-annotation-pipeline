#!/usr/bin/perl -w
use strict;
my $filename=$ARGV[0];
my $total_info;
my $timer;

open(GTH,"$filename") or die $!;
while(<GTH>){
	chomp();
	if($_=~/(\S+)\t(\S+)\t(gene)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\tID=gene(\S+);Target/){
		$timer++;
		$total_info.="$1\t$2\tmRNA\t$4\t$5\t$6\t$7\t$8\tID=gth$timer;\n";
	}elsif($_=~/(\S+)\t(\S+)\t(exon)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\tParent=gene(\S+)/){
		$total_info.="$1\t$2\tCDS\t$4\t$5\t$6\t$7\t$8\tParent=gth$timer;\n";
	}else{
		next;
	}
	}
close GTH;

print $total_info;
exit;
