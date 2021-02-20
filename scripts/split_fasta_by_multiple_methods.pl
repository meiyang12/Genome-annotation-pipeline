#! /usr/bin/perl -w  
=head1 Deception  
  
    split a fasta file to several pieces fastly.   
=head1 usage  
    perl split_fasta_by_multiple_methods.pl  <fa.file> [options]  
        -m , --method       the cut method,as follows:[1]  
                        1:keep every cut file having the same sequence numbers as far as possible;  
                        2.keep every cut file having the same size sequences as far as possible;  
                        3.cut fasta by assign sequece numbers in each cut file not the piece numbers.  
        -p , --piece        the piece number you want to get when you choose method 1 or 2.[40].  
        -c , --cut      the sequences numbers in the cut file,applyed to method 3.[100].  
        -d , --directory        output directory.[./].   
        -h , --help     screen the help information.  
  
=head1 Example  
    perl split_fasta_by_multiple_methods.pl  homo_sapiens.release-64.final.pep.fa  -m 1  -p 40  -d  ./out  
=cut  
  
use strict;  
use Getopt::Long;  
use File::Basename;  
my $method||=1;  
my $cut||=100;  
my $piece||=40;  
my $dir||="./";  
my $help;  
GetOptions(  
    "method|m:i"=>\$method,  
    "cut|c:i"=>\$cut,  
    "piece|p:i"=>\$piece,  
    "directory|d:s"=>\$dir,  
    "help|h"=>\$help,  
);  
  
  
die `pod2text $0` if(@ARGV!=1||$help);  
mkdir $dir,0755  if(! -e $dir);  
my $fa=shift;  
####################  
##### Method 1 #####  
####################  
if($method==1){  
    $fa=~/\.gz$/ ? (open IN,"gzip -cd $fa|"||die) : (open IN,$fa||die);  
    my %hash;  
    my $info;  
    map{chomp;if($_=~/^>(.*)/){$info=$1}else{$hash{$info}.="$_\n"}}<IN>;  
    close IN;  
  
    my $total=`less $fa|grep '>'|wc -l`;  
    my $count=int($total/$piece);  
    die "please assign the proper piece number!\n" unless $count;  
    for(my $i=1;$i<=$piece;$i++){  
        my $out_name=&named($fa,$i);  
        open OUT,'>', "$dir/$out_name" ||die;  
        my $flag=1;  
        while(my($k,$v)=each %hash){  
            print OUT ">$k\n$v";  
            last if($flag == $count && $i!=$piece);  
            $flag++;  
        }  
        close OUT;  
    }  
}  
####################  
##### Method 2 #####  
####################  
if($method==2){  
    $fa=~/\.gz$/ ? (open IN,"gzip -cd $fa|"||die) : (open IN,$fa||die);  
    my $total_len=0;  
    $/=">";<IN>;$/="\n";  
    while(<IN>){  
        $/=">";  
        my $seq=<IN>;  
        $/="\n";  
        $seq=~s/>|\n+//g;  
        $total_len+=length$seq;  
    }  
    close IN;  
    my $sub_len=int($total_len/$piece);  
    my $cur_len=0;  
        my $cur_content;  
        my $file_mark=1;  
    $fa=~/\.gz$/ ? (open IN,"gzip -cd $fa|"||die) : (open IN,$fa||die);  
    $/=">";<IN>;$/="\n";  
    while(<IN>){  
        chomp;  
        my $head=$_;  
        $/=">";  
        chomp(my $seq=<IN>);  
        $/="\n";  
        my $str=$seq;  
        $str=~s/\n+//g;  
        $cur_len+=length$str;  
        $cur_content.=">$head\n$seq";  
        if($cur_len >= $sub_len){  
            my $out_name=&named($fa,$file_mark);  
            open OUT,'>', "$dir/$out_name" ||die;  
            print OUT $cur_content;  
            close OUT;  
            $cur_len=0;  
            $cur_content="";  
            $file_mark++;  
        }  
    }  
    ##make the last file##  
    if($cur_len){  
        my $out_name=&named($fa,$file_mark);  
        open OUT,'>', "$dir/$out_name" ||die;  
                print OUT $cur_content;  
                close OUT;  
    }  
}  
close IN;  
####################  
##### Method 3 #####  
####################      
if($method==3){  
    $fa=~/\.gz$/ ? (open IN,"gzip -cd $fa|"||die) : (open IN,$fa||die);  
    my $i=1;  
    my $num=0;  
    my $out_name=&named($fa,$i);  
        open OUT,'>', "$dir/$out_name" ||die;  
    while(<IN>){  
                chomp;  
                if(/^>/){  
                        if($num < $cut){  
                                print OUT $_,"\n";  
                        }else{  
                                $i++;  
                                close OUT;  
                $out_name=&named($fa,$i);  
                open OUT,'>', "$dir/$out_name" ||die;  
                print OUT "$_\n";  
                $num=0;  
            }  
            $num++;  
        }else{  
            print OUT "$_\n";  
        }  
    }  
}  
close IN;  
close OUT;  
######subroutine######  
sub named{  
        my ($file,$num)=@_;  
    my $n=sprintf("%04g",$num);  
    my $name;  
    if($file=~/\.gz$/){  
                $file=~s/\.gz$//;  
                $name=basename($file).".$n";  
        }else{  
                $name=basename($file).".$n";  
        }  
    return $name;  
}  
__END__  