#!/usr/bin/env perl

use Bio::SeqIO;
use Getopt::Std;
use strict;
my %opt; 
getopts('i:o:l:t:b:f:s:x:c:nh',\%opt);
usage() if $opt{h};
my $indir=$opt{i};
my $outfile=$opt{o} ? $opt{o} : 'all.fas';
my $type=$opt{t} ? $opt{t} : 'phyloaln';
my $fill=$opt{f} ? $opt{f} : '-';
my $sep=$opt{s} ? $opt{s} : '.';
my $suffix=$opt{x} ? $opt{x} : '.fa';
my $blockfile=$opt{b};
my $nexus=$opt{n};
my @codonpos=split('',$opt{c});
my $list=$opt{l} ? $opt{l} : makelist($indir,$type,$sep,$suffix);
my $suffixq=quotemeta $suffix;
my $sepq=quotemeta $sep;
open(F,"<$list") or die "\nError: $list can't open!\n";
open(FF,">$outfile") or die "\nError: $outfile can't open!\n";
if($blockfile) {
    open(B,">$blockfile") or die "\nError: $blockfile can't open!\n";
    if($nexus) {print B "#nexus\nbegin sets;\n";}
}
while(<F>) {
chomp;
my $taxon=$_;
my $str;
my $num=0;
while(<$indir/*$suffix>) {
    my $in=Bio::SeqIO->new(-file => "$_", -format => 'fasta');
    my $gene=$_;
    $gene=~s/(([\w\.]+)\/)+//;
    $gene=~s/$suffixq//;
    my ($find,$len);
    while(my $seq=$in->next_seq()) {
    my $id=$seq->display_id;
    $len=$seq->length() unless $len;
    if($type eq 'phyloaln'&&$id=~/([^$sepq]+)$sepq/) {$id=$1;}
    elsif($type eq 'blastsearch') {$id=~s/(\|(\S+))+//;$id=~s/\_$gene//;}
    elsif($type eq 'orthograph') {
        if($id=~/^(\w+)\|(\w+)\|/) {$gene=$1;$id=$2;}
        else {die "Error: Invalid format of $id in $_!\n";}
    }
    if($taxon eq $id) {$find=1;$str.=$seq->seq;last;}
    }
    unless($find) {
        printf "No $gene in $taxon!\n";
        for(my $i=0;$i<$len;$i++) {$str.=$fill;}
    }
    if($blockfile) {
        my $start=$num+1;
        $num+=$len;
	if(@codonpos) {
		for(my $i=0;$i<@codonpos;$i++) {
			if($nexus) {print B "charset ";}
			my $cstart=$start+$codonpos[$i]-1;
			print B "$gene\_codon$codonpos[$i]\t=\t$cstart-$num\\3;\n";
		}
	}
        else {
		if($nexus) {print B "charset ";}
		print B "$gene\t=\t$start-$num;\n";
	}
    }
}
$str=~s/\s//g;
print FF "\>$taxon\n$str\n";
if($nexus) {print B "end;\n";}
close(B);
$blockfile="";
}
close(F);

sub makelist {
	my $dir=$_[0];
	my $type=$_[1];
	my $sep=$_[2];
	my $suffix=$_[3];
	my %taxon;
	my $suffixq=quotemeta $suffix;
	my $sepq=quotemeta $sep;
	if($type eq 'phyloaln') {
		while(<$dir/*$suffix>) {
    		my $in=Bio::SeqIO->new(-file => "$_", -format => 'fasta');
		my $gene=$_;
   		$gene=~s/(([\w\.]+)\/)+//;
    		$gene=~s/$suffixq//;
    		while(my $seq=$in->next_seq()) {
    		my $id=$seq->display_id;
    		if($id=~/([^$sepq]+)$sepq/) {$id=$1;}
	    	$taxon{$id}=1;
    		}
		}
	}
	elsif($type eq 'blastsearch') {
		while(<$dir/*$suffix>) {
    		my $in=Bio::SeqIO->new(-file => "$_", -format => 'fasta');
		my $gene=$_;
   		$gene=~s/(([\w\.]+)\/)+//;
    		$gene=~s/$suffixq//;
    		while(my $seq=$in->next_seq()) {
    		my $id=$seq->display_id;
    		$id=~s/(\|(\S+))+//;
    		$id=~s/\_$gene//;
	    	$taxon{$id}=1;
    		}
		}
	}
	elsif($type eq 'orthograph') {
		while(<$dir/*$suffix>) {
    		my $in=Bio::SeqIO->new(-file => "$_", -format => 'fasta');
    		while(my $seq=$in->next_seq()) {
    		my $id=$seq->display_id;
    		if($id=~/^(\w+)\|(\w+)\|/) {my $taxa=$2;$taxon{$taxa}=1;}
		else {die "Error: Invalid format of $id in $_!\n";}
    		}
		}
	}
	my $list='list';
	open(L,">$list") or die "\nError: $list can't open!\n";
	foreach my $taxa (keys %taxon) {print L "$taxa\n";}
	close(L);
	return $list;
}

sub usage {

die "
perl $0
Concatenate multiple alignments into a matrix.

Usage: 
-i   directory containing input FASTA alignment files
-o   output concatenated FASTA alignment file
-t   type of input format(phyloaln/orthograph/blastsearch, default='phyloaln', also suitable for the format with same species name in all alignments, but the name shuold not contain separate symbol)
-f   the symbol to fill the sites of absent species in the alignments(default='-')
-s   the symbol to separate the sequences name and the first space is the species name in the 'phyloaln' format(default='.')
-x   the suffix of the input FASTA alignment files(default='.fa')
-b   the block file of the positions of each alignments(default=not to output)
-n   output the block file with NEXUS format, suitable for IQ-TREE(default=no)
-c   the codon positions to be written in the block file(default=no codon position, '123' represents outputing all the three codon positions, '12' represents outputing first and second positions)
-l   the list file with all the involved species you want to be included in the output alignments, one species per line(default=automatically generated, with all species found at least once in all the alignments)
-h   this help message

Example:
connect.pl -i inputdir -o outputfile -t inputtype -f fillsymbol -s separate -x suffix -b block1file -n -c codonpos -l listfile

Written by Yu-Hao Huang (2018-2024) huangyh45\@mail3.sysu.edu.cn
";

}

