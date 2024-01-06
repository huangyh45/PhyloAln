#!/usr/bin/env perl

use FindBin qw($Bin);
use Getopt::Std;
use strict; 
my %opt=('g'=>'1','a'=>'direct','t'=>'X','n'=>'1','l'=>'alignseq.log'); 
getopts('i:o:a:g:t:m:f:n:l:h',\%opt);
usage() if $opt{h};

my $input=$opt{i};
my $output=$opt{o};
my $aligntype=$opt{a};
my $gencode=$opt{g};
my $termination=$opt{t};
my $incomplete=$opt{c};
my $ifmediate=$opt{m};
my $mafftfolder=$opt{f};
if($mafftfolder) {$mafftfolder.="/";}
our $numthreads=$opt{n};
our $logfile=$opt{l};
open(LOG,">>$logfile") or die "\nError: $logfile can't open!\n";

printdie("\nError: no input file or output file was set!\n") unless $input&&$output;
printdie("\nError: invalid alignment type was set!\n") unless $aligntype eq 'direct'||$aligntype eq 'translate'||$aligntype eq 'codon'||$aligntype eq 'complement'||$aligntype eq 'ncRNA';
printlog("\n##### Alignment begins #####\n");
if($aligntype eq 'complement') {
    my $mafftcommand=$mafftfolder."mafft";
    runcmd("$mafftcommand --thread $numthreads --adjustdirectionaccurately $input > $output");
}
elsif($aligntype eq 'ncRNA') {
    my $mafftcommand=$mafftfolder."mafft-qinsi";
    runcmd("$mafftcommand --thread $numthreads $input > $output");
}
else {
    my ($mafftin,$mafftout,$transfile,$aafile);
    if($aligntype ne 'direct') {
        $transfile=$input;
        $transfile=~s/(\.(\w+))+/\.trans\.fas/;
        my %parameter=('-i',$input,'-o',$transfile,'-g',$gencode,'-t',$termination,'-c',$incomplete,'-n',$numthreads,'-l',$logfile);
        runpl("transseq.pl",\%parameter);
        $mafftin=$transfile;
        if($aligntype eq 'codon') {
        $aafile=$output;
        $aafile=~s/(\.(\w+))+/\.aa\.fas/;
        $mafftout=$aafile;
        }
        else {$mafftout=$output;}
    }
    else {$mafftin=$input;$mafftout=$output;}
    my $mafftcommand=$mafftfolder."linsi";
    runcmd("$mafftcommand --thread $numthreads $mafftin > $mafftout");
    if($aligntype eq 'codon') {
        my %parameter=('-i',$input,'-b',$aafile,'-o',$output,'-g',$gencode,'-t',$termination,'-n',$numthreads,'-l',$logfile);
        runpl("revertransseq.pl",\%parameter);
    }
    if($ifmediate) {
    if($transfile) {unlink("$transfile") or printdie("\nError: $transfile fail to delete!\n");}
    if($aafile) {unlink("$aafile") or printdie("\nError: $aafile fail to delete!\n");}
    }
}
printlog("##### Alignment complished #####\n\n");

close(LOG);

sub runcmd {
printlog("$_[0]\n");
my $iferr=system("$_[0] 2>$logfile.alignseq.temp");
open(TEMP,"<$logfile.alignseq.temp") or printdie("\nError: $logfile.alignseq.temp can't open!\n");
while(<TEMP>) {printlog("$_");}
close(TEMP);
unlink("$logfile.alignseq.temp") or printdie("\nError: $logfile.alignseq.temp fail to delete!\n");
if($iferr) {
    my $command=$_[0];
    if($command=~/(\S+)/o) {$command=$1;}
    printdie("\nError in $command!\n");
}
}

sub runpl {
my $pl=$_[0];
my $ref=$_[1];
my %parameter=%$ref;
my $command="";
foreach my $x (keys %parameter) {
    if($parameter{$x}) {
        if($parameter{$x}=~/\s/) {$parameter{$x}="\'$parameter{$x}\'";}
        $command.=" $x $parameter{$x}";
    }
}
printlog("$Bin/$pl $command\n");
if(system("$Bin/$pl $command")) {printdie("\nError in $pl!\n");}
}

sub printlog {
print LOG $_[0];
printf $_[0];
}

sub printdie {
print LOG $_[0];
die "$_[0]\nYou can use '-h' to watch detailed help.\n";
}

sub usage {

die "
perl $0
Align sequences in a file by mafft.
Requirement: mafft

Usage: 
-i   input sequences file
-o   output sequences file
-a   type of alignment(direct/translate/codon/complement/ncRNA, default='direct', 'translate' means alignment of translation of sequences)
-g   genetic code(default=1, invertebrate mitochondrion=5)
-t   symbol of termination(default='X', mafft will clean '*')
-c   if translate incomplete codons into 'X'(default=no)
-m   if delete the intermediate files, such as translated files and aligned aa files(default=no)
-f   the folder where mafft/linsi is, if mafft/linsi had been in PATH you can ignore this parameter
-n   num threads(default=1)
-l   log file(default='alignseq.log')
-h   this help message

Example:
alignseq.pl -i inputfile -o outputfile -a aligntype -g gencode -t termination -c 1 -m 1 -f mafftfolder -n numthreads -l logfile

Written by Yu-Hao Huang (2017-2023) huangyh45\@mail2.sysu.edu.cn
";

}
