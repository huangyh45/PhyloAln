#!/usr/bin/env perl

use Bio::SeqIO;
use Getopt::Std;
use Parallel::ForkManager;
use strict; 
my %opt=('g'=>'1','t'=>'*','n'=>'1','l'=>'transseq.log'); 
getopts('i:o:g:t:c:a:n:l:h',\%opt);
usage() if $opt{h};

my $ntfile=$opt{i};
my $aafile=$opt{o};
my $gencode=$opt{g};
my $termination=$opt{t};
my $incomplete=$opt{c};
my $transall=$opt{a};
our $numthreads=$opt{n};
our $logfile=$opt{l};
open(LOG,">>$logfile") or die "\nError: $logfile can't open!\n";

printdie("\nError: no input file or output file was set!\n") unless $ntfile&&$aafile;
printlog("\n##### Sequences translation begins #####\n");
my $in=Bio::SeqIO->new(-file => "$ntfile", -format => 'fasta');
my @seqs;
my $pm = new Parallel::ForkManager($numthreads); 
$pm -> run_on_finish( sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        if($data_structure_reference) {
            my @arr=@$data_structure_reference;
            if(@arr==1) {printdie($arr[0]);}
            else {
                push(@seqs, $arr[1]);
                printlog($arr[0]);
            }
        }
    });
while(my $seq=$in->next_seq()) {
    $pm->start and next;
    my $id=$seq->display_id;
    my @frames=(0,1,2,-1,-2,-3);
    my $finalseq;
    for(my $i=0;$i<@frames;$i++) {
        my ($transseq,$seq1);
        if($frames[$i]<0) {
        $seq1=$seq->revcom;
        $transseq=$seq1->translate(-codontable_id => $gencode, -terminator => $termination, -frame => -$frames[$i]-1);
        $seq1=$seq1->trunc(-$frames[$i],$seq1->length());
        }
        else {
        $transseq=$seq->translate(-codontable_id => $gencode, -terminator => $termination, -frame => $frames[$i]);
        $seq1=$seq->trunc($frames[$i]+1,$seq->length());
        }
        my $seq2=$transseq;
        if($incomplete&&$seq2->length()*3!=$seq1->length()) {
        if($seq2->length()*3+1==$seq1->length()||$seq2->length()*3+2==$seq1->length()) {
        $seq2=Bio::Seq->new( -seq => ($transseq->seq).'X' , -id  => $id );
        }
        else {my $diestr="\nError: nt length does not match the aa length in $id!\n";$pm->finish(0,[$diestr]);}
        }
        unless($transall) {$finalseq=$seq2;last;}
        else {
        my $start=$frames[$i]<0 ? $seq->length()+$frames[$i]+1 : $frames[$i]+1;
        my $end=$frames[$i]<0 ? $start%3+1 : $seq->length()-($seq->length()-$start+1)%3;
        if($incomplete) {$end=$frames[$i]<0 ? 1 : $seq->length();}
        $finalseq=Bio::Seq->new( -seq => $seq2->seq , -id  => "$id|$start-$end|$frames[$i]")
        }
    }
    $pm->finish(0,["Translation of $id finished.\n",$finalseq]);
}
$pm->wait_all_children;
printlog("#####i Sequences translation complished #####\n\n");
close(LOG);
my $out = Bio::SeqIO->new(-file => ">$aafile", -format => 'fasta');
foreach my $seq (@seqs) {$out->write_seq($seq);}

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
Translate nucleotide sequences in a file to amino acid sequences.

Usage: 
-i   input nucleotide sequences file
-o   output amino acid sequences file
-g   genetic code(default=1, invertebrate mitochondrion=5)
-t   symbol of termination(default='*')
-c   if translate incomplete codons into 'X'(default=no)
-a   if translate all six ORF(default=no)
-n   num threads(default=1)
-l   log file(default='transseq.log')
-h   this help message

Example:
transseq.pl -i ntfile -o aafile -g gencode -t termination -c 1 -a 1 -n numthreads -l logfile

Written by Yu-Hao Huang (2017-2024) huangyh45\@mail3.sysu.edu.cn
";

}
