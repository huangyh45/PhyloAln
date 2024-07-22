#!/usr/bin/env perl

use Bio::SeqIO;
use Bio::DB::Fasta;
use Getopt::Std;
use Parallel::ForkManager;
use strict; 
my %opt=('g'=>'1','t'=>'*','n'=>'1','l'=>'revertransseq.log'); 
getopts('i:b:o:g:t:n:l:h',\%opt);
usage() if $opt{h};

my @ntfiles=split(',',$opt{i});
my $aafile=$opt{b};
my $alignfile=$opt{o};
my $gencode=$opt{g};
my $termination=$opt{t};
our $numthreads=$opt{n};
our $logfile=$opt{l};
open(LOG,">>$logfile") or die "\nError: $logfile can't open!\n";

printdie("\nError: no input file, blueprint file or output file was set!\n") unless @ntfiles&&$aafile&&$alignfile;
printlog("\n##### Sequences reverse-translation begins #####\n");
my $aln=Bio::SeqIO->new(-file => "$aafile", -format => 'fasta');
my @seqs;
my $pm = new Parallel::ForkManager($numthreads); 
$pm -> run_on_finish( sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        my @arr=@$data_structure_reference;
        printlog($arr[0]);
        if(@arr>2) {push(@seqs,$arr[2]);}
        else {printdie($arr[1]);}
    });
while(my $seq=$aln->next_seq()) {
    $pm->start and next;
    my $outstr;
    my $seqid=$seq->display_id;
    my ($db,$ntseq);
    foreach my $ntfile (@ntfiles) {
        $db= Bio::DB::Fasta->new("$ntfile");
        $ntseq=$db->get_Seq_by_id($seqid);
        if($ntseq) {last;}
    }
    unless($ntseq) {$outstr.="\nWarning: $seqid can not be extracted in $opt{i}!\n";$pm->finish(0,[$outstr]);next;}
    my $seqstr;
    my $j=1;
    my $i;
    for($i=1;$i<=$seq->length();$i++) {
        my $base=$seq->subseq($i, $i);
        if($base eq '-') {$seqstr.="---";}
        elsif($j>$ntseq->length()) {$pm->finish(0,[$outstr,"\nError: no codon can be extracted in $seqid: $i-$base!"]);}
        else {
        my $n=$j+2>$ntseq->length() ? $ntseq->length() : $j+2;
        my $codon=$ntseq->subseq($j, $n);
        my $transbase=(Bio::Tools::CodonTable->new(-id=>$gencode))->translate($codon);
        $j=$n+1;
        if(length($codon)<3&&$base eq 'X') {$seqstr.=$codon;$outstr.="\nWarning: there is an incomplete codon in $seqid: $i-$base-$transbase($codon)\n";}
        elsif($base eq $transbase) {$seqstr.=$codon;}
        elsif($transbase eq '*') {
        if($base eq $termination) {$seqstr.=$codon;}
        else {$outstr.="\nWarning: there is an unexpected termination codon in $seqid: $i-$base-$transbase($codon)\n";$i-=1;}
        if($j<=$ntseq->length()) {$outstr.="\nWarning: there is a middle termination codon in $seqid: $i-$base-$transbase($codon)\n";}
        }
        else {$pm->finish(0,[$outstr,"\nError: the codon does not match the amino acid in $seqid: $i-$base-$transbase($codon)\n"]);}
        }
    }
    if($j+2==$ntseq->length()&&(Bio::Tools::CodonTable->new(-id=>$gencode))->is_ter_codon($ntseq->subseq($j, $j+2))) {
        my $codon=$ntseq->subseq($j, $j+2);
        $outstr.="\nWarning: there is an unexpected termination codon in $seqid: $i-end-*($codon)\n";
    }
    elsif($ntseq->length()==$j||$ntseq->length()==$j+1) {
        my $codon=$ntseq->subseq($j, $ntseq->length());
        $outstr.="\nWarning: there is an incomplete codon in $seqid: $i-end-($codon)\n";
    }
    elsif($ntseq->length()!=$j-1) {
        $pm->finish(0,[$outstr,"\nError: nt length does not match the aa length in $seqid!"]);
    }
    my $alignseq = Bio::Seq->new( -seq => $seqstr,
                                   -id  => $seqid,
                                 );
    $outstr.="Reverse-translation of $seqid finished.\n";
    $pm->finish(0,[$outstr,'',$alignseq]);
}
$pm->wait_all_children;
my $out = Bio::SeqIO->new(-file => ">$alignfile", -format => 'fasta');
foreach my $seq (@seqs) {$out->write_seq($seq);}
printlog("##### Sequences reverse-translation complished #####\n\n");
close(LOG);


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
Used the aligned translated sequences in a file as blueprint to aligned nucleotide sequences, which means reverse-translation.

Usage: 
-i   input nucleotide sequences file or files(separated by ',')
-b   aligned amino acid sequences file translated by input file as blueprint
-o   output aligned nucleotide sequences file
-g   genetic code(default=1, invertebrate mitochondrion=5)
-t   symbol of termination in blueprint(default='*')
-n   num threads(default=1)
-l   log file(default='revertransseq.log')
-h   this help message

Example:
revertransseq.pl -i ntfile1,ntfile2,ntfile3 -b aafile -o alignedfile -g gencode -t termination -n numthreads -l logfile

Written by Yu-Hao Huang (2017-2024) huangyh45\@mail3.sysu.edu.cn
";

}
