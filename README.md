# PhyloAln
## PhyloAln: a convenient reference-based tool to align sequences and high-throughput reads for phylogeny and evolution in the omic era

![logo](https://github.com/huangyh45/PhyloAln/blob/main/logo.png)  

PhyloAln is a reference-based alignment tool for phylogeny and evolution. PhyloAln can directly map not only the raw reads but also assembled or translated sequences into the reference alignments, which is suitable for different omic data and skips the complex preparation in the traditional method including data assembly, gene prediction, orthology assignment and sequence alignment, with a relatively high accuracy in the aligned sites and the downstream phylogeny. It is also able to detect and remove foreign and cross contamination in the generated alignments, which is not considered in other reference-based methods, and thus improve the quality of the alignments for downstream analyses.

### Catalogue
- [Installation](#installation)
  - [1) Installation from source](#1-installation-from-source)
  - [2) Installation using Conda](#2-installation-using-conda)
- [Usage](#usage)
  - [Quick start](#quick-start)
  - [A practice using PhyloAln for phylogenomics](#a-practice-using-phyloaln-for-phylogenomics)
  - [Input](#input)
  - [Output](#output)
  - [Example commands for different data](#example-commands-for-different-data)
  - [Detailed parameters](#detailed-parameters)
  - [Limitations](#limitations)
- [Auxiliary scripts for PhyloAln and phylogenetic analyses](#auxiliary-scripts-for-phyloaln-and-phylogenetic-analyses)
  - [Script to translate sequences: transseq.pl](#transseqpl)
  - [Script to back-translate sequences: revertransseq.pl](#revertransseqpl)
  - [Script to align the sequences: alignseq.pl](#alignseqpl)
  - [Script to concatenate the alignments: connect.pl](#connectpl)
  - [Script to combine each result alignment of different PhyloAln runs: merge_seqs.py](#merge_seqspy)
  - [Script to select the sequences in the files in bulk: select_seqs.py](#select_seqspy)
  - [Script to trim the alignments based on unknown sites in bulk: trim_matrix.py](#trim_matrixpy)
  - [Script to root the phylognetic tree: root_tree.py](#root_treepy)
  - [Script to test performance of PhyloAln: test_effect.py](#test_effectpy)
- [Citation](#citation)
- [Questions & Answers](#questions--answers)
  - [How can I obtain the reference alignments and the final tree?](#how-can-i-obtain-the-reference-alignments-and-the-final-tree)
  - [Does selection of the outgroup influence detection of foreign contamination? How can I choose an appropriate outgroup?](#does-selection-of-the-outgroup-influence-detection-of-foreign-contamination-how-can-i-choose-an-appropriate-outgroup)
  - [The required memory is too large to run PhyloAln.](#the-required-memory-is-too-large-to-run-PhyloAln)
  - [The positions of sites in the reference alignments are changed in the output alignments.](#the-positions-of-sites-in-the-reference-alignments-are-changed-in-the-output-alignments)
  - [How can I assemble the paired-end reads?](#how-can-i-assemble-the-paired-end-reads)
  - [Can PhyloAln generate the alignments of multiple-copy genes for gene family analyses?](#can-phyloaln-generate-the-alignments-of-multiple-copy-genes-for-gene-family-analyses)
  
### Installation

#### 1) Installation from source
##### Requirements
- python >=3.7.4 (https://www.python.org/downloads/)
- biopython >=1.77 (https://biopython.org/wiki/Download)
- hmmer >=3.1 (http://hmmer.org/download.html)
- mafft >=7.467 (optional for the auxiliary scripts, https://mafft.cbrc.jp/alignment/software/source.html)
- ete3 >=3.1.2 (optional for the auxiliary scripts, http://etetoolkit.org/download/)
- perl >=5.26.2 (optional for the auxiliary scripts, https://www.perl.org/get.html)
- perl-bioperl >=1.7.2 (optional for the auxiliary scripts, https://github.com/bioperl/bioperl-live/blob/master/README.md)
- perl-parallel-forkmanager >=2.02 (optional for the auxiliary scripts, https://github.com/dluxhu/perl-parallel-forkmanager)

After installing these requirements, you can download PhyloAln from this GitHub repo directly from this page or using the command in your computer:  
```
git clone https://github.com/huangyh45/PhyloAln.git
```
If your computer needs execute permissions to run the programs, such as the Linux or macOS system, you should first run the command :  
```
chmod -R +x /your/PhyloAln/path/ 
```
Then, you can test if PhyloAln has been available using the commands:   
```
cd /your/PhyloAln/path/  
export PATH=$PATH:/your/PhyloAln/path/:/your/PhyloAln/path/scripts  
bash tests/run_test.sh
```
When you see "Successfully installed" at the end of the screen output, PhyloAln and all its auxiliary scripts has been successfully installed and available.  
If the test fails, you should check if the requirements have been successfully installed and executable in the current environment.  
After test, you can manually delete all the newly generated files, or run the command to delete them:  
```
rm -rf alignseq.log all.block all.fas list tests/run_test.config tests/PhyloAln_* tests/aln tests/ref/*.fas tests/ref/*.index
```
##### Full usage experience
If you have installed [IQ-TREE](http://www.iqtree.org/#download) and want to experience the usage of all the scripts with real examples through a simple phylogenomic flow after installation, you can run the commands (it will spend some minutes):  
```
bash tests/run_test.sh full
```
After running, you can manually delete all the newly generated files, or run the command to delete them:  
```
rm -rf alignseq.log all.block all.fas list tests/run_test.config tests/PhyloAln_* tests/aln tests/ref/*.fas tests/ref/*.index
```
#### 2) Installation using Conda
Download the [Conda configure file of requirements](https://github.com/huangyh45/PhyloAln/releases/download/v0.1.0/requirement.txt), and install the requirements using the command:  
```
conda install --file requirement.txt
```
If the installation spends too much time, you can try to install the requirements and all their dependencies with fixed but not latest version. Download the [Conda configure file of requirements with fixed version](https://github.com/huangyh45/PhyloAln/releases/download/v0.1.0/requirement_fix.txt), and install these requirements using the command:  
```
conda install --file requirement_fix.txt
```
Then, you can manually install PhyloAln through the above steps in [1) Installation from source](#1-installation-from-source), or download [Conda package of PhyloAln](https://github.com/huangyh45/PhyloAln/releases/download/v0.1.0/phyloaln-0.1.0-py312_0.conda.tar.bz2) and run the command:  
```
conda install phyloaln-0.1.0-py312_0.conda.tar.bz2
```
Note: This is a temporary Conda installation method. The conda package of PhyloAln on public channel is preparing and coming soon.

### Usage

#### Quick start
If you have only one reference alignment FASTA file and sequence data from only one source/species, you can use -a to input the reference alignment file, -s to input the species name and -i to input the FASTA/FASTQ sequence/read file(s), like this command:  
```
PhyloAln -a reference_alignment_file -s species -i sequence_file1 (sequence_file2) -o output_directory
```

You can also use -c to input a configure file representing information of sequence data from multiple sources/species. The configure file should be tab-separated and like this:  
```
species1  /absolute/path/sequence_file1
species2  /absolute/path/sequence_file1,/absolute/path/sequence_file2
```
If you have a directory containing multiple reference alignment FASTA files with a same suffix, you can use -d to input the directory and -x to input the suffix. The command using multiple reference alignments and multiple sources/species is like this:  
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory
```

#### A practice using PhyloAln for phylogenomics
The following practice is for phylogenomics using codon alignments of nuclear single-copy orthologous groups and 20 CPUs.
##### 1. obtain the reference orthologous sequences
You can download the reference sequences from the ortholog database (e.g., [OrthoDB](https://www.orthodb.org/), [OMA](https://omabrowser.org/oma/home/)), or perform *de novo* orthology assignment (e.g., by [OrthoFinder](https://github.com/davidemms/OrthoFinder)). The reference species should contain an outgroup for PhyloAln.
##### 2. codon alignment for each ortholog group
In this step, you can use our auxiliary script [alignseq.pl](#alignseqpl)  
Run the shell commands:  
```
mkdir aln  
for file in orthogroup/*.fa; do  
  name=`basename $file`  
  scripts/alignseq.pl -i $file -o aln/$name -a codon -n 20  
done
```
##### 3. trim the alignments
In this step, you can use the tool [trimAl](https://github.com/inab/trimal)  
Run the shell commands:  
```
mkdir ref_aln  
for file in aln/*.fa; do  
  name=`basename $file`  
  trimal -in $file -out ref_aln/$name -automated1 -keepheader  
done
```
The reference alignments have been generated here. And you can directly obtain the existing alignments as reference instead of the above three steps, for example, from the published supplementary data.
##### 4. write the configure of the species and data
The format of the configure file is TSV and like this:  
```
species1  /absolute/path/sequence_file1  
species2  /absolute/path/sequence_file1,/absolute/path/sequence_file2
```
##### 5. run PhyloAln the map the sequences/reads into the reference alignments
```
PhyloAln -d ref_aln -c config.tsv -p 20 -m codon -u outgroup
```
The output alignments can be trimmed to remove the sites with too many unknown bases using our auxiliary script [trim_matrix.py](#trim_matrixpy)
##### 6. concatenate the alignments into a supermatrix
This step can be done with our auxiliary script [connect.pl](#connectpl)  
For codon dataset, you can run it like:  
```
scripts/connect.pl -i PhyloAln_out/nt_out -f N -b all.block -n -c 123
```
For protein dataset, the command is like:  
```
scripts/connect.pl -i PhyloAln_out/aa_out -f X -b all.block -n
```
##### 7. reconstruct the phylogenetic tree
You can build the tree by [IQ-TREE](http://www.iqtree.org/#download) like this:  
```
iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 20 --prefix species_tree
```
##### 8. root the tree
You can root the tree with the outgroup using our auxiliary script [root_tree.py](#root_treepy)
```
scripts/root_tree.py species_tree.treefile species_tree.rooted.tre outgroup
```
Finally you obtain a species tree with NEWICK format here and you can then visualize it or use it in other downstream analyses.

#### Input
PhyloAln needs two types of file:  
- the alignment file(s) with FASTA format. Trimmed alignments with conservative sites are recommended. Multiple alignment files with the same suffix should be placed into a directory for input.
- the sequence/read file(s) with FASTA or FASTQ format. Compressed files ending with ".gz" are allowed. Sequence/read files from multiple sources/species should be inputed through a configure file as described in quick start.

#### Output
PhyloAln generates new alignment file(s) with FASTA format. Each output alignment is corresponding to each reference alignment file, with the aligned target sequences from the provided sequence/read file(s). If using prot, codon or dna_codon mode, the translated protein alignments will be also generated. These alignments are mainly for phylogenetic analyses and evolutionary analyses using conservative sites.

#### Example commands for different data
Notice: the following commands are only recommended according to our practice, and you can modify the options as you need.

Map the DNA sequences/reads into the DNA alignments:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -u outgroup
```
Map the protein sequences into the protein alignments:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -u outgroup -n
```
Map the DNA sequences/reads into the protein alignments:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m prot -u outgroup
```
Map the DNA sequences/reads into the codon alignments:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m codon -u outgroup
```
Map the DNA sequences/reads into large numbers of codon alignments:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m codon -u outgroup -b
```
Map the DNA assembly sequences into the codon alignments:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m codon -u outgroup -b -r
```
Map the long reads with high insertion and deletetion rates into the codon alignments (actually not recommended to use long reads with high error rates):
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m dna_codon -u outgroup
```
Map the genomic sequences/reads with intron regions into the codon alignments:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m codon -u outgroup -l 200 -f large_fasta
```
Map the sequences/reads into the codon alignments using the non-standard genetic code (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for detail), for example, the codon alignments of plastid protein-coding genes:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m codon -u outgroup -g 11
```
Map the sequences/reads into the concatenated or other long DNA alignments:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -u outgroup ---ref_split_len 1000
```

#### Detailed parameters
```
usage: PhyloAln [options] -a reference_alignment_file -s species -i fasta_file -f fasta -o output_directory  
PhyloAln [options] -d reference_alignments_directory -c config.tsv -f fastq -o output_directory  
  
A program to directly generate alignments from FASTA/FASTQ files based on reference alignments for phylogenetic analyses.  
  
optional arguments:  
  -h, --help            show this help message and exit  
  -a ALN, --aln ALN     the single reference FASTA alignment file  
  -d ALN_DIR, --aln_dir ALN_DIR  
                        the directory containing all the reference FASTA  
                        alignment files  
  -x ALN_SUFFIX, --aln_suffix ALN_SUFFIX  
                        the suffix of the reference FASTA alignment files when  
                        using "-d"(default:.fa)  
  -s SPECIES, --species SPECIES  
                        the studied species ID for the provided FASTA/FASTQ  
                        files(-i)  
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]  
                        the input FASTA/FASTQ file(s) of the single  
                        species(-s), compressed files ending with ".gz" are  
                        allowed  
  -c CONFIG, --config CONFIG  
                        the TSV file with the format of 'species  
                        sequence_file(s)(absolute path, files separated by  
                        commas)' per line for multiple species  
  -f {guess,fastq,fasta,large_fasta}, --file_format {guess,fastq,fasta,large_fasta}  
                        the file format of the provided FASTA/FASTQ files,  
                        'large_fasta' is recommended for speeding up reading  
                        the FASTA files with long sequences(e.g. genome  
                        sequences) and cannot be guessed(default:guess)  
  -o OUTPUT, --output OUTPUT  
                        the output directory containing the  
                        results(default:PhyloAln_out)  
  -p CPU, --cpu CPU     maximum threads to be totally used in parallel  
                        tasks(default:8)  
  --parallel PARALLEL   number of parallel tasks for each alignments, number  
                        of CPUs used for single alignment will be  
                        automatically calculated by '--cpu /  
                        --parallel'(default:the smaller value between number  
                        of alignments and the maximum threads to be used)  
  -m {dna,prot,codon,dna_codon}, --mol_type {dna,prot,codon,dna_codon}  
                        the molecular type of the reference  
                        alignments(default:dna, 'dna' suitable for nucleotide-  
                        to-nucleotide or protein-to-protein alignment, 'prot'  
                        suitable for protein-to-nucleotide alignment, 'codon'  
                        and 'dna_codon' suitable for codon-to-nucleotide  
                        alignment based on protein and nucleotide alignments  
                        respectively)  
  -g GENCODE, --gencode GENCODE  
                        the genetic code used in translation(default:1 = the  
                        standard code, see https://www.ncbi.nlm.nih.gov/Taxono  
                        my/Utils/wprintgc.cgi)  
  --ref_split_len REF_SPLIT_LEN  
                        If provided, split the reference alignments longer  
                        than this length into short alignments with this  
                        length, ~1000 may be recommended for concatenated  
                        alignments, and codon alignments should be devided by  
                        3  
  -l SPLIT_LEN, --split_len SPLIT_LEN  
                        If provided, split the sequences longer than this  
                        length into short sequences with this length, 200 may  
                        be recommended for long genomic reads or sequences  
  --split_slide SPLIT_SLIDE  
                        the slide to split the sequences using sliding window  
                        method(default:half of '--split_len')  
  -n, --no_reverse      not to prepare and search the reverse strand of the  
                        sequences, recommended for searching protein or CDS  
                        sequences  
  --low_mem             use a low-memory but slower mode to prepare the reads,  
                        'large_fasta' format is not supported and gz  
                        compressed files may still spend some memory  
  --hmmbuild_parameters HMMBUILD_PARAMETERS [HMMBUILD_PARAMETERS ...]  
                        the parameters when using HMMER hmmbuild for reference  
                        preparation(default:[])  
  --hmmsearch_parameters HMMSEARCH_PARAMETERS [HMMSEARCH_PARAMETERS ...]  
                        the parameters when using HMMER hmmsearch for mapping  
                        the sequences(default:[])  
  -b, --no_assemble     not to assemble the raw sequences based on overlap  
                        regions  
  --overlap_len OVERLAP_LEN  
                        minimum overlap length when assembling the raw  
                        sequences(default:30)  
  --overlap_pident OVERLAP_PIDENT  
                        minimum overlap percent identity when assembling the  
                        raw sequences(default:98.00)  
  -t, --no_out_filter   not to filter the foreign or no-signal sequences based  
                        on conservative score  
  -u OUTGROUP, --outgroup OUTGROUP  
                        the outgroup species for foreign or no-signal  
                        sequences detection(default:the first sequence in each  
                        alignment)  
  -q SEP, --sep SEP     the separate symbol between species name and gene  
                        identifier in the sequence headers of the  
                        alignments(default:.)  
  --outgroup_weight OUTGROUP_WEIGHT  
                        the weight coefficient to adjust strictness of the  
                        foreign or no-signal sequence filter, small number or  
                        decimal means ralaxed criterion (default:0.90, 1 = not  
                        adjust)  
  -r, --no_cross_species  
                        not to remove the cross contamination for multiple  
                        species  
  --cross_overlap_len CROSS_OVERLAP_LEN  
                        minimum overlap length when cross contamination  
                        detection(default:30)  
  --cross_overlap_pident CROSS_OVERLAP_PIDENT  
                        minimum overlap percent identity when cross  
                        contamination detection(default:98.00)  
  --min_exp MIN_EXP     minimum expression value when cross contamination  
                        detection(default:0.20)  
  --min_exp_fold MIN_EXP_FOLD  
                        minimum expression fold when cross contamination  
                        detection(default:2.00)  
  -w UNKNOW_SYMBOL, --unknow_symbol UNKNOW_SYMBOL  
                        the symbol representing unknown bases for missing  
                        regions(default:unknow = 'N' in nucleotide alignments  
                        and 'X' in protein alignments)  
  -z {consensus,consensus_strict,all,expression,length}, --final_seq {consensus,consensus_strict,all,expression,length}  
                        the mode to output the sequences(default:consensus,  
                        'consensus' means selecting most common bases from all  
                        sequences, 'consensus_strict' means only selecting the  
                        common bases and remaining the different bases unknow,  
                        'all' means remaining all sequences, 'expression'  
                        means the sequence with highest read counts after  
                        assembly, 'length' means sequence with longest length  
  -y, --no_ref          not to output the reference sequences  
  -v, --version         show program's version number and exit  
  
Written by Yu-Hao Huang (2023-2024) huangyh45@mail3.sysu.edu.cn
```

#### Limitations
- PhyloAln is only designed for phylogenetic analyses and evolutionary analyses with reference-based conservative sites, and thus cannot perform *de novo* assembly due to non-conservative sites and sites not covered in the reference alignments. The unmapped sites will be ignored.
- We prioritize the flexibility of PhyloAln and thus did not provide the upstream steps of collecting the reference sequences and generating the reference alignments, and the downstream phylogenetic analyses. But you can use the auxiliary scripts to help preparation and perform downstream analyses.
- In current version, we did not heavily focus on optimizing the runtime and memory usage of PhyloAln. Specially, speed and memory usage may be influenced by the numbers of the reference alignments and the numbers of the target sequences/reads. A version with optimized parallel and storage operations and optional accessories using C or other rapid languages may be developed in the future. Faster sequence search tools may also be a candidate to be integrated into PhyloAln as an option to speed up the alignments.

###  Auxiliary scripts for PhyloAln and phylogenetic analyses

#### transseq.pl
Requirements:  
- perl >=5.26.2
- perl-bioperl >=1.7.2
- perl-parallel-forkmanager >=2.02

```
perl scripts/transseq.pl  
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
  
Written by Yu-Hao Huang (2017-2024) huangyh45@mail3.sysu.edu.cn
```

#### revertransseq.pl
Requirements:  
- perl >=5.26.2
- perl-bioperl >=1.7.2
- perl-parallel-forkmanager >=2.02

```
perl scripts/revertransseq.pl  
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
  
Written by Yu-Hao Huang (2017-2024) huangyh45@mail3.sysu.edu.cn
```

#### alignseq.pl
Requirements:  
- perl >=5.26.2
- mafft >=7.467
- transseq.pl
- revertransseq.pl

```
perl scripts/alignseq.pl  
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
  
Written by Yu-Hao Huang (2017-2024) huangyh45@mail3.sysu.edu.cn
```

#### connect.pl
Requirements:  
- perl >=5.26.2
- perl-bioperl >=1.7.2

```
perl scripts/connect.pl  
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
  
Written by Yu-Hao Huang (2018-2024) huangyh45@mail3.sysu.edu.cn  
```

#### merge_seqs.py
Requirements:
- python >=3.7.4

The script can be used to merge the output alignments in different PhyloAln output directories with the same reference alignments, for example, for data of different batches.  
Usage:  
```
scripts/merge_seqs.py output_dir PhyloAln_dir1 PhyloAln_dir2 (PhyloAln_dir3 ...)
```

#### select_seqs.py
Requirements:
- python >=3.7.4

The script can be used to select or exclude a list of species (first space separated by separate_symbol from the sequence name) or sequences (the sequence name) from the sequence FASTA files with the same suffix in a directory and output the managed sequence files to a new diretory.  
Usage:  
```
scripts/select_seqs.py input_dir selected_species_or_sequences(separated by comma) output_dir fasta_suffix(default='.fa') separate_symbol(default='.') if_list_for_exclusion(default=no)
```

#### trim_matrix.py
Requirements:
- python >=3.7.4

The script can be used to trim first the colomns (sites) and/or then the rows (sequences) in the sequence matrixes in FASTA files with the same suffix in a directory based on the unknown sites and output the managed sequence files to a new diretory.  
Usage:  
```
scripts/trim_matrix.py input_dir output_dir unknown_symbol(default='X') known_number(>=1)_or_percent(<1)_for_columns(default=0.5) known_number(>=1)_or_percent(<1)_for_rows(default=0) fasta_suffix(default='.fa')
```

#### root_tree.py
Requirements:
- python >=3.7.4
- ete3 >=3.1.2

The script can be used to root the tree with NEWICK format by ETE 3 package and output the rooted NEWICK tree file.  
Usage:  
```
scripts/root_tree.py input.nwk output.nwk outgroup/outgroups(seperated by comma)
```

#### test_effect.py
Requirements:
- python >=3.7.4

The script can be used to calculate the completeness and percent identity of the alignments in FASTA files with the same suffix in a directory compared with the reference alignments in another directory, mainly for testing the effect of reference-based alignment tools, such as PhyloAln.  
Usage:  
```
scripts/test_effect.py reference_dir:ref_species_or_seq_name target_dir:target_species_or_seq_name output_tsv unknown_symbol(default='N') separate(default='.') fasta_suffix(default='.fa') selected_species_or_sequences(separated by comma)
```

### Citation
Huang Y-H, Sun Y-F, Li H, Li H-S, Pang H. PhyloAln: A Convenient Reference-Based Tool to Align Sequences and High-Throughput Reads for Phylogeny and Evolution in the Omic Era. Molecular Biology and Evolution. 2024;41(7):msae150. https://doi.org/10.1093/molbev/msae150

### Questions & Answers

#### How can I obtain the reference alignments and the final tree?
We do not provide the upstream preparation of the reference alignments and the downstream phylogenetic analyses in PhyloAln. You can mannually collect the reference sequences, align them to generate the reference alignments and build the tree. These steps are flexible as you like. The reference alignments are recommended to contain an outgroup for foreign decontamination in PhyloAln and rooting tree. A detailed practice of phylogenomics using nuclear single-copy protein-coding genes can be seen here ([A practice using PhyloAln for phylogenomics](#a-practice-using-phyloaln-for-phylogenomics)). For other types of the data, such as non-protein-coding genes or genes with non-standard genetic codes, you can collect the reference sequences from [NCBI](https://www.ncbi.nlm.nih.gov/) or other places, and additionally adjust the options of alignseq.pl and PhyloAln.
#### Does selection of the outgroup influence detection of foreign contamination? How can I choose an appropriate outgroup?
Actually, in a specific reference alignment, selection of the outgroup have minimum impact on the results through our test (see [our article](https://doi.org/10.1093/molbev/msae150) for detail). Therefore, if you are not sure which species should be the outgroup, you can tentatively not defined the outgroup, and PhyloAln will acquiescently use the first sequences in the reference alignments as the outgroup.  
But when preparing the reference alignments, it should be noticed that the evolutionary distance between the ingroups and the defined outgroup may have impact on detection of foreign contamination based on conservative score. The contamination from the species phylogenetically close to the reference species is relatively hard to be distinguished from the clean ingroup sequences, compared with the contamination from the species distinct from all the reference species, such as symbiotic bacteria of the target eukaryotic species. If the defined outgroup species is too divergent from the ingroups, a large amount of foreign contamination, especially those from species closer to the ingroups than the defined outgroup species, may not be detected and removed.   
Consequently, it should be better that the users have priori knowledge of choosing the defined outgroup when constructing or obtaining the reference alignments. In most cases, the defined outgroup in PhyloAln is recommended to be from close or sister group of the monophyletic ingroup. If several outgroup species are used for phylogenetic reconstruction, an outgroup with most divergent position with the ingroups in the known phylogeny could be selected in PhyloAln. In addition, the sensitivity of detection can be manually adjusted by setting a weight coefficient, which is default as 0.9 (see `--outgroup_weight` in [parameters](#detailed-parameters) for detail). 
#### The required memory is too large to run PhyloAln.
By default, the step to prepare the sequences/reads is in parallel and thus memory-consuming, especially when the data is large. You can try adding the option `--low_mem` to use a low-memory but slower mode to prepare the sequences/reads. In addition, decompression of the ".gz"-ended files will spend some memory. You can also try decompressing the files manually and then running PhyloAln.
#### The positions of sites in the reference alignments are changed in the output alignments.
When HMMER3 search, some non-conservative sites are deleted (e.g., gappy sites) or sometimes realigned. This has little impact on the downstream phylogenetic or evolutionary analyses. If you want to remain unchanged reference alignments or need special HMMER3 search, you can try utilizing the options `--hmmbuild_parameters` and `--hmmbuild_parameters` to control the parameters of HMMER3. For example, you can try adding the option `--hmmbuild_parameters '--symfrac' '0'` to remain the gappy sites.
#### How can I assemble the paired-end reads?
PhyloAln does not have the method to specifically assemble the paired-end reads. It only mapped all the sequences/reads into the alignments and build a consensus in the assemble and/or output steps. You can input both two paired-end read files for a single source/species (see `-i` and `-c` in [parameters](#detailed-parameters) for detail). Furthermore, if you focus on the effect of assembly using paired-end reads, you can first merge them by other tools (e.g., [fastp](https://github.com/OpenGene/fastp)) and then input the merged read files with or without unpaired read files into PhyloAln.
#### Can PhyloAln generate the alignments of multiple-copy genes for gene family analyses?
Actually, we have designed options of this possibility. You can try it like this:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m codon -u outgroup -z all --overlap_len overlap_len --overlap_pident overlap_pident
```
`-z all` represents outputing all the assembled sequences instead of consensus of them. And you can adjust `--overlap_len` and `--overlap_pident` to find a best assembly for the genes.
