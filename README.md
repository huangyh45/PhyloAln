# PhyloAln
## PhyloAln: a convenient reference-based tool to align sequences and high-throughput reads for phylogeny and evolution in the omic era

PhyloAln is a reference-based alignment tool for phylogeny and evolution. PhyloAln can directly map not only the raw reads but also assembled or translated sequences into the reference alignments, which is suitable for different omic data and skips the complex preparation in the traditional method including data assembly, gene prediction, orthology assignment and sequence alignment, with a relatively high accuracy in the aligned sites and the downstream phylogeny. It is also able to detect and remove foreign and cross contamination in the generated alignments, which is not considered in other reference-based methods, and thus improve the quality of the alignments for downstream analyses.

### Catalogue
- [Installation](#installation)
  - [1) Installation from source](#1-installation-from-source)
  - [2) Installation using Conda](#2-installation-using-conda)
- [Usage](#usage)
  - [Quick start](#quick-start)
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
- [Questions & Answers](#questions--answers)
  - [The required memory is too large to run PhyloAln.](#the-required-memory-is-too-large-to-run-PhyloAln)
  - [The positions of sites in the reference alignments are changed in the output alignments.](#the-positions-of-sites-in-the-reference-alignments-are-changed-in-the-output-alignments)
  - [Can PhyloAln generate the alignments of multiple-copy genes for gene family analyses?](#can-phyloaln-generate-the-alignments-of-multiple-copy-genes-for-gene-family-analyses)
  
### Installation

#### 1) Installation from source
##### Requirements
- python >=3.7.4
- biopython >=1.77
- hmmer >=3.1
- mafft >=7.467 (optional for the auxiliary scripts)
- ete3 >=3.1.2 (optional for the auxiliary scripts)
- perl >=5.26.2 (optional for the auxiliary scripts)
- perl-bioperl >=1.7.2 (optional for the auxiliary scripts)
- perl-parallel-forkmanager >=2.02 (optional for the auxiliary scripts)

After installing these requirements, you can download PhyloAln from this GitHub repo using:  
```
git clone https://github.com/huangyh45/PhyloAln
```

Then, you can test if PhyloAln has been available using the commands:   
```
cd /your/PhyloAln/path/  
export PATH=$PATH:/your/PhyloAln/path/:/your/PhyloAln/path/scripts  
bash tests/run_test.sh && echo "Successfully installed"
```

#### 2) Installation using Conda
```
conda install phyloaln
```

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

#### Input
PhyloAln needs two types of file:  
- the alignment file(s) with FASTA format. Trimmed alignments with conservative sites are recommended. Multiple alignment files with the same suffix should be placed into a directory for input.
- the sequence/read file(s) with FASTA or FASTQ format. Compressed files ending with ".gz" are allowed. Sequence/read files from multiple sources/species should be inputed through a configure file as described in quick start.

#### Output
PhyloAln generates new alignment file(s) with FASTA format. Each output alignment is corresponding to each reference alignment file, with the aligned target sequences from the provided sequence/read file(s). If using prot, codon or dna_codon mode, the translated protein ailgnments will be also generated. These alignments are mainly for phylogenetic analyses and evolutionary analyses using conservative sites.

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
Map the long reads with high insertion and deletetion rates into the codon alignments:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m dna_codon -u outgroup
```
Map the genomic sequences/reads with intron regions into the codon alignments:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m codon -u outgroup -l 200 -f large_fasta
```
Map the sequences/reads into the codon alignments using the non-standard genetic code (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for detail), for example, the codon alignments of plastd protein-coding genes:
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
  
Written by Yu-Hao Huang (2023) huangyh45@mail2.sysu.edu.cn
```

#### Limitations
- PhyloAln is only designed for phylogenetic analyses and evolutionary analyses with conservative sites, and thus cannot perform de novo assembly due to non-conservative sites.
- Python framework limits the processing time. Specially, speed is largely affected by the numbers of the reference alignments, the numbers of the target sequences/reads and the length of the reference alignments (main impact on the time of HMMER3 search).
- We mainly focus on the flexibility of PhyloAln and thus did not provide the upstream steps of collecting the reference sequences and generating the reference alignments, and the downstream phylogenetic analyses. But you can use the auxiliary scripts to help preparation and perform downstream analyses.

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
  
Written by Yu-Hao Huang (2017-2023) huangyh45@mail2.sysu.edu.cn
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
  
Written by Yu-Hao Huang (2017-2023) huangyh45@mail2.sysu.edu.cn
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
  
Written by Yu-Hao Huang (2017-2023) huangyh45@mail2.sysu.edu.cn
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
  
Written by Yu-Hao Huang (2018-2023) huangyh45@mail2.sysu.edu.cn  
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

### Questions & Answers

#### The required memory is too large to run PhyloAln.
By default, the step to prepare the sequences/reads is in parallel and thus memory-consuming, especailly when the data is large. You can try adding the option `--low_mem` to use a low-memory but slower mode to prepare the sequences/reads. In addition, decompression of the ".gz"-ended files will spend some memory. You can also try decompressing the files manually and then running PhyloAln.
#### The positions of sites in the reference alignments are changed in the output alignments.
When HMMER3 search, some non-conservative sites are deleted (e.g., gappy sites) or sometimes realigned. This has little impact on the downstream phylogenetic or evolutionary analyses. If you want to remain unchanged reference alignments or need special HMMER3 search, you can try utilizing the options `--hmmbuild_parameters` and `--hmmbuild_parameters` to control the parameters of HMMER3. For example, you can try adding the option `--hmmbuild_parameters '--symfrac' '0'` to remain the gappy sites.
#### Can PhyloAln generate the alignments of multiple-copy genes for gene family analyses?
Actually, we have designed options of this possibility. You can try it like this:
```
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m codon -u outgroup -z all --overlap_len overlap_len --overlap_pident overlap_pident
```
`-z all` represents outputing all the assembled sequences instead of consensus of them. And you can adjust `--overlap_len` and `--overlap_pident` to find a best output for the genes.
