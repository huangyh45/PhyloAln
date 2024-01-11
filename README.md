# PhyloAln
## PhyloAln: a convenient reference-based tool to align sequences and high-throughput reads for phylogeny and evolution in the omic era

PhyloAln is a reference-based alignment tool for phylogeny and evolution. PhyloAln can directly map not only the raw reads but also assembled or translated sequences into the reference alignments, which is suitable for different omic data and skips the complex preparation in the traditional method including data assembly, gene prediction, orthology assignment and sequence alignment, with a relatively high accuracy in the aligned sites and the downstream phylogeny. It is also able to detect and remove foreign and cross contamination in the generated alignments, which is not considered in other reference-based methods, and thus improve the quality of the alignments for downstream analyses.

### Installation

#### 1) Installation from source
##### Requirements
- python >=3.7.4
- biopython >=1.77
- hmmer >=3.1
- mafft >=7.467
- ete3 >=3.1.2
- perl >=5.26.2
- perl-bioperl >=1.7.2
- perl-parallel-forkmanager >=2.02

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
PhyloAln -d reference_alignments_directory -c config.tsv -x alignment_file_name_suffix -o output_directory -p 20 -m codon -u outgroup -l 200
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
