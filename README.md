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
