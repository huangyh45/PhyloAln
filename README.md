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
`git clone https://github.com/DessimozLab/read2tree.git`

Then, you can test if PhyloAln has been available using the commands:  
`cd /your-PhyloAln-path/`  
`export PATH=$PATH:/your-PhyloAln-path/:/your-PhyloAln-path/scripts`  
`bash tests/run_test.sh`

#### 2) Installation using Conda
`conda install phyloaln`
