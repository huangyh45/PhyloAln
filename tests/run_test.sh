set -e
dir=`pwd`
if [ ! -d "tests/aln" ]; then
	mkdir tests/aln
fi
for file in tests/ref/*.fa; do
	og=`basename $file`
	alignseq.pl -i $file -o tests/aln/$og -a codon -n 5
done
echo "DMELA	$dir/tests/DMELA_1.fq.gz,$dir/tests/DMELA_2.fq.gz" > tests/run_test.config
echo "DWILL	$dir/tests/DWILL_1.fq.gz,$dir/tests/DWILL_2.fq.gz" >> tests/run_test.config
PhyloAln -d tests/aln -x .fa -c tests/run_test.config -o tests/PhyloAln_out -p 5 -m codon -u SLEBA
connect.pl -i tests/PhyloAln_out/nt_out -f N -b all.block -n -c 123
if [ $1 ]; then
	# full usage experience of all the scripts with real examples through a simple phylogenomic flow
	PhyloAln -d tests/aln -x .fa -s DWILL2 -i $dir/tests/DWILL_1.fq.gz $dir/tests/DWILL_2.fq.gz -o tests/PhyloAln_out2 -p 5 -m codon -u SLEBA   # additionally generate result alignments for only DWILL
	test_effect.py tests/PhyloAln_out2/nt_out:DYAKU tests/PhyloAln_out/nt_out:DYAKU tests/PhyloAln_DYAKUvsDYAKU.tsv N . .fa DBUSC,DMOJA,DPSEU,DYAKU,SLEBA   # compare result aligned sequences of DYAKU in two runs. The sequences should be 100% identical.
	merge_seqs.py tests/PhyloAln_all tests/PhyloAln_out tests/PhyloAln_out2   # merge results of two runs
	select_seqs.py tests/PhyloAln_all/nt_out DBUSC,DPSEU tests/PhyloAln_all/nt_sel .fa . 1   # exclude the sequences of DBUSC and DPSEU from the result alignments
	trim_matrix.py tests/PhyloAln_all/nt_sel tests/PhyloAln_all/nt_trim N 0.5 10   # exclude the columns (sites) with known bases (not 'N') < 50% and the rows (species) with known bases (not 'N') < 10
	connect.pl -i tests/PhyloAln_all/nt_trim -f N -b all.block -n -c 123   # concatenate the alignments of five genes into a supermatrix and generate a partition file of codon positions for IQ-TREE
	iqtree -s all.fas -p all.block -m MFP+MERGE -B 1000 -T AUTO --threads-max 5 --prefix tests/PhyloAln_all/species_tree   # reconstruct phylogeny of the supermatrix by IQ-TREE
	root_tree.py tests/PhyloAln_all/species_tree.treefile tests/PhyloAln_all/species_tree.rooted.tre SLEBA   # root the species phylogeny using SLEBA (outgroup) and generate the rooted tree 'tests/PhyloAln_all/species_tree.rooted.tre'
	echo "Successfully complete running"
else
	connect.pl -i tests/PhyloAln_out/nt_out -f N -b all.block -n -c 123
	merge_seqs.py -h
	root_tree.py -h
	select_seqs.py -h
	test_effect.py -h
	trim_matrix.py -h
	echo "Successfully installed"
fi
