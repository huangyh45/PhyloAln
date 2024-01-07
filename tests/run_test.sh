dir=`pwd`
for file in tests/ref/*.fa; do
	alignseq.pl -i $file -o $file.codon_aln.fa -a codon -n 5
done
echo "DMELA	$dir/tests/DMELA_1.fq.gz,$dir/tests/DMELA_2.fq.gz" > tests/run_test.config
echo "DWILL	$dir/tests/DWILL_1.fq.gz,$dir/tests/DWILL_2.fq.gz" >> tests/run_test.config
PhyloAln -d tests/ref -x .fa.codon_aln.fa -c tests/run_test.config -o tests/PhyloAln_out -p 5 -m codon -u SLEBA
connect.pl -i tests/PhyloAln_out/nt_out -f N -b all.block -n -c 123
merge_seqs.py -h
root_tree.py -h
select_seqs.py -h
test_effect.py -h
trim_matrix.py -h
