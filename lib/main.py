#!/usr/bin/env python3

import sys
import os
import argparse
try:
	import library as lib
	import map
	import assemble as ab
except ImportError:
	import lib.library as lib
	import lib.map as map
	import lib.assemble as ab

def main(args):
	# the arguments
	parser = argparse.ArgumentParser(prog='PhyloAln', usage="%(prog)s [options] -a reference_alignment_file -s species -i fasta_file -f fasta -o output_directory\n%(prog)s [options] -d reference_alignments_directory -c config.tsv -f fastq -o output_directory", description='''A program to directly generate alignments from FASTA/FASTQ files based on reference alignments for phylogenetic analyses.''', epilog="""Written by Yu-Hao Huang (2023) huangyh45@mail2.sysu.edu.cn""")
	parser.add_argument('-a', '--aln', type=os.path.abspath, help='the single reference FASTA alignment file')
	parser.add_argument('-d', '--aln_dir', type=os.path.abspath, help='the directory containing all the reference FASTA alignment files')
	parser.add_argument('-x', '--aln_suffix', default='.fa', help='the suffix of the reference FASTA alignment files when using "-d"(default:%(default)s)')
	parser.add_argument('-s', '--species', help='the studied species ID for the provided FASTA/FASTQ files(-i)')
	parser.add_argument('-i', '--input', type=os.path.abspath, nargs='+', help='the input FASTA/FASTQ file(s) of the single species(-s), compressed files ending with ".gz" are allowed')
	parser.add_argument('-c', '--config', type=os.path.abspath, help="the TSV file with the format of 'species	sequence_file(s)(absolute path, files separated by commas)' per line for multiple species")
	parser.add_argument('-f', '--file_format', choices=['guess', 'fastq', 'fasta', 'large_fasta'], default='guess', help="the file format of the provided FASTA/FASTQ files, 'large_fasta' is recommended for speeding up reading the FASTA files with long sequences(e.g. genome sequences) and cannot be guessed(default:%(default)s)")
	parser.add_argument('-o', '--output', default='PhyloAln_out', type=os.path.abspath, help='the output directory containing the results(default:%(default)s)')
	parser.add_argument('-p', '--cpu', type=int, default=8, help="maximum threads to be totally used in parallel tasks(default:%(default)d)")
	parser.add_argument('--parallel', type=int, help="number of parallel tasks for each alignments, number of CPUs used for single alignment will be automatically calculated by '--cpu / --parallel'(default:the smaller value between number of alignments and the maximum threads to be used)")
	parser.add_argument('-m', '--mol_type', choices=['dna', 'prot', 'codon', 'dna_codon'], default='dna', help="the molecular type of the reference alignments(default:%(default)s, 'dna' suitable for nucleotide-to-nucleotide or protein-to-protein alignment, 'prot' suitable for protein-to-nucleotide alignment, 'codon' and 'dna_codon' suitable for codon-to-nucleotide alignment based on protein and nucleotide alignments respectively)")
	parser.add_argument('-g', '--gencode', type=int, default=1, help="the genetic code used in translation(default:%(default)d = the standard code, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)")
	parser.add_argument('--ref_split_len', type=int, help="If provided, split the reference alignments longer than this length into short alignments with this length, ~1000 may be recommended for concatenated alignments, and codon alignments should be devided by 3")
	parser.add_argument('-l','--split_len', type=int, help="If provided, split the sequences longer than this length into short sequences with this length, 200 may be recommended for long genomic reads or sequences")
	parser.add_argument('--split_slide', type=int, help="the slide to split the sequences using sliding window method(default:half of '--split_len')")
	parser.add_argument('-n', '--no_reverse', action='store_true', help="not to prepare and search the reverse strand of the sequences, recommended for searching protein or CDS sequences")
	parser.add_argument('--low_mem', action='store_true', help="use a low-memory but slower mode to prepare the reads, 'large_fasta' format is not supported and gz compressed files may still spend some memory")
	parser.add_argument('--hmmbuild_parameters', nargs='+', default=[], help="the parameters when using HMMER hmmbuild for reference preparation(default:%(default)s)")
	parser.add_argument('--hmmsearch_parameters', nargs='+', default=[], help="the parameters when using HMMER hmmsearch for mapping the sequences(default:%(default)s)")
	parser.add_argument('-b', '--no_assemble', action='store_true', help="not to assemble the raw sequences based on overlap regions")
	parser.add_argument('--overlap_len', type=int, default=30, help="minimum overlap length when assembling the raw sequences(default:%(default)d)")
	parser.add_argument('--overlap_pident', type=float, default=98, help="minimum overlap percent identity when assembling the raw sequences(default:%(default).2f)")
	parser.add_argument('-t', '--no_out_filter', action='store_true', help="not to filter the foreign or no-signal sequences based on conservative score")
	parser.add_argument('-u', '--outgroup', help="the outgroup species for foreign or no-signal sequences detection(default:the first sequence in each alignment)")
	parser.add_argument('-q', '--sep', default='.', help="the separate symbol between species name and gene identifier in the sequence headers of the alignments(default:%(default)s)")
	parser.add_argument('--outgroup_weight', type=float, default=0.9, help="the weight coefficient to adjust strictness of the foreign or no-signal sequence filter, small number or decimal means ralaxed criterion (default:%(default).2f, 1 = not adjust)")
	parser.add_argument('-r', '--no_cross_species', action='store_true', help="not to remove the cross contamination for multiple species")
	parser.add_argument('--cross_overlap_len', type=int, default=30, help="minimum overlap length when cross contamination detection(default:%(default)d)")
	parser.add_argument('--cross_overlap_pident', type=float, default=98, help="minimum overlap percent identity when cross contamination detection(default:%(default).2f)")
	parser.add_argument('--min_exp', type=float, default=0.2, help="minimum expression value when cross contamination detection(default:%(default).2f)")
	parser.add_argument('--min_exp_fold', type=float, default=2, help="minimum expression fold when cross contamination detection(default:%(default).2f)")
	parser.add_argument('-w', '--unknow_symbol', default='unknow', help="the symbol representing unknown bases for missing regions(default:%(default)s = 'N' in nucleotide alignments and 'X' in protein alignments)")
	parser.add_argument('-z', '--final_seq', choices=['consensus', 'consensus_strict', 'all', 'expression', 'length'], default='consensus', help="the mode to output the sequences(default:%(default)s, 'consensus' means selecting most common bases from all sequences, 'consensus_strict' means only selecting the common bases and remaining the different bases unknow, 'all' means remaining all sequences, 'expression' means the sequence with highest read counts after assembly, 'length' means sequence with longest length")
	parser.add_argument('-y', '--no_ref', action='store_true', help="not to output the reference sequences")
	parser.add_argument('-v', '--version', action='version', version="%(prog)s v0.1")
	args = parser.parse_args(args)

	# parse the alignment files
	alns = {}
	if args.ref_split_len:
		if args.aln_dir:
			print("\nError: split of the alignment is not supported for multiple alignments, please input a single alignment file through '-a' or '--aln'!")
			sys.exit(1)
		elif args.final_seq == 'all':
			print("\nError: split of the alignment is not supported to output all sequences, please choice other options to keep unqiue sequences instead!")
			sys.exit(1)
		elif args.aln:
			alns = map.split_ref(args.aln, args.ref_split_len)
	elif args.aln:
		alns['aln'] = args.aln
	elif args.aln_dir:
		for filename in os.listdir(args.aln_dir):
			if filename.endswith(args.aln_suffix):
				alns[filename.replace(args.aln_suffix, '')] = os.path.join(args.aln_dir, filename)
	if not alns:
		print("\nError: fail to find any alignment!")
		sys.exit(1)

	# parse the species data
	rawdata = {}
	if args.species and args.input:
		rawdata[args.species] = args.input
	elif args.config:
		for line in open(args.config):
			arr = line.rstrip().split("\t")
			rawdata[arr[0]] = arr[1].split(',')
	if not rawdata:
		print("\nError: fail to find any species data!")
		sys.exit(1)

	# check the unknow symbol, low-memory mode and outgroup, and parse the parallel task number and length to split
	if args.unknow_symbol != 'unknow' and len(args.unknow_symbol) > 1:
		print("\nError: the symbol representing unknown bases should be single character!")
		sys.exit(1)
	if args.low_mem and args.file_format == 'large_fasta':
		print("\nError: the format of 'large_fasta' is not supported in the low-memory mode! If you want to use the low-memory mode, you can use 'fasta' format and it will take a while!")
		sys.exit(1)
	if not args.outgroup and not args.no_out_filter:
		print("\nWarning: no outgroup was set, the first sequence in each alignment will be considered as outgroup in foreign sequence filter!")
	if args.parallel is None:
		args.parallel = min(len(alns), args.cpu)
	else:
		args.parallel = min(args.parallel, len(alns), args.cpu)
	if args.split_len is not None and args.split_slide is None:
		args.split_slide = int(args.split_len / 2)

	# create and enter the output directory, and output the splitted reference alignments
	if not os.path.isdir(args.output):
		os.makedirs(args.output)
	os.chdir(args.output)
	if not os.path.isdir('ok'):
		os.mkdir('ok')
	if args.ref_split_len:
		if not os.path.isdir('ref_split'):
			os.mkdir('ref_split')
		for aln, aln_seqstr in alns.items():
			outfile = open(os.path.join('ref_split', aln + '.fa'), 'w')
			outfile.write(aln_seqstr)
			outfile.close()
			alns[aln] = os.path.join('ref_split', aln + '.fa')

	# prepare the reference HMMs
	map.prepare_ref(alns, cpu=args.cpu, np=args.parallel, moltype=args.mol_type, gencode=args.gencode, parameters=args.hmmbuild_parameters)

	# map (by HMMER), extract, assemble the sequences and remove foreign or no-signal sequences of each species
	total_reads = {}
	assemblers = {}
	for sp, fastxs in rawdata.items():
		all_seqs, total_reads[sp] = map.map_reads(alns, sp, fastxs, file_format=args.file_format, cpu=args.cpu, np=args.parallel, moltype=args.mol_type, gencode=args.gencode, split_len=args.split_len, split_slide=args.split_slide, no_reverse=args.no_reverse, low_mem=args.low_mem, parameters=args.hmmsearch_parameters)
		all_hmmres = map.extract_reads(alns, sp, all_seqs, moltype=args.mol_type, split_len=args.split_len, low_mem=args.low_mem)
		del all_seqs
		assemblers[sp] = ab.generate_assembly_mp(alns, sp, all_hmmres, np = args.parallel, moltype=args.mol_type, gencode=args.gencode, no_assemble=args.no_assemble, overlap_len=args.overlap_len, overlap_pident=args.overlap_pident, no_out_filter=args.no_out_filter, outgroup=args.outgroup, sep=args.sep, outgroup_weight=args.outgroup_weight, final_seq=args.final_seq)
	
	# cross decontamination and output
	for group_name in alns.keys():
		assemblers[group_name] = {}
		for sp in rawdata.keys():
			assemblers[group_name][sp] = assemblers[sp][group_name]
	for sp in rawdata.keys():
		assemblers.pop(sp)
	ab.cross_and_output_mp(alns.keys(), list(rawdata.keys()), assemblers, total_reads, np = args.parallel, moltype=args.mol_type, gencode=args.gencode, no_assemble=args.no_assemble, no_cross_species=args.no_cross_species, min_overlap=args.cross_overlap_len, min_pident=args.cross_overlap_pident, min_exp=args.min_exp, min_fold=args.min_exp_fold, unknow=args.unknow_symbol, final_seq=args.final_seq, no_ref=args.no_ref, sep=args.sep)

	# concatenate the output alignments if split the reference alignment
	if args.ref_split_len:
		if args.unknow_symbol == 'unknow':
			fill = 'N'
		else:
			fill = args.unknow_symbol
		ab.concatenate('nt_out', alns, fill)
		if args.mol_type != 'dna':
			if args.unknow_symbol == 'unknow':
				fill = 'X'
			else:
				fill = args.unknow_symbol
			ab.concatenate('aa_out', alns, fill)

if __name__ == "__main__":
    main(sys.argv[1:])
