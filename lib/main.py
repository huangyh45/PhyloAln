#!/usr/bin/env python3

import sys
import os
import argparse
try:
	import aln as al
except ImportError:
	import lib.aln as al

def parameter(x):
	return str(x).strip()

def main(args):
	# the arguments
	parser = argparse.ArgumentParser(prog='PhyloAln', usage="%(prog)s [options] -a reference_alignment_file -s species -i fasta_file -f fasta -o output_directory\n%(prog)s [options] -d reference_alignments_directory -c config.tsv -f fastq -o output_directory", description="A program to directly generate multiple sequence alignments from FASTA/FASTQ files based on reference alignments for phylogenetic analyses.\nCitation: Huang Y-H, Sun Y-F, Li H, Li H-S, Pang H. 2024. MBE. 41(7):msae150. https://doi.org/10.1093/molbev/msae150", epilog="Written by Yu-Hao Huang (2023-2025) huangyh45@mail3.sysu.edu.cn", formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-a', '--aln', type=os.path.abspath, help='the single reference FASTA alignment file, or STO(for search mode: hmmer-sto, mmseqs-sto or infernal-sto)/HMM(for search mode: hmmer-hmm or hmmer-db)/CM(for search mode: infernal-cm or infernal-db) file')
	parser.add_argument('-d', '--aln_dir', type=os.path.abspath, help='the directory containing all the reference FASTA alignment files or STO/HMM/CM files')
	parser.add_argument('-x', '--aln_suffix', default='.fa', help='the suffix of the reference FASTA alignment files when using "-d"(default:%(default)s)')
	parser.add_argument('-s', '--species', help='the studied species ID for the provided FASTA/FASTQ files(-i)')
	parser.add_argument('-i', '--input', type=os.path.abspath, nargs='+', help='the input FASTA/FASTQ file(s) of the single species(-s), compressed files ending with ".gz" are allowed')
	parser.add_argument('-c', '--config', type=os.path.abspath, help="the TSV file with the format of 'species	sequence_file(s)(absolute path, files separated by commas)' per line for multiple species")
	parser.add_argument('-f', '--file_format', choices=['guess', 'fastq', 'fasta'], default='guess', help="the file format of the provided FASTA/FASTQ files(default:%(default)s)")
	parser.add_argument('-o', '--output', default='PhyloAln_out', type=os.path.abspath, help='the output directory containing the results(default:%(default)s)')
	parser.add_argument('-p', '--cpu', type=int, default=8, help="maximum threads to be totally used in parallel tasks(default:%(default)d)")
	parser.add_argument('--parallel', type=int, help="number of parallel tasks for each alignments, number of CPUs used for single alignment will be automatically calculated by '--cpu / --parallel'(default:the smaller value between number of alignments and the maximum threads to be used)")
	parser.add_argument('-e', '--mode', choices=['dna2reads', 'prot2reads', 'codon2reads', 'fast_dna2reads', 'fast_prot2reads', 'fast_codon2reads', 'dna2trans', 'prot2trans', 'codon2trans', 'dna2genome', 'prot2genome', 'codon2genome', 'rna2rna', 'prot2prot', 'codon2codon', 'gene_dna2dna', 'gene_rna2rna', 'gene_codon2codon', 'gene_codon2dna', 'gene_prot2prot', 'gene_prot2dna', 'gene_dna2genome', 'gene_codon2genome', 'gene_prot2genome'], help="the common mode to automatically set the parameters for easy use(**NOTICE: if you manually set those parameters, the parameters you set will be ignored and covered! See https://github.com/huangyh45/PhyloAln/blob/main/README.md#example-commands-for-different-data-and-common-mode-for-easy-use for detailed parameters)")
	parser.add_argument('-m', '--mol_type', choices=['dna', 'prot', 'codon', 'dna_codon'], default='dna', help="the molecular type of the reference alignments(default:%(default)s, 'dna' suitable for nucleotide-to-nucleotide or protein-to-protein alignment, 'prot' suitable for protein-to-nucleotide alignment, 'codon' and 'dna_codon' suitable for codon-to-nucleotide alignment based on protein and nucleotide alignments respectively)")
	parser.add_argument('-g', '--gencode', type=int, default=1, help="the genetic code used in translation(default:%(default)d = the standard code, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)")
	parser.add_argument('--ref_split_len', type=int, help="If provided, split the reference alignments longer than this length into short alignments with this length, ~1000 may be recommended for concatenated alignments, and codon alignments should be able to be devided by 3")
	parser.add_argument('-l','--split_len', type=int, help="If provided, split the sequences longer than this length into short sequences with this length, 200 may be recommended for long genomic reads or sequences")
	parser.add_argument('--split_slide', type=int, help="the slide to split the sequences using sliding window method(default:half of '--split_len')")
	parser.add_argument('-j', '--search_mode', choices=['mmseqs', 'mmseqs-sto', 'hmmer', 'hmmer-sto', 'hmmer-hmm', 'hmmer-db', 'infernal-sto', 'infernal-cm', 'infernal-db'], default='hmmer', help="the mode to set the search tool and reference alignment format(default:%(default)s, the Infernal search modes are experimental, sometimes failing to read the search results through biopython)")
	parser.add_argument('-n', '--no_reverse', action='store_true', help="not to prepare and search the reverse strand of the sequences, recommended for searching protein or CDS sequences")
	parser.add_argument('--codon', action='store_true', help="to prepare and search only the first position of the sequences, recommended for searching codon sequences")
	parser.add_argument('--merge_len', type=int, default=250, help="merge the redundant reads/sequences with the length not larger than this value to speed up search and mapping, which will consume some memory(default:%(default)d, set it as 0 to disable this step and speed up file reading, especially when no redundant sequences are believed, such as the assembled sequences or the long sequences without being splitted)")
	parser.add_argument('--mmseqs_convertmsa_parameters', type=parameter, nargs='+', default=['--identifier-field', '0'], help="the parameters when using MMseqs2 convertmsa for reference preparation, with the format of ' --xxx' of each parameter, in which space is required(default:[' --identifier-field', '0'])")
	parser.add_argument('--mmseqs_msa2profile_parameters', type=parameter, nargs='+', default=['--match-mode', '1'], help="the parameters when using MMseqs2 msa2profile for reference preparation, with the format of ' --xxx' of each parameter, in which space is required(default:[' --match-mode', '1'])")
	parser.add_argument('--mmseqs_esearch_parameters', type=parameter, nargs='+', default=['--min-length', '12', '--strand', '2', '-e', '0.1'], help="the parameters when using MMseqs2 easy-search for mapping the sequences, with the format of ' --xxx' of each parameter, in which space is required(default:[' --min-length', '12', ' --strand', '2', ' -e', '0.1'])")
	parser.add_argument('--hmmbuild_parameters', type=parameter, nargs='+', default=[], help="the parameters when using HMMER hmmbuild for reference preparation, with the format of ' --xxx' of each parameter, in which space is required(default:%(default)s)")
	parser.add_argument('--hmmsearch_parameters', type=parameter, nargs='+', default=[], help="the parameters when using HMMER hmmsearch/hmmscan for mapping the sequences, with the format of ' --xxx' of each parameter, in which space is required(default:%(default)s)")
	parser.add_argument('--hmmemit_parameters', type=parameter, nargs='+', default=['-N', '15', '-a'], help="the parameters when using HMMER hmmemit for reference preparation, with the format of ' --xxx' of each parameter, in which space is required(default:[' -N', '15', ' -a'])")
	parser.add_argument('--cmbuild_parameters', type=parameter, nargs='+', default=['-F'], help="the parameters when using Infernal cmbuild for reference preparation, with the format of ' --xxx' of each parameter, in which space is required(default:[' -F'])")
	parser.add_argument('--cmcalibrate_parameters', type=parameter, nargs='+', default=[], help="the parameters when using Infernal cmcalibrate for mapping the sequences, with the format of ' --xxx' of each parameter, in which space is required(default:%(default)s)")
	parser.add_argument('--cmsearch_parameters', type=parameter, nargs='+', default=[], help="the parameters when using Infernal cmsearch/cmscan for mapping the sequences, with the format of ' --xxx' of each parameter, in which space is required(default:%(default)s)")
	parser.add_argument('--cmemit_parameters', type=parameter, nargs='+', default=['-N', '15', '-a', '--dna'], help="the parameters when using Infernal cmemit for reference preparation, with the format of ' --xxx' of each parameter, in which space is required(default:[' -N', '15', ' -a', ' --dna'])")
	parser.add_argument('--emit_sample_num', type=int, default=10, help="the number of times to run hmmemit/cmemit to generate a consensus sequence for constructing reference alignments(default:%(default)d)")
	parser.add_argument('--infernal_no_replace', action='store_true', help="not to replace Us with Ts when reading the Infernal search results for mapping the RNA sequences with Us")
	parser.add_argument('--trim_pos', type=int, default=5, help="maximum terminal positions to be considered to be trimmed at the start and end based on the reference alignments, this length is the number of amino acids when the '-m' is 'prot' or 'codon'(default:%(default)d)")
	parser.add_argument('--extend_pos', type=int, default=0, help="maximim positions to try to extend at the start and end based on the reference alignments when terminal positions are not trimmed, this length is the number of amino acids when the '-m' is 'prot' or 'codon'(default:%(default)d = not to extend)")
	parser.add_argument('--pos_freq', type=float, default=0.2, help="minimum position frequency based on the reference alignment when trim and extend the terminal positions(default:%(default).2f)")
	parser.add_argument('-b', '--no_assemble', action='store_true', help="not to assemble the raw sequences based on overlap regions")
	parser.add_argument('--overlap_len', type=int, default=30, help="minimum overlap length when assembling the raw sequences(default:%(default)d)")
	parser.add_argument('--overlap_pident', type=float, default=98, help="minimum overlap percent identity when assembling the raw sequences(default:%(default).2f)")
	parser.add_argument('-t', '--no_out_filter', action='store_true', help="not to filter the foreign or no-signal sequences based on conservative score")
	parser.add_argument('-u', '--outgroup', nargs='+', default=[], help="the outgroup species for foreign or no-signal sequences detection(default:all the sequences in the alignments with all sequences as ingroups)")
	parser.add_argument('--ingroup', nargs='+', default=[], help="the ingroup species for score calculation in foreign or no-signal sequences detection(default:all the sequences when all sequences are set as outgroups; all other sequences except the outgroups)")
	parser.add_argument('-q', '--sep', default='.', help="the separate symbol between species name and gene identifier in the sequence headers of the alignments(default:%(default)s)")
	parser.add_argument('--outgroup_weight', type=float, default=0.9, help="the weight coefficient to adjust strictness of the foreign or no-signal sequence filter, small number or decimal means ralaxed criterion (default:%(default).2f, 1 = not adjust)")
	parser.add_argument('--intron_len', type=int, help="If provided, try to connect the target fragments with the distance within this value, 20000 may be recommended for long genomic reads or sequences with introns")
	parser.add_argument('-r', '--no_cross_species', action='store_true', help="not to remove the cross contamination for multiple species")
	parser.add_argument('--cross_overlap_len', type=int, default=30, help="minimum overlap length when cross contamination detection(default:%(default)d)")
	parser.add_argument('--cross_overlap_pident', type=float, default=98, help="minimum overlap percent identity when cross contamination detection(default:%(default).2f)")
	parser.add_argument('--min_exp', type=float, default=0.2, help="minimum expression value when cross contamination detection(default:%(default).2f)")
	parser.add_argument('--min_exp_fold', type=float, default=5, help="minimum expression fold when cross contamination detection(default:%(default).2f)")
	parser.add_argument('-w', '--unknow_symbol', default='N', help="the symbol representing unknown bases (nucleotides or amino acids dependent on the input FASTA/FASTQ files) for missing regions(default:%(default)s)")
	parser.add_argument('--unknow_prot', default='X', help="the symbol representing unknown translated amino acids for missing regions(default:%(default)s)")
	parser.add_argument('-z', '--final_seq', choices=['consensus', 'consensus_strict', 'all', 'expression', 'length'], default='consensus', help="the mode to output the sequences(default:%(default)s, 'consensus' means selecting most common bases from all sequences, 'consensus_strict' means only selecting the common bases and remaining the different bases unknow, 'all' means remaining all sequences, 'expression' means the sequence with highest read counts after assembly, 'length' means sequence with longest length")
	parser.add_argument('-y', '--no_ref', action='store_true', help="not to output the reference sequences")
	parser.add_argument('-k', '--keep_seqid', action='store_true', help="keep original sequence IDs in the output alignments instead of renaming them based on the species ID, not recommended when the output mode is 'consensus'/'consensus_strict' or the assembly step is on")
	parser.add_argument('--info_max_seqs', type=int, default=5000, help="maximum target sequence number of each species to record target sequence information, smaller number means using less memory, it will be ignored when '-k'/'--keep_seqid' is on(default:%(default)d)")
	parser.add_argument('-v', '--version', action='version', version="%(prog)s v1.1.0")
	args = parser.parse_args(args)

	# automatically set the parameters when mode is set for easy use and check the search mode
	if args.mode is not None:
		if args.mode in ['rna2rna', 'prot2prot', 'codon2codon']:
			args.no_reverse = True
			args.no_assemble = True
			args.no_cross_species = True
			args.merge_len = 0
			if args.mode == 'codon2codon':
				args.mol_type = 'codon'
				args.codon = True
			else:
				args.mol_type = 'dna'
				if args.mode == 'prot2prot':
					args.unknow_symbol = 'X'
		elif args.mode.startswith('gene_'):
			args.no_assemble = True
			args.no_cross_species = True
			args.merge_len = 0
			args.final_seq = 'all'
			args.keep_seqid = True
			args.unknow_symbol = '-'
			args.unknow_prot = '-'
			if args.mode.endswith('2dna'):
				args.no_reverse = False
			elif args.mode.endswith('2genome'):
				args.no_reverse = False
				args.split_len = 200
				args.intron_len = 20000
			else:
				args.no_reverse = True
				if args.mode.endswith('2codon'):
					args.codon = True
			if args.mode in ['gene_prot2dna', 'gene_prot2genome']:
				args.mol_type = 'prot'
			elif args.mode.startswith('gene_codon2'):
				args.mol_type = 'codon'
			else:
				args.mol_type = 'dna'
			#if args.mode == 'gene_prot2prot':
			#	args.unknow_symbol = 'X'
		else:
			if 'dna2' in args.mode:
				args.mol_type = 'dna'
			elif 'prot2' in args.mode:
				args.mol_type = 'prot'
			elif 'codon2' in args.mode:
				args.mol_type = 'codon'
			if args.mode.endswith('reads') and not args.mode.startswith('fast_'):
				args.no_assemble = False
			else:
				args.no_assemble = True
			if args.mode.endswith('reads'):
				args.no_cross_species = False
			else:
				args.no_cross_species = True
				args.merge_len = 0
				if args.mode.endswith('2genome'):
					args.split_len = 200
	if args.search_mode.startswith('infernal'):
		if args.no_reverse:
			print("\nWarning: Infernal will search both strands itself even if reverse strand search is disabled! If you want to only search the positive strand, please include this parameter in '--cmsearch_parameters': ' --toponly'!")
		else:
			print("\nWarning: reverse strand search is automatically disabled, because Infernal will search both strands itself!")
			args.no_reverse = True
		if args.split_len:
			print("\nWarning: target sequences are not recommended to be splitted before Infernal search, because Infernal is able to search the whole genome sequences and is always used to search targets without introns!")
	elif args.search_mode.startswith('mmseqs'):
		if args.no_reverse:
			print("\nWarning: MMseqs2 will search both strands itself even if reverse strand search is disabled! If you want to only search the positive strand, please try to include or change these parameters in '--mmseqs_esearch_parameters': ' --strand' '1'!")
		else:
			print("\nWarning: reverse strand search is automatically disabled, because MMseqs2 will search both strands itself!")
			args.no_reverse = True
		if args.split_len:
			print("\nError: target sequences are not currently supported to be splitted before MMseqs2 search, because MMseqs2 is able to search the whole genome sequences!")
			sys.exit(1)
	if args.search_mode in ['hmmer-hmm', 'hmmer-db', 'infernal-cm', 'infernal-db'] and not args.no_ref:
		print("\nWarning: reference sequence output is recommended to be disabled through '-y'/'--no_ref', because HMM/CM models have no specific reference sequences!")
	if args.search_mode not in ['hmmer', 'mmseqs'] and args.mol_type == 'codon':
		print(f"\nError: codon reference is not supported in the search mode '{args.search_mode}', please choice the search mode 'hmmer' or 'mmseqs' for codon reference instead!")
		sys.exit(1)

	# parse the alignment files
	alns = {}
	if args.ref_split_len:
		if args.aln_dir:
			print("\nError: split of the alignment is not supported for multiple alignments, please input a single alignment file through '-a' or '--aln'!")
			sys.exit(1)
		elif args.final_seq == 'all':
			print("\nError: split of the alignment is not supported to output all sequences, please choice other options to keep unqiue sequences instead!")
			sys.exit(1)
		elif args.search_mode not in ['hmmer', 'mmseqs']:
			print("\nError: split of the alignment is not supported for STO/HMM/CM files, please choice other search modes instead!")
			sys.exit(1)
		elif args.aln:
			alns = al.split_ref(args.aln, args.ref_split_len)
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

	# check the unknow symbol, codon search, outgroup, and intron length to predict genes, as well as parse the parallel task number and length to split
	if len(args.unknow_symbol) > 1 or len(args.unknow_prot) > 1:
		print("\nError: the symbol representing unknown bases should be single character!")
		sys.exit(1)
	if args.codon and not args.no_reverse:
		print("\nWarning: reverse strand search is automatically disabled in codon search!")
		args.no_reverse = True
	if not args.outgroup and not args.no_out_filter:
		print("\nWarning: no outgroup was set, all the sequences in the alignments will be considered as outgroups with all sequences as ingroups in foreign sequence filter!")
	if args.intron_len and not (args.no_assemble and not args.final_seq.startswith('consensus')):
		print("\nWarning: prediction considering introns is not supported for sequences or reads required to be assembled or used to build a consensus, because it will cause difficult conditions and error, please disable prediction with introns (not set '--intron_len') or assembly step (set '-b'/'--no_assemble') and consensus step (set '-z'/'--final_seq' as 'all'/'expression'/'length')!")
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

	# prepare the reference and target files
	al.prepare_ref(alns, args)
	if args.search_mode.startswith('mmseqs'):
		if not args.mol_type.startswith('dna') and args.gencode != 1:
			args.mmseqs_esearch_parameters.extend(['--translation-table', str(args.gencode)])
		if args.codon:
			if '--strand' in args.mmseqs_esearch_parameters:
				to_remove = args.mmseqs_esearch_parameters.index('--strand')
				args.mmseqs_esearch_parameters.pop(to_remove + 1)
				args.mmseqs_esearch_parameters.pop(to_remove)
			args.mmseqs_esearch_parameters.extend(['--forward-frames', '1', '--strand', '1'])
	else:
		al.prepare_target(rawdata, args)

	if args.search_mode.endswith('-db'):
		# map (by HMMER hmmscan / Infernal cmscan) and extract the target reads/sequences when searching HMMER or Infernal database
		al.map_reads(args)
		alns = al.extract_reads('map_reads.txt', 'extract_reads', args.search_mode, args.infernal_no_replace)
	elif args.search_mode.startswith('mmseqs'):
		# map by MMseqs2 and extract the target reads/sequences
		al.mmseqs2(alns, rawdata, args)
	elif not os.path.isdir('map'):
		os.mkdir('map')

	# map and extract the target reads/sequences (when not to search HMMER or Infernal database), assemble the sequences and remove the contamination of each alignments
	al.assemble_mp(alns, rawdata, args)

	# concatenate the output alignments if split the reference alignment
	if args.ref_split_len:
		al.concatenate('nt_out', alns, args.unknow_symbol)
		if args.mol_type != 'dna':
			al.concatenate('aa_out', alns, args.unknow_prot)

if __name__ == "__main__":
	main(sys.argv[1:])

