#!/usr/bin/env python3

import sys
import os
import shutil
import warnings
from Bio import SearchIO, AlignIO, BiopythonWarning
from Bio.Seq import reverse_complement, translate
try:
	import library as lib
except ImportError:
	import lib.library as lib

# split the reference alignments into short alignments
def split_ref(alnfile, split_len):
	print("\nSplitting the reference alignment into short alignments...")
	seqs = lib.read_fastx(alnfile, 'fasta')
	aln_len = len(list(seqs.values())[0])
	alns = {}
	i = 0
	while i * split_len < aln_len:
		alns['aln_' + str(i)] = ''
		end = i * split_len + split_len
		# for final alignments and avoid too short length < 30
		if end + 30 > aln_len:
			end = aln_len
		for seqid, seqstr in seqs.items():
			alns['aln_' + str(i)] += ">{}\n{}\n".format(seqid, seqstr[i*split_len:end])
		if end >= aln_len:
			break
		i += 1
	return alns

# translate the (gappy) sequences in a FASTA file
def trans_seq(filename, output, gencode=1, dna_codon_unknow=None):
	seqs = lib.read_fastx(filename, 'fasta')
	outfile = open(output, 'w')
	for seqid, seqstr in seqs.items():
		tran_str = ''
		i = 0
		while i < len(seqstr):
			if seqstr[i:i+3] == '---':
				tran_str += '-'
			elif dna_codon_unknow:
				try:
					tran_str += translate(seqstr[i:i+3], table=gencode)
				except:
					tran_str += dna_codon_unknow
			else:
				tran_str += translate(seqstr[i:i+3], table=gencode)
			i += 3
		outfile.write(">{}\n{}\n".format(seqid, tran_str))
	outfile.close()

# construct HMM/CM file by hmmbuild/cmbuild
def hmmbuild(alnfile, group_name, cpu=1, tool='hmmbuild', suffix='.hmm', moltype='dna', gencode=1, parameters=[], calibrate_parameters=[]):
	if tool == 'mmseqs':
		if moltype == 'codon':
			trans_seq(alnfile, os.path.join('ref_hmm', group_name + '.aa.fas'), gencode=gencode)
			return 0, group_name, os.path.join('ref_hmm', group_name + '.aa.fas')
		else:
			return 0, group_name, alnfile
	log = open(os.path.join('ref_hmm', group_name + '.log'), 'w')
	if moltype == 'codon':
		trans_seq(alnfile, os.path.join('ref_hmm', group_name + '.aa.fas'), gencode=gencode)
		cmd = [tool, '-O', os.path.join('ref_hmm', group_name + '.sto'), '-n', group_name]
		if tool == 'hmmbuild':
			cmd.extend(['--cpu', str(cpu)])
		cmd.extend(parameters)
		cmd.extend([os.path.join('ref_hmm', group_name + suffix), os.path.join('ref_hmm', group_name + '.aa.fas')])
	else:
		cmd = [tool, '-O', os.path.join('ref_hmm', group_name + '.sto'), '-n', group_name]
		if tool == 'hmmbuild':
			cmd.extend(['--cpu', str(cpu)])
		cmd.extend(parameters)
		cmd.extend([os.path.join('ref_hmm', group_name + suffix), alnfile])
	ifcomplish = lib.runcmd(cmd, log, stdout=False)
	if ifcomplish:
		if tool == 'cmbuild':
			cmd = ['cmcalibrate', '--cpu', str(cpu)]
			cmd.extend(calibrate_parameters)
			cmd.append(os.path.join('ref_hmm', group_name + '.cm'))
			ifcomplish = lib.runcmd(cmd, log, stdout=False)
		log.close()
		if ifcomplish:
			return 0, group_name
		else:
			return 1, group_name
	else:
		log.close()
		return 1, group_name

# emit alignment from HMM/CM file by hmmemit/cmemit
def hmmemit(alnfile, group_name, cpu=1, tool='hmmemit', suffix='.hmm', parameters=[], sample_num=10, fetch=None):
	log = open(os.path.join('ref_hmm', group_name + '.log'), 'w')
	if fetch:
		cmd = [fetch, '-o', os.path.join('ref_hmm', group_name + suffix), alnfile, group_name]
		ifcomplish = lib.runcmd(cmd, log, stdout=False)
		if not ifcomplish:
			log.close()
			return 1, group_name
	else:
		outfile = open(os.path.join('ref_hmm', group_name + suffix), 'w')
		for line in open(alnfile):
			outfile.write(line)
		outfile.close()
	seqs = {}
	for i in range(sample_num):
		cmd = [tool, '-o', os.path.join('ref_hmm', group_name + '.sto')]
		cmd.extend(parameters)
		cmd.append(os.path.join('ref_hmm', group_name + suffix))
		ifcomplish = lib.runcmd(cmd, log, stdout=False)
		if not ifcomplish:
			log.close()
			return 1, group_name
		seqs['sample' + str(i+1)] = ''
		alignment = AlignIO.read(os.path.join('ref_hmm', group_name + '.sto'), "stockholm")
		cons_pos = alignment._per_col_annotations['reference_annotation']
		for j in range(alignment.get_alignment_length()):
			if cons_pos[j] != '.':
				seq_comp = {}
				for seq in alignment:
					base = str(seq.seq)[j].upper().replace('.', '-').replace('~', '-')
					seq_comp.setdefault(base, 0)
					seq_comp[base] += 1
				seqs['sample' + str(i+1)] += max(seq_comp, key=seq_comp.get)
	outfile = open(os.path.join('ref_hmm', group_name + '.ref.fas'), 'w')
	for seqid, seqstr in seqs.items():
		outfile.write(f">{seqid}\n{seqstr}\n")
	outfile.close()
	os.remove(os.path.join('ref_hmm', group_name + '.sto'))
	log.close()
	return 0, group_name

# prepare the reference files
def prepare_ref(alns, args):
	np = min(args.parallel, len(alns), args.cpu)
	ncpu = int(args.cpu / np)

	print("\nPreparing the reference alignments...")
	if os.path.isfile(os.path.join('ok', 'prepare_ref.ok')):
		print("\nUsing the existing reference files in directory 'ref_hmm'")
	else:
		if not os.path.isdir('ref_hmm'):
			os.mkdir('ref_hmm')
		if not args.search_mode.endswith('-db'):
			args_list = []
			if args.search_mode == 'hmmer-hmm':
				kwds = {'tool': 'hmmemit', 'suffix': '.hmm', 'parameters': args.hmmemit_parameters, 'sample_num': args.emit_sample_num}
				print("\nEmitting alignments from HMMs...")
			elif args.search_mode == 'infernal-cm':
				kwds = {'tool': 'cmemit', 'suffix': '.cm', 'parameters': args.cmemit_parameters, 'sample_num': args.emit_sample_num}
				print("\nEmitting alignments from CMs...")
			elif args.search_mode.startswith('hmmer'):
				kwds = {'tool': 'hmmbuild', 'suffix': '.hmm', 'moltype': args.mol_type, 'gencode': args.gencode, 'parameters': args.hmmbuild_parameters}
				print("\nBuilding HMMs for mapping...")
			elif args.search_mode.startswith('infernal'):
				kwds = {'tool': 'cmbuild', 'suffix': '.cm', 'moltype': args.mol_type, 'gencode': args.gencode, 'parameters': args.cmbuild_parameters, 'calibrate_parameters': args.cmcalibrate_parameters}
				print("\nBuilding CMs for mapping...")
				lib.check_programs(['cmcalibrate'])
			elif args.search_mode.startswith('mmseqs'):
				kwds = {'tool': 'mmseqs', 'moltype': args.mol_type, 'gencode': args.gencode}
				print("\nBuilding profile for mapping...")
			lib.check_programs([kwds['tool']])
			usedcpu = ncpu * np
			for group_name, alnfile in alns.items():
				if usedcpu < args.cpu:
					args_list.append((alnfile, group_name, ncpu+1))
					usedcpu += 1
				else:
					args_list.append((alnfile, group_name, ncpu))
			if kwds['tool'].endswith('build') or kwds['tool'] == 'mmseqs':
				iferrors = lib.run_mp(hmmbuild, args_list, np, kwds=kwds)
			else:
				iferrors = lib.run_mp(hmmemit, args_list, np, kwds=kwds)
			if kwds['tool'] == 'mmseqs':
				if args.search_mode == 'mmseqs':
					aln_format = 'fasta'
				else:
					aln_format = 'stockholm'
				outfile = open('first2alnid.tsv', 'w')
				alignments = []
				for iferror in iferrors:
					alignment = AlignIO.read(iferror[2], aln_format)
					if args.mol_type != 'codon':
						AlignIO.write([alignment], os.path.join('ref_hmm', iferror[1] + '.ref.fas'), "fasta")
					alignments.append(alignment)
					for record in alignment:
						outfile.write(f"{record.id}\t{iferror[1]}\n")
						break
				outfile.close()
				AlignIO.write(alignments, os.path.join('ref_hmm', 'ref.sto'), "stockholm")
				del alignments
				log = open(os.path.join('ref_hmm', 'mmseqs_ref.log'), 'w')
				cmd = ['mmseqs', 'convertmsa', os.path.join('ref_hmm', 'ref.sto'), os.path.join('ref_hmm', 'ref_msadb')]
				cmd.extend(args.mmseqs_convertmsa_parameters)
				lib.runcmd(cmd, log, error=True)
				cmd = ['mmseqs', 'msa2profile', os.path.join('ref_hmm', 'ref_msadb'), os.path.join('ref_hmm', 'ref_profile'), '--threads', str(args.cpu)]
				cmd.extend(args.mmseqs_msa2profile_parameters)
				lib.runcmd(cmd, log, error=True)
				log.close()
			else:
				errors = []
				for iferror in iferrors:
					if iferror[0] == 1:
						errors.append(iferror[1])
				if errors:
					print("\nError in {} commands: {}".format(kwds['tool'], ', '.join(errors)))
					sys.exit(1)
		else:
			if args.search_mode == 'hmmer-db':
				hmmsuffix = '.hmm'
			else:
				hmmsuffix = '.cm'
			outfile = open(os.path.join('ref_hmm', 'ref' + hmmsuffix), 'w')
			for group_name, alnfile in alns.items():
				for line in open(alnfile):
					outfile.write(line)
			outfile.close()
			if args.search_mode == 'hmmer-db':
				lib.check_programs(['hmmpress'])
				log = open(os.path.join('ref_hmm', 'hmmpress.log'), 'w')
				cmd = ['hmmpress', '-f', os.path.join('ref_hmm', 'ref.hmm')]
			else:
				lib.check_programs(['cmpress'])
				log = open(os.path.join('ref_hmm', 'cmpress.log'), 'w')
				cmd = ['cmpress', '-F', os.path.join('ref_hmm', 'ref.cm')]
			lib.runcmd(cmd, log, error=True)
			log.close()
		open(os.path.join('ok', 'prepare_ref.ok'), 'w').close()

# convert a raw FASTA file to a reversed and/or translated FASTA file
def raw2target(raw, target, moltype='dna', gencode=1, no_reverse=False, codon=False):
	outfile = open(target, 'w')
	# if the file has size > 20Mb, use Bio.SeqIO index_db to read it
	if os.path.getsize(raw) > 20000000:
		seqs = lib.read_fastx(raw, 'fasta_db')
	else:
		seqs = lib.read_fastx(raw, 'fasta')
	for seqid, seqstr0 in seqs.items():
		if moltype.startswith('dna'):
			outfile.write(">{}\n{}\n".format(seqid, seqstr0))
			if not no_reverse:
				outfile.write(">{}_rev\n{}\n".format(seqid, reverse_complement(seqstr0)))
		else:
			# supress the warnings due to the last incomplete codons when translate the sequences
			with warnings.catch_warnings():
				warnings.simplefilter('ignore', BiopythonWarning)
				for j in [1,2,3]:
					seqstr = seqstr0[(j-1):]
					outfile.write(">{}_pos{}\n{}\n".format(seqid, j, translate(seqstr, table=gencode)))
					if codon:
						break
				if not no_reverse:
					seqstr0 = reverse_complement(seqstr0)
					for j in [1,2,3]:
						seqstr = seqstr0[(j-1):]
						outfile.write(">{}_rev_pos{}\n{}\n".format(seqid, j, translate(seqstr, table=gencode)))
	outfile.close()

# if there are more than 1000 raw target files, merge them into 100 files and enlarge size of each file
def check_file_num(i, each_size):
	if i < 1000:
		return i, each_size
	for i in range(100, 1000):
		outfile = open('all.raw.fasta.' + str(i % 100), 'a')
		for line in open('all.raw.fasta.' + str(i)):
			outfile.write(line)
		outfile.close()
		os.remove('all.raw.fasta.' + str(i))
	return 100, 10 * each_size

# prepare the target FASTA file: convert the raw FASTQ/FASTA files to (translated) FASTA file
def prepare_target(rawdata, args):
	raw_fasta = 'all.raw.fasta'
	target_fasta = 'all.target.fasta'
	if os.path.isfile(os.path.join('ok', "prepare_target.ok")):
		print("\nUsing the existing target FASTA file")
	else:
		print("\nPreparing the target FASTA file...")
		total_counts = {}
		to_merge = {}
		each_size = 10000000
		data_size = 0
		i = 0
		outfile = open(raw_fasta + '.' + str(i), 'w')
		for sp, fastxs in rawdata.items():
			total_counts[sp] = 0
			for j in range(len(fastxs)):
				fastx_iter, file_format = lib.read_fastx(fastxs[j], args.file_format, return_iter=True)
				seqid = None
				seqstr = None
				line_num = 0
				split_start = 0
				for line in fastx_iter:
					line = line.rstrip()
					if file_format == 'fastq' and line.startswith('@') and line_num % 4 == 0:
						seqid = line.replace('@', '', 1).replace(' ', '_')
					elif file_format == 'fasta' and line.startswith('>'):
						if seqid:
							if seqstr:
								if args.split_len:
									tagid = f"{seqid}_PhAlTag_{sp}_fastx{j+1}_split_{split_start+1}_{split_start+len(seqstr)}"
								else:
									tagid = f"{seqid}_PhAlTag_{sp}_fastx{j+1}"
								# merge the redundant reads/sequences not longer than args.merge_len
								if len(seqstr) <= args.merge_len and to_merge.get(seqstr) is not None:
									to_merge[seqstr].append(tagid)
								else:
									outfile.write(f">{tagid}\n{seqstr}\n")
									data_size += len(seqstr)
									# output the sequences to a new file if this file has data >= 10Mb
									if data_size >= each_size:
										outfile.close()
										data_size = 0
										i += 1
										i, each_size = check_file_num(i, each_size)
										outfile = open(raw_fasta + '.' + str(i), 'w')
									if len(seqstr) <= args.merge_len:
										to_merge[seqstr] = [tagid]
							total_counts[sp] += 1
							split_start = 0
						seqid = line.lstrip('>').split()[0]
						seqstr = ''
						# directly output the sequences when not split and merge the sequences to speed up the file reading
						if args.merge_len == 0 and not args.split_len:
							if data_size >= each_size:
								outfile.close()
								data_size = 0
								i += 1
								i, each_size = check_file_num(i, each_size)
								outfile = open(raw_fasta + '.' + str(i), 'w')
							outfile.write(f">{seqid}_PhAlTag_{sp}_fastx{j+1}\n")
					elif seqid:
						if file_format == 'fastq':
							seqstr = line
							if args.split_len:
								split_start = 0
								while len(seqstr) > args.split_len:
									seqfrag = seqstr[0:args.split_len]
									tagid = f"{seqid}_PhAlTag_{sp}_fastx{j+1}_split_{split_start+1}_{split_start+args.split_len}"
									if len(seqfrag) <= args.merge_len and to_merge.get(seqfrag) is not None:
										to_merge[seqfrag].append(tagid)
									else:
										outfile.write(f">{tagid}\n{seqfrag}\n")
										data_size += args.split_len
										if data_size >= each_size:
											outfile.close()
											data_size = 0
											i += 1
											i, each_size = check_file_num(i, each_size)
											outfile = open(raw_fasta + '.' + str(i), 'w')
										if args.split_len <= args.merge_len:
											to_merge[seqfrag] = [tagid]
									split_start += args.split_slide
									seqstr = seqstr[args.split_slide:]
								tagid = f"{seqid}_PhAlTag_{sp}_fastx{j+1}_split_{split_start+1}_{split_start+len(seqstr)}"
							else:
								tagid = f"{seqid}_PhAlTag_{sp}_fastx{j+1}"
							if len(seqstr) <= args.merge_len and to_merge.get(seqstr) is not None:
								to_merge[seqstr].append(tagid)
							else:
								outfile.write(f">{tagid}\n{seqstr}\n")
								data_size += len(seqstr)
								if data_size >= each_size:
									outfile.close()
									data_size = 0
									i += 1
									i, each_size = check_file_num(i, each_size)
									outfile = open(raw_fasta + '.' + str(i), 'w')
								if len(seqstr) <= args.merge_len:
									to_merge[seqstr] = [tagid]
							total_counts[sp] += 1
							seqid = None
						# directly output the sequences when not split and merge the sequences to spped up the file reading
						elif args.merge_len == 0 and not args.split_len:
							outfile.write(line + "\n")
							data_size += len(line)
						else:
							seqstr += line
							if args.split_len:
								while len(seqstr) > args.split_len:
									seqfrag = seqstr[0:args.split_len]
									tagid = f"{seqid}_PhAlTag_{sp}_fastx{j+1}_split_{split_start+1}_{split_start+args.split_len}"
									if len(seqfrag) <= args.merge_len and to_merge.get(seqfrag) is not None:
										to_merge[seqfrag].append(tagid)
									else:
										outfile.write(f">{tagid}\n{seqfrag}\n")
										data_size += args.split_len
										if data_size >= each_size:
											outfile.close()
											data_size = 0
											i += 1
											i, each_size = check_file_num(i, each_size)
											outfile = open(raw_fasta + '.' + str(i), 'w')
										if args.split_len <= args.merge_len:
											to_merge[seqfrag] = [tagid]
									split_start += args.split_slide
									seqstr = seqstr[args.split_slide:]
					line_num += 1
				if file_format == 'fasta' and seqstr:
					if args.split_len:
						tagid = f"{seqid}_PhAlTag_{sp}_fastx{j+1}_split_{split_start+1}_{split_start+len(seqstr)}"
					else:
						tagid = f"{seqid}_PhAlTag_{sp}_fastx{j+1}"
					if len(seqstr) <= args.merge_len and to_merge.get(seqstr) is not None:
						to_merge[seqstr].append(seqid)
					else:
						outfile.write(f">{tagid}\n{seqstr}\n")
						data_size += len(seqstr)
						if data_size >= each_size:
							outfile.close()
							data_size = 0
							i += 1
							i, each_size = check_file_num(i, each_size)
							outfile = open(raw_fasta + '.' + str(i), 'w')
						if len(seqstr) <= args.merge_len:
							to_merge[seqstr] = [tagid]
					total_counts[sp] += 1
		outfile.close()
		# delete the final output file if it is empty
		if os.path.getsize(raw_fasta + '.' + str(i)) == 0:
			os.remove(raw_fasta + '.' + str(i))
			i = i - 1
		outfile = open('total_counts.tsv', 'w')
		for sp, count in total_counts.items():
			outfile.write(f"{sp}\t{count}\n")
		outfile.close()
		outfile = open('merged_redundance.txt', 'w')
		for seqids in to_merge.values():
			if len(seqids) > 1:
				outfile.write(' '.join(seqids) + "\n")
		outfile.close()
		del total_counts
		del to_merge

		print(f"Packed the sequences from totally {len(rawdata)} species into {i+1} parts and preparing in multiprocess...")
		args_list = []
		kwds = {'moltype': args.mol_type, 'gencode': args.gencode, 'no_reverse': args.no_reverse, 'codon': args.codon}
		for j in range(i + 1):
			args_list.append((raw_fasta + '.' + str(j), target_fasta + '.' + str(j)))
		lib.run_mp(raw2target, args_list, args.cpu, kwds=kwds)
		rawfile = open(raw_fasta, 'w')
		outfile = open(target_fasta, 'w')
		for j in range(i + 1):
			for line in open(raw_fasta + '.' + str(j)):
				rawfile.write(line)
			for line in open(target_fasta + '.' + str(j)):
				outfile.write(line)
			os.remove(raw_fasta + '.' + str(j))
			os.remove(target_fasta + '.' + str(j))
			# delete the index of Bio.SeqIO index_db if it exists
			if os.path.exists(raw_fasta + '.' + str(j) + '.idx'):
				os.remove(raw_fasta + '.' + str(j) + '.idx')
		rawfile.close()
		outfile.close()
		open(os.path.join('ok', 'prepare_target.ok'), 'w').close()

# map the reads to the HMMs/CMs of alignments
def map_reads(args):
	if os.path.isfile(os.path.join('ok', "map_reads.ok")):
		print("\nUsing the existing search result file")
	else:
		print("\nMapping the reads to reference alignments...")
		log = open('map_reads.log', 'w')
		if args.search_mode.startswith('hmmer'):
			lib.check_programs(['hmmscan'])
			cmd = ['hmmscan', '-o', 'map_reads.txt', '--tblout', 'map_reads.tbl', '--cpu', str(args.cpu)]
			cmd.extend(args.hmmsearch_parameters)
			cmd.extend([os.path.join('ref_hmm', 'ref.hmm'), 'all.target.fasta'])
		elif args.search_mode.startswith('infernal'):
			lib.check_programs(['cmscan'])
			cmd = ['cmscan', '-o', 'map_reads.txt', '--tblout', 'map_reads.tbl', '--cpu', str(args.cpu)]
			cmd.extend(args.cmsearch_parameters)
			cmd.extend([os.path.join('ref_hmm', 'ref.cm'), 'all.target.fasta'])
		lib.runcmd(cmd, log, error=True)
		log.close()
		os.remove('all.target.fasta')
		open(os.path.join('ok', "map_reads.ok"), 'w').close()

# extract the target reads from the raw data
def extract_reads(hmmresult, step_name, search_mode='hmmer', infernal_no_replace=False, log=None):
	rep = {}
	if os.path.isfile(os.path.join('ok', step_name + ".ok")):
		if log:
			log.write("\nUsing the existing extracted read files\n")
		else:
			print("\nUsing the existing extracted read files")
		if search_mode.endswith('-db'):
			for file in os.listdir('map'):
				if file.endswith('.hit_info.tsv'):
					rep[file.split('.hit_info.tsv')[0]] = 1
	else:
		if search_mode.endswith('-db'):
			print("\nExtracting the target reads...")
			if os.path.isdir('map'):
				shutil.rmtree('map')
			os.mkdir('map')
		else:
			open(os.path.join('map', step_name.split('extract_reads_')[-1] + '.targets.fas'), 'w').close()
			open(os.path.join('map', step_name.split('extract_reads_')[-1] + '.hit_info.tsv'), 'w').close()
			if log:
				log.write("\nExtracting the target reads...\n")
		last_query = 'nothing'
		outfile = None
		seqs = lib.read_fastx('all.raw.fasta', 'fasta_db', return_iter=True)
		if search_mode.startswith('hmmer'):
			search_format = 'hmmer3-text'
		elif search_mode.startswith('infernal'):
			search_format = 'infernal-text'
		for qresult in SearchIO.parse(hmmresult, search_format):
			for hit in qresult:
				for HSP in hit:
					for HSPfrag in HSP:
						if search_mode.endswith('-db'):
							query_id = HSPfrag.hit_id
							hit_id = HSPfrag.query_id
							query_start = HSPfrag.hit_start
							query_end = HSPfrag.hit_end
							hit_start = HSPfrag.query_start
							hit_end = HSPfrag.query_end
							query_seq = str(HSPfrag.hit.seq).upper().replace('.', '-')
							hit_seq = str(HSPfrag.query.seq).upper().replace('.', '-')
						else:
							query_id = HSPfrag.query_id
							hit_id = HSPfrag.hit_id
							query_start = HSPfrag.query_start
							query_end = HSPfrag.query_end
							hit_start = HSPfrag.hit_start
							hit_end = HSPfrag.hit_end
							query_seq = str(HSPfrag.query.seq).upper().replace('.', '-')
							hit_seq = str(HSPfrag.hit.seq).upper().replace('.', '-')
						if query_id != last_query:
							if outfile:
								outfile.close()
								fasfile.close()
							outfile = open(os.path.join('map', query_id + '.hit_info.tsv'), 'a')
							fasfile = open(os.path.join('map', query_id + '.targets.fas'), 'a')
							rep.setdefault(query_id, {})
							last_query = query_id
						strand = '+'
						start_pos = None
						if search_mode.startswith('hmmer'):
							query_start += 1
							hit_start += 1
							if hit_id.endswith('_pos1') or hit_id.endswith('_pos2') or hit_id.endswith('_pos3'):
								start_pos = int(hit_id[-1]) - 1
								hit_id = hit_id[:-5]
							if hit_id.endswith('_rev'):
								hit_id = hit_id[:-4]
								strand = 'rev'
						elif search_mode.startswith('infernal'):
							if HSPfrag.hit_strand == -1:
								strand = '-'
						if search_mode.startswith('infernal') and not infernal_no_replace:
							query_seq = query_seq.replace('U', 'T')
							hit_seq = hit_seq.replace('U', 'T')
						outfile.write(f"{query_id}\t{hit_id}\t{query_start}\t{query_end}\t{hit_start}\t{hit_end}\t{strand}\t{start_pos}\t{query_seq}\t{hit_seq}\n")
						if rep[query_id].get(hit_id) is None:
							fasfile.write(f">{hit_id}\n{str(seqs[hit_id].seq)}\n")
							rep[query_id][hit_id] = 1
		if outfile:
			outfile.close()
			fasfile.close()
		seqs.close()
		if search_mode.endswith('-db'):
			os.remove('all.raw.fasta')
			os.remove('all.raw.fasta.idx')
		open(os.path.join('ok', step_name + '.ok'), 'w').close()
	return rep

# run MMseqs2 easy-search to map the reads
def mmseqs2(alns, rawdata, args):
	if os.path.isfile(os.path.join('ok', 'extract_reads.ok')):
		print("\nUsing the existing extracted read files")
		return 0
	for sp, fastxs in rawdata.items():
		for i in range(len(fastxs)):
			if os.path.isfile(os.path.join('ok', f"map_reads_{sp}_fastx{i+1}.ok")):
				print(f"\nUsing the existing MMseqs2 result file of target file {i+1} of {sp}")
			else:
				print(f"\nMapping the sequences from target file {i+1} of {sp}...")
				lib.check_programs(['mmseqs'])
				log = open(f"map_reads_{sp}_fastx{i+1}.log", 'w')
				cmd = ['mmseqs', 'easy-search', fastxs[i], os.path.join('ref_hmm', 'ref_profile'), f"map_reads_{sp}_fastx{i+1}.mmseqs.tsv", 'tmp', '--format-output', 'query,target,qstart,qend,tstart,tend,qaln,taln', '--threads', str(args.cpu)]
				cmd.extend(args.mmseqs_esearch_parameters)
				lib.runcmd(cmd, log, error=True)
				log.close()
				open(os.path.join('ok', f"map_reads_{sp}_fastx{i+1}.ok"), 'w').close()
	print("\nExtracting the target reads...")
	if os.path.isdir('map'):
		shutil.rmtree('map')
	os.mkdir('map')
	if os.path.exists('temp.idx'):
		os.remove('temp.idx')
	first2alnid = {}
	for line in open('first2alnid.tsv'):
		arr = line.rstrip().split("\t")
		first2alnid[arr[0]] = arr[1]
	if args.mol_type.startswith('dna'):
		frame = 'None'
	else:
		frame = 0
	total_counts = {}
	for sp, fastxs in rawdata.items():
		total_counts[sp] = 0
		for i in range(len(fastxs)):
			rep = {}
			for alnid in alns.keys():
				rep[alnid] = {}
			for line in open(f"map_reads_{sp}_fastx{i+1}.mmseqs.tsv"):
				arr = line.rstrip("\n").split("\t")
				alnid = first2alnid[arr[1]]
				outfile = open(os.path.join('map', alnid + '.hit_info.tsv'), 'a')
				outfile.write(f"{alnid}\t{arr[0]}_PhAlTag_{sp}_fastx{i+1}\t{arr[4]}\t{arr[5]}\t{arr[2]}\t{arr[3]}\tdirect\t{frame}\t{arr[7]}\t{arr[6]}\n")
				outfile.close()
				rep[alnid][arr[0]] = 1
			if fastxs[i].endswith('.gz'):
				seqs = lib.read_fastx(fastxs[i], args.file_format)
				for alnid, repinfo in rep.items():
					fasfile = open(os.path.join('map', alnid + '.targets.fas'), 'a')
					for seqid in repinfo.keys():
						fasfile.write(f">{seqid}_PhAlTag_{sp}_fastx{i+1}\n{seqs[seqid]}\n")
					fasfile.close()
				total_counts[sp] += len(seqs)
				del seqs
			else:
				seqs, file_format = lib.read_fastx(fastxs[i], args.file_format, return_iter=True)
				seqs = lib.read_fastx(fastxs[i], file_format + '_db', return_iter=True, temp_idx=True)
				for alnid, repinfo in rep.items():
					fasfile = open(os.path.join('map', alnid + '.targets.fas'), 'a')
					for seqid in repinfo.keys():
						fasfile.write(f">{seqid}_PhAlTag_{sp}_fastx{i+1}\n{str(seqs[seqid].seq)}\n")
					fasfile.close()
				total_counts[sp] += len(seqs)
				seqs.close()
				os.remove('temp.idx')
	for alnid in alns.keys():
		if not os.path.exists(os.path.join('map', alnid + '.hit_info.tsv')):
			open(os.path.join('map', alnid + '.targets.fas'), 'w').close()
			open(os.path.join('map', alnid + '.hit_info.tsv'), 'w').close()
	outfile = open('total_counts.tsv', 'w')
	for sp, count in total_counts.items():
		outfile.write(f"{sp}\t{count}\n")
	outfile.close()
	open('merged_redundance.txt', 'w').close()
	open(os.path.join('ok', 'extract_reads.ok'), 'w').close()

# map the reads to the HMM/CM (when not database) or emit an alignments from HMM/CM database, assemble the sequences and remove the contamination of an alignment
def assemble(group_name, cpu, codon_aln=None, outdir='.', tool='assemble.py', gencode=1, dna_codon_unknow=None, hmmtool='hmmsearch', parameters=[], sample_num=10, infernal_no_replace=False):
	log = open(os.path.join('map', group_name + '.assemble.log'), 'w')
	if hmmtool.endswith('search'):
		if os.path.isfile(os.path.join('ok', f"{hmmtool}_{group_name}.ok")):
			log.write(f"\nUsing the existing search result file of {group_name}\n")
		else:
			log.write(f"\nMapping the reads to reference alignment of {group_name}...\n")
			if hmmtool.startswith('hmm'):
				search_mode = 'hmmer'
				search_suffix = '.hmm'
			else:
				search_mode = 'infernal-sto'
				search_suffix = '.cm'
			mlog = open(os.path.join('map', f"{group_name}.{hmmtool}.log"), 'w')
			cmd = [hmmtool, '-o', os.path.join('map', f"{group_name}.{hmmtool}.txt"), '--tblout', os.path.join('map', f"{group_name}.{hmmtool}.tbl"), '--cpu', str(cpu)]
			cmd.extend(parameters)
			cmd.extend([os.path.join('ref_hmm', group_name + search_suffix), 'all.target.fasta'])
			ifcomplish = lib.runcmd(cmd, mlog, stdout=False)
			mlog.close()
			if not ifcomplish:
				log.write(f"\nError in {hmmtool} command, please check '{os.path.join('map', group_name + '.' + hmmtool + '.log')}'!\n")
				log.close()
				return 1, group_name
			open(os.path.join('ok', f"{hmmtool}_{group_name}.ok"), 'w').close()
		extract_reads(os.path.join('map', f"{group_name}.{hmmtool}.txt"), 'extract_reads_' + group_name, search_mode, infernal_no_replace, log)
	elif hmmtool.endswith('emit'):
		if os.path.isfile(os.path.join('ok', f"{hmmtool}_{group_name}.ok")):
			log.write(f"\nUsing the existing emitted reference alignment of {group_name} from HMM/CM database\n")
		else:
			log.write(f"\nEmitting reference alignments of {group_name} from HMM/CM database...\n")
			if hmmtool.startswith('hmm'):
				suffix = 'hmm'
			else:
				suffix = 'cm'
			iferror, nothing = hmmemit(os.path.join('ref_hmm', 'ref.hmm'), group_name, cpu, hmmtool, '.' + suffix, parameters, sample_num, fetch=suffix+'fetch')
			if iferror:
				log.write(f"\nError in {hmmtool} command, please check '{os.path.join('ref_hmm', group_name + '.log')}'!\n")
				log.close()
				return 1, group_name
			open(os.path.join('ok', f"{hmmtool}_{group_name}.ok"), 'w').close()
	if hmmtool == 'mmseqs':
		if codon_aln:
			codon_aln = lib.read_fastx(codon_aln, 'fasta')
			ref_file = os.path.join('ref_hmm', group_name + '.aa.fas')
		else:
			ref_file = os.path.join('ref_hmm', group_name + '.ref.fas')
		try:
			ref_aln = lib.read_fastx(ref_file, 'fasta')
			new_ref = {}
			codon_ref = {}
			for seqid in ref_aln.keys():
				ref_aln[seqid] = ref_aln[seqid].upper().replace('.', '-').replace('~', '-')
				new_ref[seqid] = ''
				codon_ref[seqid] = ''
			if '--match-mode' in parameters:
				match_mode = parameters[parameters.index('--match-mode') + 1]
			else:
				match_mode = '0'
			if '--match-ratio' in parameters:
				match_num = float(parameters[parameters.index('--match-ratio') + 1])
			else:
				match_num = 0.5
			match_num = match_num * len(ref_aln)
			first_seqstr = list(ref_aln.values())[0]
			for i in range(len(first_seqstr)):
				ifdel = False
				if match_mode == '0':
					if first_seqstr[i] == '-':
						ifdel = True
				else:
					nogap_num = 0
					for seqstr in ref_aln.values():
						if seqstr[i] != '-':
							nogap_num += 1
					# MMSeqs seems to remain the columns with non-gaps + 1 (consensus sequence)?
					if nogap_num + 1 < match_num:
						ifdel = True
				if not ifdel:
					for seqid in ref_aln.keys():
						new_ref[seqid] += ref_aln[seqid][i]
						if codon_aln:
							codon_ref[seqid] += codon_aln[seqid][(3*i):(3*i+3)]
			outfile = open(os.path.join('map', group_name + '.ref.fas'), 'w')
			for seqid, seqstr in new_ref.items():
				outfile.write(f">{seqid}\n{seqstr}\n")
			outfile.close()
			if codon_aln:
				outfile = open(os.path.join('map', group_name + '.codon.ref.fas'), 'w')
				for seqid, seqstr in codon_ref.items():
					outfile.write(f">{seqid}\n{seqstr.upper()}\n")
				outfile.close()
		except:
			log.write(f"\nError: fail to convert '{ref_file}' into '{os.path.join('map', group_name + '.ref.fas')}'!\n")
			log.close()
			return 1, group_name
	elif os.path.exists(os.path.join('ref_hmm', group_name + '.ref.fas')):
		outfile = open(os.path.join('map', group_name + '.ref.fas'), 'w')
		for line in open(os.path.join('ref_hmm', group_name + '.ref.fas')):
			outfile.write(line)
		outfile.close()
	else:
		try:
			codon_ref = {}
			if codon_aln:
				codon_aln0 = lib.read_fastx(codon_aln, 'fasta')
				# remove the sites of all gaps
				codon_aln = {}
				for seqid in codon_aln0.keys():
					codon_aln[seqid] = ''
				for i in range(len(list(codon_aln0.values())[0])):
					ifallgaps = True
					for seqstr in codon_aln0.values():
						if seqstr[i] != '-':
							ifallgaps = False
					if not ifallgaps:
						for seqid, seqstr in codon_aln0.items():
							codon_aln[seqid] += seqstr[i]
			alignment = AlignIO.read(os.path.join('ref_hmm', group_name + '.sto'), "stockholm")
			cons_pos = alignment._per_col_annotations['reference_annotation']
			outfile = open(os.path.join('map', group_name + '.ref.fas'), 'w')
			for seq in alignment:
				seqstr0 = str(seq.seq)
				seqstr = ''
				codon_ref[seq.id] = ''
				for i in range(len(seqstr0)):
					if cons_pos[i] != '.':
						seqstr += seqstr0[i]
						if codon_aln:
							codon_ref[seq.id] += codon_aln[seq.id][(3*i):(3*i+3)]
				seqstr = seqstr.upper().replace('.', '-').replace('~', '-')
				if hmmtool.startswith('cm') and not infernal_no_replace:
					seqstr = seqstr.replace('U', 'T')
				outfile.write(f">{seq.id}\n{seqstr}\n")
			outfile.close()
			if codon_aln:
				outfile = open(os.path.join('map', group_name + '.codon.ref.fas'), 'w')
				for seqid, seqstr in codon_ref.items():
					outfile.write(f">{seqid}\n{seqstr.upper()}\n")
				outfile.close()
		except:
			log.write(f"\nError: fail to convert '{os.path.join('ref_hmm', group_name + '.sto')}' into '{os.path.join('map', group_name + '.ref.fas')}'!\n")
			log.close()
			return 1, group_name
	cmd = [tool, group_name, outdir, str(cpu)]
	ifcomplish = lib.runcmd(cmd, log, stdout=False)
	if ifcomplish:
		if dna_codon_unknow:
			try:
				trans_seq(os.path.join('nt_out', group_name + '.fa'), os.path.join('aa_out', group_name + '.fa'), gencode, dna_codon_unknow)
			except:
				log.write(f"\nError: fail to translate '{os.path.join('nt_out', group_name + '.fa')}' into '{os.path.join('aa_out', group_name + '.fa')}'!\n")
				log.close()
				return 1, group_name
		log.close()
		return 0, group_name
	else:
		log.close()
		return 1, group_name

# map the reads to the HMMs/CMs of alignments (when not database) or emit reference alignments from HMM/CM database, assemble the sequences and remove the contamination of each alignment in multiprocess
def assemble_mp(alns, rawdata, args):
	outfile = open('parameters.config', 'w')
	for args_key, args_value in vars(args).items():
		outfile.write(f"{args_key}:\t{args_value}\n")
	outfile.write("species:\t['" + ', '.join(rawdata.keys()) + "']\n")
	outfile.close()
	if not os.path.isdir('stat_info'):
		os.mkdir('stat_info')
	if not os.path.isdir('nt_out'):
		os.mkdir('nt_out')
	if args.mol_type != 'dna' and not os.path.isdir('aa_out'):
		os.mkdir('aa_out')

	if args.search_mode.endswith('-db'):
		print("\nGenerating the reference-based alignments through assembling, removing foreign sequences and cross contamination, and printing the new alignments...")
	else:
		print("\nGenerating the reference-based alignments through searching, assembling, removing foreign sequences and cross contamination, and printing the new alignments...")
	np = args.parallel
	ncpu = int(args.cpu / np)
	args_list = []
	kwds = {'outdir': args.output, 'tool': lib.check_assemble()}
	if args.mol_type == 'dna_codon':
		kwds['gencode'] = args.gencode
		kwds['dna_codon_unknow'] = args.unknow_prot
	usedcpu = ncpu * np
	if args.search_mode.endswith('-db'):
		getsize_dir = 'map'
		getsize_suffix = '.hit_info.tsv'
	else:
		getsize_dir = 'ref_hmm'
		if os.path.exists(os.path.join('ref_hmm', list(alns.keys())[0] + '.ref.fas')):
			getsize_suffix = '.ref.fas'
		elif os.path.exists(os.path.join('ref_hmm', list(alns.keys())[0] + '.sto')):
			getsize_suffix = '.sto'
		else:
			getsize_suffix = '.aa.fas'
	for group_name, alnfile in sorted(alns.items(), key=lambda x : os.path.getsize(os.path.join(getsize_dir, x[0] + getsize_suffix)), reverse=True):
		if usedcpu < args.cpu:
			if args.mol_type == 'codon' and not args.no_ref:
				args_list.append((group_name, ncpu+1, alnfile))
			else:
				args_list.append((group_name, ncpu+1))
			usedcpu += 1
		else:
			if args.mol_type == 'codon' and not args.no_ref:
				args_list.append((group_name, ncpu, alnfile))
			else:
				args_list.append((group_name, ncpu))
	if args.search_mode.startswith('mmseqs'):
		kwds['hmmtool'] = 'mmseqs'
		kwds['parameters'] = args.mmseqs_msa2profile_parameters
	elif args.search_mode.endswith('-db'):
		if args.search_mode.startswith('hmmer'):
			kwds['hmmtool'] = 'hmmemit'
			kwds['parameters'] = args.hmmemit_parameters
			kwds['sample_num'] = args.emit_sample_num
			lib.check_programs(['hmmfetch'])
		else:
			kwds['hmmtool'] = 'cmemit'
			kwds['parameters'] = args.cmemit_parameters
			kwds['sample_num'] = args.emit_sample_num
			lib.check_programs(['cmfetch'])
	else:
		if os.path.exists('all.raw.fasta'):
			# build the index of Bio.SeqIO index_db
			lib.read_fastx('all.raw.fasta', 'fasta_db', return_iter=True).close()
		if args.search_mode.startswith('hmmer'):
			kwds['hmmtool'] = 'hmmsearch'
			kwds['parameters'] = args.hmmsearch_parameters
		else:
			kwds['hmmtool'] = 'cmsearch'
			kwds['parameters'] = args.cmsearch_parameters
			kwds['infernal_no_replace'] = args.infernal_no_replace
	lib.check_programs([kwds['hmmtool']])
	iferrors = lib.run_mp(assemble, args_list, np, kwds=kwds)
	errors = []
	for iferror in iferrors:
		if iferror[0] == 1:
			errors.append(iferror[1])
	if errors:
		print("\nError in assemble: {}\nplease check map/alignment_name.assemble.log!".format(', '.join(errors)))
		sys.exit(1)
	if not args.search_mode.startswith('mmseqs') and not args.search_mode.endswith('-db') and os.path.exists('all.raw.fasta'):
		os.remove('all.target.fasta')
		os.remove('all.raw.fasta')
		os.remove('all.raw.fasta.idx')

# concatenate the short alignments if split the reference alignements
def concatenate(dirname, alns, fill='-'):
	seqs = {}
	total_len = 0
	for aln in alns.keys():
		aln_seqs = lib.read_fastx(os.path.join(dirname, aln + '.fa'), 'fasta')
		for seqid, seqstr in aln_seqs.items():
			if seqid.endswith('.' + aln + '.1'):
				seqid = seqid.replace('.' + aln + '.1', '')
			if seqs.get(seqid) is None:
				seqs[seqid] = fill * total_len + seqstr
			else:
				seqs[seqid] += seqstr
			aln_len = len(seqstr)
		for seqid in seqs.keys():
			if aln_seqs.get(seqid) is None and aln_seqs.get(seqid + '.' + aln + '.1') is None:
				seqs[seqid] += fill * aln_len
		total_len += aln_len
	outfile = open(os.path.join(dirname, 'aln.concatenated.fa'), 'w')
	for seqid, seqstr in seqs.items():
		outfile.write(">{}\n{}\n".format(seqid, seqstr))
	outfile.close()
