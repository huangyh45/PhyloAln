#!/usr/bin/env python3

import sys
import os
import warnings
from Bio import SearchIO, BiopythonWarning
from Bio.Seq import reverse_complement, translate
try:
    import library as lib
except ImportError:
    import lib.library as lib

# split the reference alignments into short alignments
def split_ref(alnfile, split_len):
	print("Splitting the reference alignment into short alignments...")
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
		i += 1
	return alns

# output the splited sequences to different FASTA file by single CPU
def output_fasta_percpu(seq_list, output_fasta, fastx_num=1, moltype='dna', gencode=1, outfile0=None, split_len=None, split_slide=None, no_reverse=False, low_mem_iter=None, low_mem_format=None):
	if outfile0 is None:
		outfile = open(output_fasta, 'w')
	else:
		outfile = outfile0
	if low_mem_iter is not None:
		# low-memory mode: directly read the file instead of reading the store
		seqid = None
		seqstr = None
		count = 0
		line_num = 0
		for line in low_mem_iter:
			line = line.rstrip()
			if low_mem_format == 'fastq' and line.startswith('@') and line_num % 4 == 0:
				seqid = line.replace('@', '', 1).replace(' ', '_')
			elif low_mem_format == 'fasta' and line.startswith('>'):
				if seqid:
					# prepare the single sequence
					output_fasta_percpu([(seqid, seqstr)], output_fasta, fastx_num, moltype, gencode, outfile, split_len=split_len, split_slide=split_slide, no_reverse=no_reverse)
					count += 1
				arr = line.split(" ")
				seqid = arr[0].lstrip('>')
				seqstr = ''
			elif seqid:
				if low_mem_format == 'fastq':
					# prepare the single sequence
					output_fasta_percpu([(seqid, line)], output_fasta, fastx_num, moltype, gencode, outfile, split_len=split_len, split_slide=split_slide, no_reverse=no_reverse)
					count += 1
					seqid = None
				else:
					seqstr += line
			line_num += 1
		if low_mem_format == 'fasta' and seqstr:
			output_fasta_percpu([(seqid, seqstr)], output_fasta, fastx_num, moltype, gencode, outfile, split_len=split_len, split_slide=split_slide, no_reverse=no_reverse)
			count += 1
		# low-memory mode for this FASTA/FASTQ file ends here and return the read count
		return count
	for seqid, seqstr0 in seq_list:
		if split_len:
			# split the sequences into short sequences with length of split_len
			i = 0
			while i + split_len < len(seqstr0):
				output_fasta_percpu([('_'.join([seqid, 'split', str(i+1), str(i+split_len)]), seqstr0[i:i+split_len])], output_fasta, fastx_num, moltype, gencode, outfile, no_reverse=no_reverse)
				i += split_slide
			output_fasta_percpu([('_'.join([seqid, 'split', str(i+1), str(len(seqstr0))]), seqstr0[i:])], output_fasta, fastx_num, moltype, gencode, outfile, no_reverse=no_reverse)
		elif moltype.startswith('dna'):
			outfile.write(">{}_fastx{}\n{}\n".format(seqid, fastx_num, seqstr0))
			if not no_reverse:
				outfile.write(">{}_fastx{}_rev\n{}\n".format(seqid, fastx_num, reverse_complement(seqstr0)))
		else:
			# supress the warnings due to the last incomplete codons when translate the sequences
			with warnings.catch_warnings():
				warnings.simplefilter('ignore', BiopythonWarning)
				for j in [1,2,3]:
					seqstr = seqstr0[(j-1):]
					outfile.write(">{}_fastx{}_pos{}\n{}\n".format(seqid, fastx_num, j, translate(seqstr, table=gencode)))
				if not no_reverse:
					seqstr0 = reverse_complement(seqstr0)
					for j in [1,2,3]:
						seqstr = seqstr0[(j-1):]
						outfile.write(">{}_fastx{}_pos{}rev\n{}\n".format(seqid, fastx_num, j, translate(seqstr, table=gencode)))
	if outfile0 is None:
		outfile.close()

# convert the raw FASTQ/FASTA files to (translated) FASTA format
def fastx2fasta(fastxs, fasta, file_format='guess', cpu=8, moltype='dna', gencode=1, split_len=None, split_slide=None, no_reverse=False, low_mem=False, output=True):
	if output:
		outfile = open(fasta, 'w')
	all_seqs = []
	total_count = 0
	for i in range(len(fastxs)):
		if low_mem:
			fastx_iter, file_format = lib.read_fastx(fastxs[i], file_format, low_mem=True)
			all_seqs.append([fastxs[i], file_format])
			seqs = {}
		else:
			seqs = lib.read_fastx(fastxs[i], file_format)
			all_seqs.append(seqs)
			total_count += len(seqs)
		if output:
			if not low_mem and cpu > 1 and (split_len or len(seqs) > 10000):
				# run in multiprocess when too much sequences
				print("Binning the reads in '{}' into {} parts and preparing in multiprocess...".format(fastxs[i], cpu))
				nseq = int(len(seqs) / cpu)
				if nseq * cpu < len(seqs):
					nseq += 1
				kwds = {'fastx_num': i+1, 'moltype': moltype, 'gencode': gencode, 'split_len': split_len, 'split_slide': split_slide, 'no_reverse': no_reverse}
				args_list = []
				for j in range(cpu-1):
					args_list.append((list(seqs.items())[(nseq*j):(nseq*(j+1))], fasta + '.' + str(j)))
				args_list.append((list(seqs.items())[(nseq*(cpu-1)):], fasta + '.' + str(cpu-1)))
				lib.run_mp(output_fasta_percpu, args_list, cpu, kwds=kwds)
				for j in range(cpu):
					for line in open(fasta + '.' + str(j)):
						outfile.write(line)
					os.remove(fasta + '.' + str(j))
			elif low_mem:
				# output using a low-memory method and obtain the read count
				total_count += output_fasta_percpu([], fasta, i+1, moltype, gencode, outfile, split_len=split_len, split_slide=split_slide, no_reverse=no_reverse, low_mem_iter=fastx_iter, low_mem_format=file_format)
			else:
				# directly output by single cpu
				output_fasta_percpu(list(seqs.items()), fasta, i+1, moltype, gencode, outfile, split_len=split_len, split_slide=split_slide, no_reverse=no_reverse)
	if output:
		outfile.close()
	return all_seqs, total_count

# construct HMM file by hmmbuild
def hmmbuild(alnfile, group_name, cpu=1, moltype='dna', gencode=1, parameters=[]):
	log = open(os.path.join('ref_hmm', group_name + '.log'), 'w')
	if moltype == 'codon':
		lib.trans_seq(alnfile, os.path.join('ref_hmm', group_name + '.aa.fas'), gencode=gencode)
		cmd = ['hmmbuild', '-O', os.path.join('ref_hmm', group_name + '.sto'), '--cpu', str(cpu)]
		cmd.extend(parameters)
		cmd.extend([os.path.join('ref_hmm', group_name + '.hmm'), os.path.join('ref_hmm', group_name + '.aa.fas')])
	else:
		cmd = ['hmmbuild', '-O', os.path.join('ref_hmm', group_name + '.sto'), '--cpu', str(cpu)]
		cmd.extend(parameters)
		cmd.extend([os.path.join('ref_hmm', group_name + '.hmm'), alnfile])
	ifcomplish = lib.runcmd(cmd, log, stdout=False)
	log.close()
	if ifcomplish:
		return 0, group_name
	else:
		return 1, group_name

# HMM search by HMMER
def hmmsearch(group_name, species, cpu=1, parameters=[]):
	if not os.path.isfile(os.path.join('ok', "map_{}_{}.ok".format(species, group_name))):
		fasta = species + '.temp.fasta'
		log = open(os.path.join('map_' + species, group_name + '.log'), 'w')
		cmd = ['hmmsearch', '-o', os.path.join('map_' + species, group_name + '.txt'), '--tblout', os.path.join('map_' + species, group_name + '.tbl'), '--cpu', str(cpu)]
		cmd.extend(parameters)
		cmd.extend([os.path.join('ref_hmm', group_name + '.hmm'), fasta])
		ifcomplish = lib.runcmd(cmd, log, stdout=False)
		log.close()
		if ifcomplish:
			open(os.path.join('ok', "map_{}_{}.ok".format(species, group_name)), 'w').close()
			return 0, group_name
		else:
			return 1, group_name
	return 0, group_name

# prepare the HMMs of alignments
def prepare_ref(alns, cpu=8, np=8, moltype='dna', gencode=1, parameters=[]):
	lib.check_programs(['hmmbuild'])
	np = min(np, len(alns), cpu)
	ncpu = int(cpu / np)
	
	print("\nPreparing the reference alignments...")
	if os.path.isfile(os.path.join('ok', 'prepare_alignments.ok')):
		print("\nUsing the existing hmm files in directory 'ref_hmm'")
	else:
		print("\nBuilding HMMs for mapping...")
		if not os.path.isdir('ref_hmm'):
			os.mkdir('ref_hmm')
		args_list = []
		kwds = {'moltype': moltype, 'gencode': gencode, 'parameters': parameters}
		usedcpu = ncpu * np
		for group_name, alnfile in alns.items():
			if usedcpu < cpu:
				args_list.append((alnfile, group_name, ncpu+1))
				usedcpu += 1
			else:
				args_list.append((alnfile, group_name, ncpu))
		iferrors = lib.run_mp(hmmbuild, args_list, np, kwds=kwds)
		errors = []
		for iferror in iferrors:
			if iferror[0] == 1:
				errors.append(iferror[1])
		if errors:
			print("\nError in hmmbuild commands: {}".format(', '.join(errors)))
			sys.exit(1)
		open(os.path.join('ok', 'prepare_alignments.ok'), 'w').close()

# map the reads to the HMMs of alignments
def map_reads(alns, species, fastxs, file_format='guess', cpu=8, np=8, moltype='dna', gencode=1, split_len=None, split_slide=None, no_reverse=False, low_mem=False, parameters=[]):
	lib.check_programs(['hmmsearch'])
	if os.path.isfile(os.path.join('ok', "prepare_{}.ok".format(species))):
		print("\nUsing the existing temp FASTA file of {}".format(species))
		all_seqs, total_reads = fastx2fasta(fastxs, species + '.temp.fasta', file_format, cpu=cpu, moltype=moltype, gencode=gencode, split_len=split_len, split_slide=split_slide, no_reverse=no_reverse, low_mem=low_mem, output=False)
	else:
		print("\nPreparing the reads of {}...".format(species))
		all_seqs, total_reads = fastx2fasta(fastxs, species + '.temp.fasta', file_format, cpu=cpu, moltype=moltype, split_len=split_len, split_slide=split_slide, no_reverse=no_reverse, low_mem=low_mem, gencode=gencode)
		open(os.path.join('ok', "prepare_{}.ok".format(species)), 'w').close()

	if not os.path.isdir('map_' + species):
		os.mkdir('map_' + species)
	np = min(np, len(alns), cpu)
	ncpu = int(cpu / np)
	print("\nMapping the reads of {} to reference alignments...".format(species))
	args_list = []
	kwds = {'parameters': parameters}
	usedcpu = ncpu * np
	for group_name in alns.keys():
		if usedcpu < cpu:
			args_list.append((group_name, species, ncpu+1))
			usedcpu += 1
		else:
			args_list.append((group_name, species, ncpu))
	iferrors = lib.run_mp(hmmsearch, args_list, np, kwds=kwds)
	errors = []
	for iferror in iferrors:
		if iferror[0] == 1:
			errors.append(iferror[1])
	if errors:
		print("\nError in hmmsearch commands: {}".format(', '.join(errors)))
		sys.exit(1)
	if os.path.exists(species + '.temp.fasta'):
		os.remove(species + '.temp.fasta')
	return all_seqs, total_reads

# read the information of HMMER results
def read_hmmer(hmmtxt):
	hmmresults = []
	qresults =  SearchIO.parse(hmmtxt, 'hmmer3-text')
	for qresult in qresults:
		for hit in qresult:
			for HSP in hit:
				for HSPfrag in HSP:
					hmmresults.append(HSPfrag)
	return hmmresults

# extract the target reads from the raw data
def extract_reads(alns, species, all_seqs, moltype='dna', split_len=None, low_mem=False):
	all_hmmres = {}
	target_seqids = {}
	for group_name in alns.keys():
		hmmresults = read_hmmer(os.path.join('map_' + species, group_name + '.txt'))
		all_hmmres[group_name] = hmmresults
		outfile = open(os.path.join('map_' + species, group_name + '.targets.fa'), 'w')
		for hmmresult in hmmresults:
			if moltype.startswith('dna'):
				hitid = hmmresult.hit_id
				if hitid.endswith('_rev'):
					hitid = hitid.replace('_rev', '')
				seqid = '_'.join(hitid.split('_')[:-1])
				fastx_num = hitid.split('_')[-1].replace('fastx', '')
			else:
				seqid = '_'.join(hmmresult.hit_id.split('_')[:-2])
				fastx_num = hmmresult.hit_id.split('_')[-2].replace('fastx', '')
			if low_mem:
				if target_seqids.get(fastx_num) is None:
					target_seqids[fastx_num] = {}
				if split_len:
					seqid0, start_end = seqid.split('_split_')
					if target_seqids[fastx_num].get(seqid0) is None:
						target_seqids[fastx_num][seqid0] = {}
					if target_seqids[fastx_num][seqid0].get(group_name) is None:
						target_seqids[fastx_num][seqid0][group_name] = []
					target_seqids[fastx_num][seqid0][group_name].append(start_end.split('_'))
				else:
					if target_seqids[fastx_num].get(seqid) is None:
						target_seqids[fastx_num][seqid] = []
					target_seqids[fastx_num][seqid].append(group_name)
			elif split_len:
				seqid, start_end = seqid.split('_split_')
				start, end = start_end.split('_')
				outfile.write(">{}\n{}\n".format('_'.join([seqid, 'split', start, end, 'fastx' + fastx_num]), all_seqs[int(fastx_num)-1][seqid][int(start)-1:int(end)]))
			else:
				outfile.write(">{}\n{}\n".format(seqid + '_fastx' + fastx_num, all_seqs[int(fastx_num)-1][seqid]))
		outfile.close()
	if low_mem:
		# low-memory mode: extract the reads from the file instead of the store
		for fastx_num, seqids in target_seqids.items():
			fastx_iter, file_format = lib.read_fastx(all_seqs[int(fastx_num)-1][0], all_seqs[int(fastx_num)-1][1], low_mem=True)
			seqid = None
			seqstr = None
			line_num = 0
			for line in fastx_iter:
				line = line.rstrip()
				if file_format == 'fastq' and line.startswith('@') and line_num % 4 == 0:
					seqid = line.replace('@', '', 1).replace(' ', '_')
				elif file_format == 'fasta' and line.startswith('>'):
					if seqid is not None and seqids.get(seqid) is not None:
						if split_len:
							for group_name, start_ends in seqids[seqid].items():
								outfile = open(os.path.join('map_' + species, group_name + '.targets.fa'), 'a')
								for start_end in start_ends:
									outfile.write(">{}\n{}\n".format('_'.join([seqid, 'split', '_'.join(start_end), 'fastx' + fastx_num]), seqstr[int(start_end[0])-1:int(start_end[1])]))
								outfile.close()
						else:
							for group_name in seqids[seqid]:
								outfile = open(os.path.join('map_' + species, group_name + '.targets.fa'), 'a')
								outfile.write(">{}\n{}\n".format(seqid + '_fastx' + fastx_num, seqstr))
								outfile.close()
					arr = line.split(" ")
					seqid = arr[0].lstrip('>')
					seqstr = ''
				elif seqid:
					if file_format == 'fastq':
						if seqids.get(seqid) is not None:
							if split_len:
								for group_name, start_ends in seqids[seqid].items():
									outfile = open(os.path.join('map_' + species, group_name + '.targets.fa'), 'a')
									for start_end in start_ends:
										outfile.write(">{}\n{}\n".format('_'.join([seqid, 'split', '_'.join(start_end), 'fastx' + fastx_num]), line[int(start_end[0])-1:int(start_end[1])]))
									outfile.close()
							else:
								for group_name in seqids[seqid]:
									outfile = open(os.path.join('map_' + species, group_name + '.targets.fa'), 'a')
									outfile.write(">{}\n{}\n".format(seqid + '_fastx' + fastx_num, line))
									outfile.close()
						seqid = None
					else:
						seqstr += line
				line_num += 1
			if file_format == 'fasta' and seqstr and seqids.get(seqid) is not None: 
				if split_len:
					for group_name, start_ends in seqids[seqid].items():
						outfile = open(os.path.join('map_' + species, group_name + '.targets.fa'), 'a')
						for start_end in start_ends:
							outfile.write(">{}\n{}\n".format('_'.join([seqid, 'split', '_'.join(start_end), 'fastx' + fastx_num]), seqstr[int(start_end[0])-1:int(start_end[1])]))
						outfile.close()
				else:
					for group_name in seqids[seqid]:
						outfile = open(os.path.join('map_' + species, group_name + '.targets.fa'), 'a')
						outfile.write(">{}\n{}\n".format(seqid + '_fastx' + fastx_num, seqstr))
						outfile.close()
	return all_hmmres
