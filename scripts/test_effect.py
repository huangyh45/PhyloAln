#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

def read_fasta(fasta, target, sep='.', select_list=None):
	seqs = {}
	if select_list:
		seqid = None
		for line in open(fasta):
			line = line.rstrip()
			if line.startswith('>'):
				arr = line.split(" ")
				seqid = arr[0].lstrip('>')
				if seqid in select_list or seqid.split(sep)[0] in select_list:
					seqs[seqid] = ''
			elif seqs.get(seqid) is not None:
				seqs[seqid] += line.upper()
	seqid = None
	seqstr = ''
	for line in open(fasta):
		line = line.rstrip()
		if line.startswith('>'):
			if seqstr:
				break
			arr = line.split(" ")
			seqid = arr[0].lstrip('>')
			if seqid != target and not seqid.startswith(target + sep):
				seqid = None
		elif seqid:
			seqstr += line.upper()
	# debug
	# print(seqid, seqstr)
	# print(seqs)
	# input()
	return seqstr, seqs

if len(sys.argv) == 1 or sys.argv[1] == '-h':
    print("Usage: {} reference_dir:ref_species_or_seq_name target_dir:target_species_or_seq_name output_tsv unknown_symbol(default='N') separate(default='.') fasta_suffix(default='.fa') selected_species_or_sequences(separated by comma)".format(sys.argv[0]))
    sys.exit(0)
if len(sys.argv) <= 3:
    print("Error: options < 3!\nUsage: {} reference_dir:ref_species_or_seq_name target_dir:target_species_or_seq_name output_tsv unknown_symbol(default='N') separate(default='.') fasta_suffix(default='.fa') selected_species_or_sequences(separated by comma)".format(sys.argv[0]))
    sys.exit(1)
unknow = 'N'
if len(sys.argv) > 4:
	unknow = sys.argv[4]
sep = '.'
if len(sys.argv) > 5:
	sep = sys.argv[5]
suffix = '.fa'
if len(sys.argv) > 6:
	suffix = sys.argv[6]
select_list = None
if len(sys.argv) > 7:
	select_list = sys.argv[7].split(',')

ref_dir, ref_sp = sys.argv[1].split(':')
target_dir, target_sp = sys.argv[2].split(':')
outfile = open(sys.argv[3], 'w')
outfile.write("file\treference length\ttarget length\tnident\tcompleteness\tpident\n")
total_rlen = 0
total_tlen = 0
total_nident = 0
files = os.listdir(ref_dir)
for filename in files:
	if not filename.endswith(suffix):
		continue
	ref_seq, ref_sel = read_fasta(os.path.join(ref_dir, filename), ref_sp, sep, select_list)
	target_seq, target_sel = read_fasta(os.path.join(target_dir, filename), target_sp, sep, select_list)
	indexes = []
	if select_list:
		i = 0
		j = 0
		while i < len(list(ref_sel.values())[0]):
			if j >= len(list(target_sel.values())[0]):
				if_match = False
			else:
				if_match = True
				for sel in ref_sel.keys():
					# to use species
					# sp = sel.split(sep)[0]
					if ref_sel[sel][i] != target_sel[sel][j]:
						if_match = False
						break
			if if_match:
				j += 1
			else:
				indexes.append(i)
			i += 1
		print(filename)
		str_list = list(ref_seq)
		for i in reversed(indexes):
			str_list.pop(i)
		ref_seq = ''.join(str_list)
	start = None
	for i in range(len(ref_seq)):
		if ref_seq[i] not in [unknow, '-']:
			start = i
			break
	end = None
	for i in range(len(ref_seq)-1, -1, -1):
		if ref_seq[i] not in [unknow, '-']:
			end = i
			break
	if start is None or end is None:
		# continue
		print("\nError: invalid reference in '{}': '{}'!\n".format(filename, ref_seq))
		sys.exit(1)
	if not target_seq:
		rlen = 0
		for i in range(start, end+1):
			if ref_seq[i] != unknow:
				rlen += 1
		total_rlen += rlen
		outfile.write("{}\t{}\t0\t0\t0\t0\n".format(filename, rlen))
		continue
	rlen = 0
	tlen = 0
	nident = 0
	# print(filename, start, end)
	for i in range(start, end+1):
		if ref_seq[i] == unknow:
			continue
		rlen += 1
		if target_seq[i] == unknow:
			continue
		tlen += 1
		if target_seq[i] == ref_seq[i]:
			nident += 1
	total_rlen += rlen
	total_tlen += tlen
	total_nident += nident
	outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(filename, rlen, tlen, nident, int(10000*tlen/rlen+0.5)/100, int(10000*nident/tlen+0.5)/100))
outfile.write("total\t{}\t{}\t{}\t{}\t{}\n".format(total_rlen, total_tlen, total_nident, int(10000*total_tlen/total_rlen+0.5)/100, int(10000*total_nident/total_tlen+0.5)/100))
outfile.close()
