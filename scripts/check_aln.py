#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

def read_fasta(fasta):
	seqs = {}
	seqid = ''
	for line in open(fasta):
		line = line.rstrip()
		if line.startswith('>'):
			arr = line.split(" ")
			seqid = arr[0].lstrip('>')
			seqs[seqid] = ''
		elif seqs.get(seqid) is not None:
			seqs[seqid] += line.upper()
	return seqs

if len(sys.argv) == 1 or sys.argv[1] == '-h':
    print("Usage: {} input_dir output_dir(default='none') aver_freq_per_site(default=0.75) gap_symbol(default='-') start_end_no_gap_number(>=1)_or_percent(<1)(default=0.6) fasta_suffix(default='.fa')".format(sys.argv[0]))
    sys.exit(0)
if len(sys.argv) < 2:
    print("Error: options < 1!\nUsage: {} input_dir output_dir(default='none') aver_freq_per_site(default=0.75) gap_symbol(default='-') start_end_no_gap_number(>=1)_or_percent(<1)(default=0.6) fasta_suffix(default='.fa')".format(sys.argv[0]))
    sys.exit(1)
outdir = None
if len(sys.argv) > 2:
	if sys.argv[2].lower() != 'none':
		outdir = sys.argv[2]
freq = 0.75
if len(sys.argv) > 3:
	freq = float(sys.argv[3])
gap = '-'
if len(sys.argv) > 4:
	gap = sys.argv[4]
pgap = 0.6
if len(sys.argv) > 5:
	pgap = float(sys.argv[5])
suffix = '.fa'
if len(sys.argv) > 6:
	suffix = sys.argv[6]

if outdir:
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

files = os.listdir(sys.argv[1])
for filename in files:
	if not filename.endswith(suffix):
		continue
	seqs = read_fasta(os.path.join(sys.argv[1], filename))
	unaln = []
	start = None
	end = None
	if pgap < 1:
		ngap = pgap * len(seqs)
	else:
		ngap = pgap
	for i in range(len(list(seqs.values())[0])):
			n = 0
			for seqstr in seqs.values():
				if seqstr[i] != gap:
					n += 1
			if n >= ngap:
				if start is None:
					start = i
				end = i
	print(f"\n{filename}, length: {len(list(seqs.values())[0])}, valid regions: {start}-{end}")
	if start is None or end is None:
		for seqid in seqs.keys():
			print(f"{seqid}, None-None: -Inf < {freq}, unaligned!")
			unaln.append(seqid)
	else:
		scores = {}
		for i in range(start, end + 1):
			scores[i] = {}
			for seqstr in seqs.values():
				scores[i].setdefault(seqstr[i], 0)
				scores[i][seqstr[i]] += 1
		for seqid, seqstr in seqs.items():
			seq_start = None
			seq_end = None
			for i in range(start, end + 1):
				if seqstr[i] != gap:
					if seq_start is None:
						 seq_start = i
					seq_end = i
			if seq_start is None or seq_end is None:
				print(f"{seqid}, {seq_start}-{seq_end}: -Inf < {freq}, unaligned!")
				unaln.append(seqid)
				continue
			aver_freq = 0
			for i in range(seq_start, seq_end + 1):
				aver_freq += (scores[i][seqstr[i]] / len(seqs))
			aver_freq = aver_freq / (seq_end - seq_start + 1)
			if aver_freq < freq:
				print(f"{seqid}, {seq_start}-{seq_end}: {aver_freq} < {freq}, unaligned!")
				unaln.append(seqid)
			else:
				print(f"{seqid}, {seq_start}-{seq_end}: {aver_freq}")
	if outdir:
		outfile = open(os.path.join(outdir, filename), 'w')
		for seqid, seqstr in seqs.items():
			if seqid not in unaln:
				outfile.write(">{}\n{}\n".format(seqid, seqstr))
		outfile.close()
