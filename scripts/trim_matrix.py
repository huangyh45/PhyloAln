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
			seqs[seqid] += line
	return seqs

if len(sys.argv) == 1 or sys.argv[1] == '-h':
    print("Usage: {} input_dir output_dir unknown_symbol(default='X') known_number(>=1)_or_percent(<1)_for_columns(default=0.5) known_number(>=1)_or_percent(<1)_for_rows(default=0) fasta_suffix(default='.fa')".format(sys.argv[0]))
    sys.exit(0)
if len(sys.argv) < 3:
    print("Error: options < 2!\nUsage: {} input_dir output_dir unknown_symbol(default='X') known_percent_for_columns(default=50) known_percent_for_rows(default=0) fasta_suffix(default='.fa')".format(sys.argv[0]))
    sys.exit(1)
unknow = 'X'
if len(sys.argv) > 3:
	unknow = sys.argv[3]
pcol = 0.5
if len(sys.argv) > 4:
	pcol = float(sys.argv[4])
prow = 0
if len(sys.argv) > 5:
	prow = float(sys.argv[5])
suffix = '.fa'
if len(sys.argv) > 6:
	suffix = sys.argv[6]

if not os.path.isdir(sys.argv[2]):
	os.mkdir(sys.argv[2])

files = os.listdir(sys.argv[1])
for filename in files:
	if not filename.endswith(suffix):
		continue
	seqs = read_fasta(os.path.join(sys.argv[1], filename))
	if pcol > 0:
		if pcol < 1:
			ncol = pcol * len(seqs)
		else:
			ncol = pcol
		newseqs = {}
		for seqid in seqs.keys():
			newseqs[seqid] = ''
		for i in range(len(list(seqs.values())[0])):
			n = 0
			for seqstr in seqs.values():
				if seqstr[i] != unknow:
					n += 1
			if n >= ncol:
				for seqid, seqstr in seqs.items():
					newseqs[seqid] += seqstr[i]
		seqs = newseqs
	if prow > 0:
		if prow < 1:
			nrow = prow * len(list(seqs.values())[0])
		else:
			nrow = prow
		for seqid, seqstr in list(seqs.items()):
			n = 0
			for base in seqstr:
				if base != unknow:
					n += 1
			if n < nrow:
				print("Removing {} from {}: known sites {}/{} = {}%".format(seqid, filename, n, len(seqstr), int(n / len(seqstr) * 10000 + 0.5) / 100))
				seqs.pop(seqid)
	outfile = open(os.path.join(sys.argv[2], filename), 'w')
	for seqid, seqstr in seqs.items():
		outfile.write(">{}\n{}\n".format(seqid, seqstr))
	outfile.close()

