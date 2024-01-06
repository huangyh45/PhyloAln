#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

if len(sys.argv) == 1 or sys.argv[1] == '-h':
    print("Usage: {} output_dir PhyloAln_dir1 PhyloAln_dir2 (PhyloAln_dir3 ...)".format(sys.argv[0]))
    sys.exit(0)
if len(sys.argv) < 3:
    print("Error: options < 2!\nUsage: {} output_dir PhyloAln_dir1 PhyloAln_dir2 (PhyloAln_dir3 ...)".format(sys.argv[0]))
    sys.exit(1)

os.mkdir(sys.argv[1])
if os.path.isdir(os.path.join(sys.argv[2], 'aa_out')):
	os.mkdir(os.path.join(sys.argv[1], 'aa_out'))
if os.path.isdir(os.path.join(sys.argv[2], 'nt_out')):
	os.mkdir(os.path.join(sys.argv[1], 'nt_out'))

aafiles = {}
ntfiles = {}
if os.path.isdir(os.path.join(sys.argv[2], 'aa_out')):
	for dirname in sys.argv[2:]:
		filenames = os.listdir(os.path.join(dirname, 'aa_out'))
		for filename in filenames:
			if not filename.endswith('.fa'):
				continue
			if aafiles.get(filename) is None:
				aafiles[filename] = []
			aafiles[filename].append(os.path.join(dirname, 'aa_out', filename))
if os.path.isdir(os.path.join(sys.argv[2], 'nt_out')):
	for dirname in sys.argv[2:]:
		filenames = os.listdir(os.path.join(dirname, 'nt_out'))
		for filename in filenames:
			if not filename.endswith('.fa'):
				continue
			if ntfiles.get(filename) is None:
				ntfiles[filename] = []
			ntfiles[filename].append(os.path.join(dirname, 'nt_out', filename))

for aafile, filenames in aafiles.items():
	seqs = {}
	for filename in filenames:
		seqid = ''
		for line in open(filename):
			line = line.rstrip()
			if line.startswith('>'):
				arr = line.split(" ")
				seqid = arr[0].lstrip('>')
				seqs[seqid] = ''
			elif seqs.get(seqid) is not None:
				seqs[seqid] += line
	outfile = open(os.path.join(sys.argv[1], 'aa_out', aafile), 'w')
	for seqid, seqstr in seqs.items():
		outfile.write(">{}\n{}\n".format(seqid, seqstr))
	outfile.close()
for ntfile, filenames in ntfiles.items():
	seqs = {}
	for filename in filenames:
		seqid = ''
		for line in open(filename):
			line = line.rstrip()
			if line.startswith('>'):
				arr = line.split(" ")
				seqid = arr[0].lstrip('>')
				seqs[seqid] = ''
			elif seqs.get(seqid) is not None:
				seqs[seqid] += line
	outfile = open(os.path.join(sys.argv[1], 'nt_out', ntfile), 'w')
	for seqid, seqstr in seqs.items():
		outfile.write(">{}\n{}\n".format(seqid, seqstr))
	outfile.close()

