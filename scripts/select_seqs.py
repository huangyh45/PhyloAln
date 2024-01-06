#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

if len(sys.argv) == 1 or sys.argv[1] == '-h':
    print("Usage: {} input_dir selected_species_or_sequences(separated by comma) output_dir fasta_suffix(default='.fa') separate_symbol(default='.') if_list_for_exclusion(default=no)".format(sys.argv[0]))
    sys.exit(0)
if len(sys.argv) <= 3:
    print("Error: options < 3!\nUsage: {} input_dir selected_species_or_sequences(separated by comma) output_dir fasta_suffix(default='.fa') separate_symbol(default='.') if_list_for_exclusion(default=no)".format(sys.argv[0]))
    sys.exit(1)
select_list = sys.argv[2].split(',')
suffix = '.fa'
if len(sys.argv) > 4:
	suffix = sys.argv[4]
sep = '.'
if len(sys.argv) > 5:
	sep = sys.argv[5]
if_select = True
if len(sys.argv) > 6:
	if sys.argv[6] == 'no':
		print("Error: ambiguous option! If you want to select the species or sequences in the list instead of to exclude them, please not input the last option!")
		sys.exit(1)
	if_select = False


if not os.path.isdir(sys.argv[3]):
	os.mkdir(sys.argv[3])

files = os.listdir(sys.argv[1])
for filename in files:
	if not filename.endswith(suffix):
		continue
	outfile = open(os.path.join(sys.argv[3], filename), 'w')
	if_output = False
	for line in open(os.path.join(sys.argv[1], filename)):
		if line.startswith('>'):
			seqid = line.lstrip('>').rstrip().split(' ')[0]
			if if_select:
				if_output = False
			else:
				if_output = True
			for sp in select_list:
				if if_select:
					if seqid == sp or seqid.startswith(sp + sep):
						if_output = True
						break
				else:
					if seqid == sp or seqid.startswith(sp + sep):
						if_output = False
						break
		if if_output:
			outfile.write(line)
	outfile.close()

