#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
from ete3 import Tree

if len(sys.argv) == 1 or sys.argv[1] == '-h':
	print("Usage: {} input.nwk output.nwk seq/seqs(separated by comma)_in_clade1_for_deletion( seq/seqs_in_clade2_for_deletion ...)".format(sys.argv[0]))
	sys.exit(0)
if len(sys.argv) < 4:
	print("Error: options < 3!\nUsage: {} input.nwk output.nwk seq/seqs(separated by comma)_in_clade1_for_deletion( seq/seqs_in_clade2_for_deletion ...)".format(sys.argv[0]))
	sys.exit(1)

tree = Tree(sys.argv[1])
leafids = []
for leaf in tree:
	leafids.append(leaf.name)
dists = {}
for seqids in sys.argv[3:]:
	if ',' in seqids:
		clade = tree.get_common_ancestor(seqids.split(','))
		for leaf in clade:
			leafids.remove(leaf.name)
	else:
		clade = tree.search_nodes(name=seqids)[0]
		leafids.remove(seqids)
	parent = clade.up
	if len(parent.children) == 2:
		sis_leaves = []
		for leaf in parent:
			if leaf not in clade and leaf.name != clade.name:
				sis_leaves.append(leaf.name)
		dists[','.join(sis_leaves)] = parent.dist + parent.children[0].dist + parent.children[1].dist - clade.dist
	#clade.detach()
tree.prune(leafids)
for seqids, dist in dists.items():
	if ',' in seqids:
		tree.get_common_ancestor(seqids.split(',')).dist = dist
	else:
		tree.search_nodes(name=seqids)[0].dist = dist 
tree.write(outfile=sys.argv[2])
