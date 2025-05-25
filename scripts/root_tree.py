#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
from ete3 import Tree

if len(sys.argv) == 1 or sys.argv[1] == '-h':
	print("Usage: {} input.nwk output.nwk outgroup/outgroups(default=the midpoint outgroup, separated by comma)".format(sys.argv[0]))
	sys.exit(0)
if len(sys.argv) < 3:
	print("Error: options < 2!\nUsage: {} input.nwk output.nwk outgroup/outgroups(default=the midpoint outgroup, separated by comma)".format(sys.argv[0]))
	sys.exit(1)

tree = Tree(sys.argv[1])
if len(sys.argv) > 3:
	outgroup = sys.argv[3]
	if ',' in outgroup:
		outgroup = tree.get_common_ancestor(outgroup.split(','))
else:
	outgroup = tree.get_midpoint_outgroup()
tree.set_outgroup(outgroup)
tree.write(outfile=sys.argv[2])
