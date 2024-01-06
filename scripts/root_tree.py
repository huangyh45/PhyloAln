#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
from ete3 import Tree

if len(sys.argv) == 1 or sys.argv[1] == '-h':
	print("Usage: {} input.nwk output.nwk outgroup/outgroups(seperated by comma)".format(sys.argv[0]))
	sys.exit(0)
if len(sys.argv) <= 3:
	print("Error: options < 3!\nUsage: {} input.nwk output.nwk outgroup/outgroups(seperated by comma)".format(sys.argv[0]))
	sys.exit(1)

tree = Tree(sys.argv[1])
outgroup = sys.argv[3]
if ',' in outgroup:
	outgroup = tree.get_common_ancestor(outgroup.split(','))
tree.set_outgroup(outgroup)
tree.write(outfile=sys.argv[2])
