#!/usr/bin/env python3

import sys
import os
from copy import deepcopy
from Bio.Seq import reverse_complement, translate
import library as lib

# class to process the information of the read or sequence hits
class read_hit:

	# construct the aligned site-level information of the hit
	def __init__(self, info, seqids, seqstr0, site_composition, start_trim=0, end_trim=0, extend_pos=0, pos_freq=0.2, gencode=1, split_len=None, ignore_info=False):
		self.qstart = int(info[0])
		self.qend = int(info[1])
		self.rep = len(seqids)
		self.count = {}
		self.seqids = {}
		self.unmap = {}
		if split_len is None:
			for seqid in seqids:
				self.seqids.setdefault(seqid, [])
				self.seqids[seqid].append([1, len(seqstr0)])
				self.count[seqid] = 1
				self.unmap.setdefault(seqid, {})
		else:
			for seqid in seqids:
				split_start, split_end = seqid.split('_split_')[-1].split('_')
				seqid = '_split_'.join(seqid.split('_split_')[:-1])
				self.seqids.setdefault(seqid, [])
				self.seqids[seqid].append([int(split_start), int(split_end)])
				self.count[seqid] = 1
				self.unmap.setdefault(seqid, {})
		self.seq_comp = {}
		self.base = {}
		pos = 0
		j = 0
		if info[5] == 'None':
			# alignment without codon (dna or dna_codon mode)
			if info[4] == 'direct':
				if int(info[2]) > int(info[3]):
					temp = info[3]
					info[3] = info[2]
					info[2] = temp
					info[4] = '-'
				else:
					info[4] = '+'
			if info[4] == '-':
				seqstr = reverse_complement(seqstr0[(int(info[2])-1):int(info[3])])
				for seqid, SEs in self.seqids.items():
					for i in range(len(SEs)):
						self.seqids[seqid][i] = [SEs[i][0]+int(info[2])-1, SEs[i][0]+int(info[3])-1, '-']
			elif info[4] == 'rev':
				seqstr = reverse_complement(seqstr0)
				seqstr = seqstr[(int(info[2])-1):int(info[3])]
				for seqid, SEs in self.seqids.items():
					for i in range(len(SEs)):
						self.seqids[seqid][i] = [SEs[i][1]-int(info[3])+1, SEs[i][1]-int(info[2])+1, '-']
			else:
				seqstr = seqstr0[(int(info[2])-1):int(info[3])]
				for seqid, SEs in self.seqids.items():
					for i in range(len(SEs)):
						self.seqids[seqid][i] = [SEs[i][0]+int(info[2])-1, SEs[i][0]+int(info[3])-1, '+']
			for i in range(len(seqstr)):
				while info[7][pos] == '-':
					self.seq_comp[int(info[0]) + j] = {'-': self.rep}
					self.base[int(info[0]) + j] = '-'
					j += 1
					pos += 1
				if info[6][pos] == '-':
					for seqid, SEs in self.seqids.items():
						for SE in SEs:
							if SE[2] == '+':
								self.unmap[seqid][SE[0] + i] = 1
							else:
								self.unmap[seqid][SE[1] - i] = 1
				else:
					self.seq_comp[int(info[0]) + j] = {seqstr[i]: self.rep}
					self.base[int(info[0]) + j] = seqstr[i]
					j += 1
				pos += 1
			if start_trim == 0:
				i = 1
				bases = []
				while i <= extend_pos and self.qstart - i > 0:
					if info[4] == '-':
						if int(info[3]) + i > len(seqstr0):
							break
						base = reverse_complement(seqstr0[int(info[3]) + i - 1])
					else:
						if int(info[2]) - i < 1:
							break
						if info[4] == 'rev':
							base = reverse_complement(seqstr0)[int(info[2]) - i - 1]
						else:
							base = seqstr0[int(info[2]) - i - 1]
					pos = self.qstart - i
					if site_composition[pos].get(base, 0) / sum(site_composition[pos].values()) < pos_freq:
						break
					else:
						bases.append(base)
					i += 1
				if bases:
					extend = len(bases)
					for seqid, SEs in self.seqids.items():
						for i in range(len(SEs)):
							if self.seqids[seqid][i][2] == '+':
								self.seqids[seqid][i][0] -= extend
							else:
								self.seqids[seqid][i][1] += extend
					for i in range(1, extend+1):
						pos = self.qstart - i
						self.seq_comp[pos] = {bases[i-1]: self.rep}
						self.base[pos] = bases[i-1]
					self.qstart = self.qstart - extend
			if end_trim == 0:
				i = 1
				bases = []
				aln_len = max(site_composition.keys())
				while i <= extend_pos and self.qend + i <= aln_len:
					if info[4] == '-':
						if int(info[2]) - i < 1:
							break
						base = reverse_complement(seqstr0[int(info[2]) - i - 1])
					else:
						if int(info[3]) + i > len(seqstr0):
							break
						if info[4] == 'rev':
							base = reverse_complement(seqstr0)[int(info[3]) + i - 1]
						else:
							base = seqstr0[int(info[3]) + i - 1]
					pos = self.qend + i
					if site_composition[pos].get(base, 0) / sum(site_composition[pos].values()) < pos_freq:
						break
					else:
						bases.append(base)
					i += 1
				if bases:
					extend = len(bases)
					for seqid, SEs in self.seqids.items():
						for i in range(len(SEs)):
							if self.seqids[seqid][i][2] == '+':
								self.seqids[seqid][i][1] += extend
							else:
								self.seqids[seqid][i][0] -= extend
					for i in range(1, extend+1):
						pos = self.qend + i
						self.seq_comp[pos] = {bases[i-1]: self.rep}
						self.base[pos] = bases[i-1]
					self.qend = self.qend + extend
		else:
			# codon alignment
			if info[4] == 'direct':
				if int(info[2]) > int(info[3]):
					temp = info[3]
					info[3] = info[2]
					info[2] = temp
					strand = '-'
					seqstr = reverse_complement(seqstr0[(int(info[2])-1):int(info[3])])
				else:
					strand = '+'
					seqstr = seqstr0[(int(info[2])-1):int(info[3])]
				for seqid, SEs in self.seqids.items():
					for i in range(len(SEs)):
						self.seqids[seqid][i] = [SEs[i][0]+int(info[2])-1, SEs[i][0]+int(info[3])-1, strand]
			elif info[4] == '-':
				seqstr = reverse_complement(seqstr0[(3*int(info[2])+int(info[5])-3):(3*int(info[3])+int(info[5]))])
				for seqid, SEs in self.seqids.items():
					for i in range(len(SEs)):
						self.seqids[seqid][i] = [SEs[i][0]+3*int(info[2])+int(info[5])-3, SEs[i][0]+3*int(info[3])+int(info[5])-1, '-']
			elif info[4] == 'rev':
				seqstr = reverse_complement(seqstr0)
				seqstr = seqstr[(3*int(info[2])+int(info[5])-3):(3*int(info[3])+int(info[5]))]
				for seqid, SEs in self.seqids.items():
					for i in range(len(SEs)):
						self.seqids[seqid][i] = [SEs[i][1]-3*int(info[3])-int(info[5])+1, SEs[i][1]-3*int(info[2])-int(info[5])+3, '-']
			else:
				seqstr = seqstr0[(3*int(info[2])+int(info[5])-3):(3*int(info[3])+int(info[5]))]
				for seqid, SEs in self.seqids.items():
					for i in range(len(SEs)):
						self.seqids[seqid][i] = [SEs[i][0]+3*int(info[2])+int(info[5])-3, SEs[i][0]+3*int(info[3])+int(info[5])-1, '+']
			# debug
			#print(seqids, info, len(seqstr), seqstr)
			for i in range(0, len(seqstr), 3):
				while info[7][pos] == '-':
					self.seq_comp[int(info[0]) + j] = {1: {'-': self.rep}, 2: {'-': self.rep}, 3: {'-': self.rep}}
					self.base[int(info[0]) + j] = '-'
					j += 1
					pos += 1
				if info[6][pos] == '-':
					for seqid, SEs in self.seqids.items():
						for SE in SEs:
							for k in range(3):
								if SE[2] == '+':
									self.unmap[seqid][SE[0] + i + k] = 1
								else:
									self.unmap[seqid][SE[1] - i - k] = 1
				else:
					self.seq_comp[int(info[0]) + j] = {1: {seqstr[i]: self.rep}, 2: {seqstr[i+1]: self.rep}, 3: {seqstr[i+2]: self.rep}}
					self.base[int(info[0]) + j] = info[7][pos]
					j += 1
				pos += 1
			if start_trim == 0:
				i = 1
				bases = []
				codons = []
				while i <= extend_pos and self.qstart - i > 0:
					if info[4] == 'direct':
						if strand == '-':
							if int(info[3]) + 3*i > len(seqstr0):
								break
							codon = reverse_complement(seqstr0[(int(info[3])+3*i-3):(int(info[3])+3*i)])
						else:
							if int(info[2]) - 3*i < 1:
								break
							codon = seqstr0[(int(info[2])-3*i-1):(int(info[2])-3*i+2)]
					elif info[4] == '-':
						if 3*int(info[3]) + int(info[5]) + 3*i > len(seqstr0):
							break
						codon = reverse_complement(seqstr0[(3*int(info[3])+int(info[5])+3*i-3):(3*int(info[3])+int(info[5])+3*i)])
					else:
						if int(info[2]) - i < 1:
							break
						if info[4] == 'rev':
							codon = reverse_complement(seqstr0)[(3*int(info[2])+int(info[5])-3*i-3):(3*int(info[2])+int(info[5])-3*i)]
						else:
							codon = seqstr0[(3*int(info[2])+int(info[5])-3*i-3):(3*int(info[2])+int(info[5])-3*i)]
					base = translate(codon, table=gencode)
					pos = self.qstart - i
					if site_composition[pos].get(base, 0) / sum(site_composition[pos].values()) < pos_freq:
						break
					else:
						codons.append(codon)
						bases.append(base)
					i += 1
				if bases:
					extend = len(bases)
					for seqid, SEs in self.seqids.items():
						for i in range(len(SEs)):
							if self.seqids[seqid][i][2] == '+':
								self.seqids[seqid][i][0] -= 3*extend
							else:
								self.seqids[seqid][i][1] += 3*extend
					for i in range(1, extend+1):
						pos = self.qstart - i
						self.seq_comp[pos] = {1: {codons[i-1][0]: self.rep}, 2: {codons[i-1][1]: self.rep}, 3: {codons[i-1][2]: self.rep}}
						self.base[pos] = bases[i-1]
					self.qstart = self.qstart - extend
					# debug
					#print('extend start', extend, vars(self))
					#input()
			if end_trim == 0:
				i = 1
				bases = []
				codons = []
				aln_len = max(site_composition.keys())
				while i <= extend_pos and self.qend + i <= aln_len:
					if info[4] == 'direct':
						if strand == '-':
							if int(info[2]) - 3*i < 1:
								break
							codon = reverse_complement(seqstr0[(int(info[2])-3*i-1):(int(info[2])-3*i+2)])
						else:
							if int(info[3]) + 3*i > len(seqstr0):
								break
							codon = seqstr0[(int(info[3])+3*i-3):(int(info[3])+3*i)]
					elif info[4] == '-':
						if int(info[2]) - i < 1:
							break
						codon = reverse_complement(seqstr0[(3*int(info[2])+int(info[5])-3*i-3):(3*int(info[2])+int(info[5])-3*i)])
					else:
						if 3*int(info[3]) + int(info[5]) + 3*i > len(seqstr0):
							break
						if info[4] == 'rev':
							codon = reverse_complement(seqstr0)[(3*int(info[3])+int(info[5])+3*i-3):(3*int(info[3])+int(info[5])+3*i)]
						else:
							codon = seqstr0[(3*int(info[3])+int(info[5])+3*i-3):(3*int(info[3])+int(info[5])+3*i)]
					base = translate(codon, table=gencode)
					pos = self.qend + i
					if site_composition[pos].get(base, 0) / sum(site_composition[pos].values()) < pos_freq:
						break
					else:
						codons.append(codon)
						bases.append(base)
					i += 1
				if bases:
					extend = len(bases)
					for seqid, SEs in self.seqids.items():
						for i in range(len(SEs)):
							if self.seqids[seqid][i][2] == '+':
								self.seqids[seqid][i][1] += 3*extend
							else:
								self.seqids[seqid][i][0] -= 3*extend
					for i in range(1, extend+1):
						pos = self.qend + i
						self.seq_comp[pos] = {1: {codons[i-1][0]: self.rep}, 2: {codons[i-1][1]: self.rep}, 3: {codons[i-1][2]: self.rep}}
						self.base[pos] = bases[i-1]
					self.qend = self.qend + extend
					# debug
					#print('extend end', extend, vars(self))
					#input()
		# when the hits are too much, ignore the target sequence information
		if ignore_info:
			self.count = {'too much': sum(self.count.values())}
			self.seqids = {'too much': [[0, 0, '.']]}
			self.unmap = {}

	# directly merge the hits into consensus
	def consensus(self, info, seqids, seqstr0, site_composition, start_trim=0, end_trim=0, extend_pos=0, pos_freq=0.2, mol_type='dna', gencode=1, split_len=None):
		qrange = [self.qstart, self.qend]
		if int(info[0]) < self.qstart:
			self.qstart = int(info[0])
		if int(info[1]) > self.qend:
			self.qend = int(info[1])
		repnum = len(seqids)
		seqinfo = {}
		if split_len is None:
			for seqid in seqids:
				seqinfo.setdefault(seqid, [])
				seqinfo[seqid].append([1, len(seqstr0)])
		else:
			for seqid in seqids:
				split_start, split_end = seqid.split('_split_')[-1].split('_')
				seqid = '_split_'.join(seqid.split('_split_')[:-1])
				seqinfo.setdefault(seqid, [])
				seqinfo[seqid].append([int(split_start), int(split_end)])
		pos = 0
		j = 0
		if info[5] == 'None':
			# alignment without codon (dna or dna_codon mode)
			if info[4] == 'direct':
				if int(info[2]) > int(info[3]):
					temp = info[3]
					info[3] = info[2]
					info[2] = temp
					info[4] = '-'
				else:
					info[4] = '+'
			if info[4] == '-':
				seqstr = reverse_complement(seqstr0[(int(info[2])-1):int(info[3])])
				if self.count.get('too much') is None:
					for seqid, SEs in seqinfo.items():
						for i in range(len(SEs)):
							seqinfo[seqid][i] = [SEs[i][0]+int(info[2])-1, SEs[i][0]+int(info[3])-1, '-']
			elif info[4] == 'rev':
				seqstr = reverse_complement(seqstr0)
				seqstr = seqstr[(int(info[2])-1):int(info[3])]
				if self.count.get('too much') is None:
					for seqid, SEs in seqinfo.items():
						for i in range(len(SEs)):
							seqinfo[seqid][i] = [SEs[i][1]-int(info[3])+1, SEs[i][1]-int(info[2])+1, '-']
			else:
				seqstr = seqstr0[(int(info[2])-1):int(info[3])]
				if self.count.get('too much') is None:
					for seqid, SEs in seqinfo.items():
						for i in range(len(SEs)):
							seqinfo[seqid][i] = [SEs[i][0]+int(info[2])-1, SEs[i][0]+int(info[3])-1, '+']
			for i in range(len(seqstr)):
				while info[7][pos] == '-':
					self.seq_comp.setdefault(int(info[0]) + j, {})
					self.seq_comp[int(info[0]) + j].setdefault('-', 0)
					self.seq_comp[int(info[0]) + j]['-'] += repnum
					j += 1
					pos += 1
				if info[6][pos] != '-':
					self.seq_comp.setdefault(int(info[0]) + j, {})
					self.seq_comp[int(info[0]) + j].setdefault(seqstr[i], 0)
					self.seq_comp[int(info[0]) + j][seqstr[i]] += repnum
					j += 1
				pos += 1
			if start_trim == 0:
				i = 1
				bases = []
				while i <= extend_pos and int(info[0]) - i > 0:
					if info[4] == '-':
						if int(info[3]) + i > len(seqstr0):
							break
						base = reverse_complement(seqstr0[int(info[3]) + i - 1])
					else:
						if int(info[2]) - i < 1:
							break
						if info[4] == 'rev':
							base = reverse_complement(seqstr0)[int(info[2]) - i - 1]
						else:
							base = seqstr0[int(info[2]) - i - 1]
					pos = int(info[0]) - i
					if site_composition[pos].get(base, 0) / sum(site_composition[pos].values()) < pos_freq:
						break
					else:
						bases.append(base)
					i += 1
				if bases:
					extend = len(bases)
					if self.count.get('too much') is None:
						for seqid, SEs in seqinfo.items():
							for i in range(len(SEs)):
								if seqinfo[seqid][i][2] == '+':
									seqinfo[seqid][i][0] -= extend
								else:
									seqinfo[seqid][i][1] += extend
					for i in range(1, extend+1):
						pos = int(info[0]) - i
						self.seq_comp.setdefault(pos, {})
						self.seq_comp[pos].setdefault(bases[i-1], 0)
						self.seq_comp[pos][bases[i-1]] += repnum
					if int(info[0]) - extend < self.qstart:
						self.qstart = int(info[0]) - extend
			if end_trim == 0:
				i = 1
				bases = []
				aln_len = max(site_composition.keys())
				while i <= extend_pos and int(info[1]) + i <= aln_len:
					if info[4] == '-':
						if int(info[2]) - i < 1:
							break
						base = reverse_complement(seqstr0[int(info[2]) - i - 1])
					else:
						if int(info[3]) + i > len(seqstr0):
							break
						if info[4] == 'rev':
							base = reverse_complement(seqstr0)[int(info[3]) + i - 1]
						else:
							base = seqstr0[int(info[3]) + i - 1]
					pos = int(info[1]) + i
					if site_composition[pos].get(base, 0) / sum(site_composition[pos].values()) < pos_freq:
						break
					else:
						bases.append(base)
					i += 1
				if bases:
					extend = len(bases)
					if self.count.get('too much') is None:
						for seqid, SEs in seqinfo.items():
							for i in range(len(SEs)):
								if seqinfo[seqid][i][2] == '+':
									seqinfo[seqid][i][1] += extend
								else:
									seqinfo[seqid][i][0] -= extend
					for i in range(1, extend+1):
						pos = int(info[1]) + i
						self.seq_comp.setdefault(pos, {})
						self.seq_comp[pos].setdefault(bases[i-1])
						self.seq_comp[pos][bases[i-1]] += repnum
					if int(info[1]) + extend > self.qend:
						self.qend = int(info[1]) + extend
		else:
			# codon alignment
			if info[4] == 'direct':
				if int(info[2]) > int(info[3]):
					temp = info[3]
					info[3] = info[2]
					info[2] = temp
					strand = '-'
					seqstr = reverse_complement(seqstr0[(int(info[2])-1):int(info[3])])
				else:
					strand = '+'
					seqstr = seqstr0[(int(info[2])-1):int(info[3])]
				if self.count.get('too much') is None:
					for seqid, SEs in seqinfo.items():
						for i in range(len(SEs)):
							seqinfo[seqid][i] = [SEs[i][0]+int(info[2])-1, SEs[i][0]+int(info[3])-1, strand]
			elif info[4] == '-':
				seqstr = reverse_complement(seqstr0[(3*int(info[2])+int(info[5])-3):(3*int(info[3])+int(info[5]))])
				if self.count.get('too much') is None:
					for seqid, SEs in seqinfo.items():
						for i in range(len(SEs)):
							seqinfo[seqid][i] = [SEs[i][0]+3*int(info[2])+int(info[5])-3, SEs[i][0]+3*int(info[3])+int(info[5])-1, '-']
			elif info[4] == 'rev':
				seqstr = reverse_complement(seqstr0)
				seqstr = seqstr[(3*int(info[2])+int(info[5])-3):(3*int(info[3])+int(info[5]))]
				if self.count.get('too much') is None:
					for seqid, SEs in seqinfo.items():
						for i in range(len(SEs)):
							seqinfo[seqid][i] = [SEs[i][1]-3*int(info[3])-int(info[5])+1, SEs[i][1]-3*int(info[2])-int(info[5])+3, '-']
			else:
				seqstr = seqstr0[(3*int(info[2])+int(info[5])-3):(3*int(info[3])+int(info[5]))]
				if self.count.get('too much') is None:
					for seqid, SEs in seqinfo.items():
						for i in range(len(SEs)):
							seqinfo[seqid][i] = [SEs[i][0]+3*int(info[2])+int(info[5])-3, SEs[i][0]+3*int(info[3])+int(info[5])-1, '+']
			for i in range(0, len(seqstr), 3):
				while info[7][pos] == '-':
					self.seq_comp.setdefault(int(info[0]) + j, {1: {}, 2: {}, 3: {}})
					for k in [1, 2, 3]:
						self.seq_comp[int(info[0]) + j][k].setdefault('-', 0)
						self.seq_comp[int(info[0]) + j][k]['-'] += repnum
					j += 1
					pos += 1
				if info[6][pos] != '-':
					self.seq_comp.setdefault(int(info[0]) + j, {1: {}, 2: {}, 3: {}})
					for k in [1, 2, 3]:
						self.seq_comp[int(info[0]) + j][k].setdefault(seqstr[i+k-1], 0)
						self.seq_comp[int(info[0]) + j][k][seqstr[i+k-1]] += repnum
					j += 1
				pos += 1
			if start_trim == 0:
				i = 1
				bases = []
				codons = []
				while i <= extend_pos and int(info[0]) - i > 0:
					if info[4] == 'direct':
						if strand == '-':
							if int(info[3]) + 3*i > len(seqstr0):
								break
							codon = reverse_complement(seqstr0[(int(info[3])+3*i-3):(int(info[3])+3*i)])
						else:
							if int(info[2]) - 3*i < 1:
								break
							codon = seqstr0[(int(info[2])-3*i-1):(int(info[2])-3*i+2)]
					elif info[4] == '-':
						if 3*int(info[3]) + int(info[5]) + 3*i > len(seqstr0):
							break
						codon = reverse_complement(seqstr0[(3*int(info[3])+int(info[5])+3*i-3):(3*int(info[3])+int(info[5])+3*i)])
					else:
						if int(info[2]) - i < 1:
							break
						if info[4] == 'rev':
							codon = reverse_complement(seqstr0)[(3*int(info[2])+int(info[5])-3*i-3):(3*int(info[2])+int(info[5])-3*i)]
						else:
							codon = seqstr0[(3*int(info[2])+int(info[5])-3*i-3):(3*int(info[2])+int(info[5])-3*i)]
					base = translate(codon, table=gencode)
					pos = int(info[0]) - i
					if site_composition[pos].get(base, 0) / sum(site_composition[pos].values()) < pos_freq:
						break
					else:
						codons.append(codon)
						bases.append(base)
					i += 1
				if bases:
					extend = len(bases)
					if self.count.get('too much') is None:
						for seqid, SEs in seqinfo.items():
							for i in range(len(SEs)):
								if seqinfo[seqid][i][2] == '+':
									seqinfo[seqid][i][0] -= 3*extend
								else:
									seqinfo[seqid][i][1] += 3*extend
					for i in range(1, extend+1):
						pos = int(info[0]) - i
						self.seq_comp.setdefault(pos, {1: {}, 2: {}, 3: {}})
						for k in [1, 2, 3]:
							self.seq_comp[pos][k].setdefault(codons[i-1][k-1], 0)
							self.seq_comp[pos][k][codons[i-1][k-1]] += repnum
					if int(info[0]) - extend < self.qstart:
						self.qstart = int(info[0]) - extend
					# debug
					#print('extend start', extend, vars(self))
					#input()
			if end_trim == 0:
				i = 1
				bases = []
				codons = []
				aln_len = max(site_composition.keys())
				while i <= extend_pos and int(info[1]) + i <= aln_len:
					if info[4] == 'direct':
						if strand == '-':
							if int(info[2]) - 3*i < 1:
								break
							codon = reverse_complement(seqstr0[(int(info[2])-3*i-1):(int(info[2])-3*i+2)])
						else:
							if int(info[3]) + 3*i > len(seqstr0):
								break
							codon = seqstr0[(int(info[3])+3*i-3):(int(info[3])+3*i)]
					elif info[4] == '-':
						if int(info[2]) - i < 1:
							break
						codon = reverse_complement(seqstr0[(3*int(info[2])+int(info[5])-3*i-3):(3*int(info[2])+int(info[5])-3*i)])
					else:
						if 3*int(info[3]) + int(info[5]) + 3*i > len(seqstr0):
							break
						if info[4] == 'rev':
							codon = reverse_complement(seqstr0)[(3*int(info[3])+int(info[5])+3*i-3):(3*int(info[3])+int(info[5])+3*i)]
						else:
							codon = seqstr0[(3*int(info[3])+int(info[5])+3*i-3):(3*int(info[3])+int(info[5])+3*i)]
					base = translate(codon, table=gencode)
					pos = int(info[1]) + i
					if site_composition[pos].get(base, 0) / sum(site_composition[pos].values()) < pos_freq:
						break
					else:
						codons.append(codon)
						bases.append(base)
					i += 1
				if bases:
					extend = len(bases)
					if self.count.get('too much') is None:
						for seqid, SEs in seqinfo.items():
							for i in range(len(SEs)):
								if seqinfo[seqid][i][2] == '+':
									seqinfo[seqid][i][1] += 3*extend
								else:
									seqinfo[seqid][i][0] -= 3*extend
					for i in range(1, extend+1):
						pos = int(info[1]) + i
						self.seq_comp.setdefault(pos, {1: {}, 2: {}, 3: {}})
						for k in [1, 2, 3]:
							self.seq_comp[pos][k].setdefault(codons[i-1][k-1], 0)
							self.seq_comp[pos][k][codons[i-1][k-1]] += repnum
					if int(info[1]) + extend > self.qend:
						self.qend = int(info[1]) + extend
					# debug
					#print('extend end', extend, vars(self))
					#input()
		if self.count.get('too much') is None:
			for seqid, SEs in seqinfo.items():
				self.seqids.setdefault(seqid, [])
				if split_len:
					# merge the splitted fragments if possible
					if qrange[0] < int(info[0]) or (qrange[0] == int(info[0]) and qrange[1] < int(info[1])):
						ifrev = False
					else:
						ifrev = True
					for SE in SEs:
						ifmerge = False
						for j in range(len(self.seqids[seqid])):
							if self.seqids[seqid][j][2] != SE[2] or (mol_type != 'dna' and (SE[0] - self.seqids[seqid][j][0]) % 3 != 0):
								continue
							elif (SE[2] == '+' and not ifrev) or (SE[2] == '-' and ifrev):
								if SE[0] >= self.seqids[seqid][j][0] and SE[0] <= self.seqids[seqid][j][1] + 1 and SE[1] >= self.seqids[seqid][j][1]:
									self.seqids[seqid][j][1] = SE[1]
									ifmerge = True
									break
							else:
								if SE[1] <= self.seqids[seqid][j][1] and SE[1] + 1 >= self.seqids[seqid][j][0] and SE[0] <= self.seqids[seqid][j][0]:
									self.seqids[seqid][j][0] = SE[0]
									ifmerge = True
									break
						if not ifmerge:
							self.seqids[seqid].append(deepcopy(SE))
				else:
					self.seqids[seqid].extend(deepcopy(SEs))
				self.count[seqid] = 1
		else:
			self.count['too much'] += len(seqinfo)

	# map from reference position to target position
	def map_pos(self, seqid, mol_type='dna'):
		self.pos = {}
		self.qpos = {}
		if self.seqids[seqid][0][2] == '-':
			step = -1
			i = self.seqids[seqid][0][1]
		else:
			step = 1
			i = self.seqids[seqid][0][0]
		for j, info in self.seq_comp.items():
			if mol_type == 'dna':
				self.pos[j] = None
				if self.seq_comp[j].get('-') is None:
					self.pos[j] = i
					self.qpos[i] = j
					i += step
			else:
				self.pos[j] = {1: None, 2: None, 3: None}
				for k in [1, 2, 3]:
					if self.seq_comp[j][k].get('-') is None:
						self.pos[j][k] = i
						self.qpos[i] = j
						i += step
			while self.unmap[seqid].get(i) is not None:
				i += step

	# calculate pseudo RPKM of the hit using read counts
	def addRPKM(self, total_count):
		self.pseudoRPKM = 1000 * 1000000 * sum(self.count.values()) / total_count / (len(self.seq_comp) - sum(1 for value in self.base.values() if value == '-'))

	# print the sequence from the hit information
	def print_seq(self, seq, protseq=None):
		for pos, base in self.seq_comp.items():
			if protseq:
				for i in [1, 2, 3]:
					seq[3*pos+i-4] = list(base[i].keys())[0]
				protseq[pos-1] = self.base[pos]
			else:
				seq[pos-1] = self.base[pos]
		return seq, protseq

# calculate the site composition of the ingroup alignments and the conservative score of the outgroup sequence(s)
def get_ref_score(refseqs, outgroup=[], ingroup=[], sep='.', unknow='N'):
	site_composition = {}
	outgroup_score = {}
	outgroup_aln = {}
	ingroup_aln = {}
	valid_range = {}
	if outgroup:
		for spid in outgroup:
			for seqid in refseqs.keys():
				if seqid == spid or seqid.startswith(spid + sep):
					outgroup_aln[seqid] = refseqs[seqid]
					break
	else:
		outgroup_aln = deepcopy(refseqs)
	if ingroup:
		for spid in ingroup:
			for seqid in refseqs.keys():
				if seqid == spid or seqid.startswith(spid + sep):
					ingroup_aln[seqid] = refseqs[seqid]
					break
	else:
		if len(outgroup_aln) == len(refseqs):
			ingroup_aln = deepcopy(refseqs)
		else:
			for seqid, seqstr in refseqs.items():
				if outgroup_aln.get(seqid) is None:
					ingroup_aln[seqid] = seqstr
	gap_symbols = ['-', unknow]
	for seqid, seqstr in refseqs.items():
		valid_range[seqid] = []
		gap_num = 0
		no_gap_num = 0
		for i in range(len(seqstr)):
			if seqstr[i] in gap_symbols:
				gap_num += 1
				if gap_num > 5:
					if no_gap_num > 5:
						valid_range[seqid].append([i - gap_num - no_gap_num + 1, i - gap_num])
					no_gap_num = 0
			else:
				if no_gap_num > 0 and gap_num > 0 and gap_num <= 5:
					no_gap_num += gap_num
				no_gap_num += 1
				gap_num = 0
		if no_gap_num > 5:
			valid_range[seqid].append([i - gap_num - no_gap_num + 1, i - gap_num])
		# debug
		#print(seqid, valid_range[seqid])
		#input()
	out_scores = {}
	for seqid in outgroup_aln.keys():
		outgroup_score[seqid] = {}
		out_scores[seqid] = []
	for i in range(len(list(refseqs.values())[0])):
		site_composition[i+1] = {}
		ifvalid = {}
		for seqid, seqstr in ingroup_aln.items():
			ifvalid[seqid] = False
			for vrange in valid_range[seqid]:
				if i >= vrange[0] and i <= vrange[1]:
					ifvalid[seqid] = True
					break
			if not ifvalid[seqid]:
				continue
			base = seqstr[i]
			if site_composition[i+1].get(base) is None:
				site_composition[i+1][base] = 1
			else:
				site_composition[i+1][base] += 1
		for seqid, seqstr in outgroup_aln.items():
			if ifvalid.get(seqid) is None:
				ifvalid[seqid] = False
				for vrange in valid_range[seqid]:
					if i >= vrange[0] and i <= vrange[1]:
						ifvalid[seqid] = True
						break
			if ifvalid[seqid]:
				outgroup_score[seqid][i+1] = site_composition[i+1].get(seqstr[i], 0)
				out_scores[seqid].append(outgroup_score[seqid][i+1] / sum(site_composition[i+1].values()))
	# calculate the average frequency of the outgroups as an alternative outgroup score if all outgroup sequences are unavaliable
	for seqid in outgroup_aln.keys():
		out_scores[seqid] = sum(out_scores[seqid]) / len(out_scores[seqid])
	aver_freq = min(out_scores.values())
	outgroup_score['PhyloAln_Alt'] = {}
	for i in range(len(list(refseqs.values())[0])):
		outgroup_score['PhyloAln_Alt'][i+1] = sum(site_composition[i+1].values()) * aver_freq
	return site_composition, outgroup_score

# quickly remove the foreign hits not to be assembled based on conservative scores and trim the terminal positions
def pretreat_hits(info, site_composition, outgroup_score, no_assemble=False, outgroup_weight=1, trim_pos=5, pos_freq=0.2):
	# remove the accidental short hits generated by HMMER
	if len(info[9]) <= 3:
		return None, 0, 0
	if outgroup_score and no_assemble:
		score = 0
		out_score = {}
		for seqid, seqscore in outgroup_score.items():
			if seqscore.get(int(info[2])) is None or seqscore.get(int(info[3])) is None:
				continue
			out_score[seqid] = 0
			for pos in range(int(info[2])+1, int(info[3])):
				if seqscore.get(pos) is None:
					del out_score[seqid]
					break
		# use the alternative outgroup score if all outgroups are unavailable
		if not out_score:
			out_score['PhyloAln_Alt'] = 0
		i = 0
		j = 0
		for i in range(len(info[9])):
			if info[8][i] != '-':
				pos = int(info[2]) + j
				score += site_composition[pos].get(info[9][i], 0)
				j += 1
				for seqid in out_score.keys():
					out_score[seqid] += outgroup_score[seqid][pos]
		# debug
		#print(info, score, out_score)
		#input()
		if score < min(out_score.values()) * outgroup_weight:
			return None, 0, 0

	start_trim = 0
	i = 0
	j = 0
	while start_trim < trim_pos and i < len(info[9]):
		if info[8][i] != '-':
			pos = int(info[2]) + start_trim
			if site_composition[pos].get(info[9][i], 0) / sum(site_composition[pos].values()) < pos_freq:
				start_trim += 1
			else:
				break
		if info[9][i] != '-':
			j += 1
		i += 1
	if i == len(info[9]):
		return None, 0, 0
	# continue to trim when next target position is a gap
	while info[9][i] == '-':
		start_trim += 1
		i += 1
	info[2] = int(info[2]) + start_trim
	info[8] = info[8][i:]
	info[9] = info[9][i:]
	if info[6] == 'direct':
		if info[7] == 'None':
			step = 1
		else:
			step = 3
		if int(info[4]) > int(info[5]):
			info[4] = int(info[4]) - (j * step)
		else:
			info[4] = int(info[4]) + (j * step)
	elif info[6] == '-':
		info[5] = int(info[5]) - j
	else:
		info[4] = int(info[4]) + j
	end_trim = 0
	i = len(info[9]) - 1
	j = 0
	while end_trim < trim_pos and i >= 0:
		if info[8][i] != '-':
			pos = int(info[3]) - end_trim
			if site_composition[pos].get(info[9][i], 0) / sum(site_composition[pos].values()) < pos_freq:
				end_trim += 1
			else:
				break
		if info[9][i] != '-':
			j += 1
		i -= 1
	if i < 0:
		return None, 0, 0
	# continue to trim when next target position is a gap
	while info[9][i] == '-':
		end_trim += 1
		i -= 1
	info[3] = int(info[3]) - end_trim
	info[8] = info[8][:(i+1)]
	info[9] = info[9][:(i+1)]
	if info[6] == 'direct':
		if info[7] == 'None':
			step = 1
		else:
			step = 3
		if int(info[4]) > int(info[5]):
			info[5] = int(info[5]) + (j * step)
		else:
			info[5] = int(info[5]) - (j * step)
	elif info[6] == '-':
		info[4] = int(info[4]) + j
	else:
		info[5] = int(info[5]) - j
	# debug
	#print('trim', start_trim, end_trim, info)
	#input()
	return info, start_trim, end_trim

# calculate the most common base of the site
def most_common(bases):
	return max(bases, key=bases.get)

# calculate the identical site number between two hits
def calculate_ident(bases1, bases2):
	total_bases1 = sum(bases1.values())
	total_bases2 = sum(bases2.values())
	ident = 0
	for base in set(bases1.keys()) & set(bases2.keys()):
		ident += bases1[base] / total_bases1 * bases2[base] / total_bases2
	return ident

# assemble or merge two hits
def merge_hits(hit1, hit2, mol_type='dna', split_len=None, fill_gap=False):
	newhit = deepcopy(hit1)
	if hit2.qend > newhit.qend:
		newhit.qend = hit2.qend
	if newhit.count.get('too much') is None:
		for seqid, SEs in hit2.seqids.items():
			newhit.seqids.setdefault(seqid, [])
			if split_len:
				# merge the splitted fragments if possible
				for SE in SEs:
					ifmerge = False
					for j in range(len(newhit.seqids[seqid])):
						if newhit.seqids[seqid][j][2] != SE[2] or (mol_type != 'dna' and (SE[0] - newhit.seqids[seqid][j][0]) % 3 != 0):
							continue
						elif SE[2] == '+':
							if SE[0] >= newhit.seqids[seqid][j][0] and SE[0] <= newhit.seqids[seqid][j][1] + 1 and SE[1] >= newhit.seqids[seqid][j][1]:
								newhit.seqids[seqid][j][1] = SE[1]
								ifmerge = True
								break
						else:
							if SE[1] <= newhit.seqids[seqid][j][1] and SE[1] + 1 >= newhit.seqids[seqid][j][0] and SE[0] <= newhit.seqids[seqid][j][0]:
								newhit.seqids[seqid][j][0] = SE[0]
								ifmerge = True
								break
					if not ifmerge:
						newhit.seqids[seqid].append(deepcopy(SE))
			else:
				newhit.seqids[seqid].extend(deepcopy(SEs))
			newhit.count[seqid] = 1
	else:
		newhit.count['too much'] += hit2.count['too much']
	for j, info in hit2.seq_comp.items():
		if mol_type == 'dna':
			newhit.seq_comp.setdefault(j, {})
			for base, count in info.items():
				newhit.seq_comp[j].setdefault(base, 0)
				newhit.seq_comp[j][base] += count
		else:
			newhit.seq_comp.setdefault(j, {1: {}, 2: {}, 3:{}})
			for k in [1, 2, 3]:
				for base, count in info[k].items():
					newhit.seq_comp[j][k].setdefault(base, 0)
					newhit.seq_comp[j][k][base] += count
	if fill_gap and hit2.qstart - hit1.qend > 1:
		if mol_type == 'dna':
			for i in range(hit1.qend+1, hit2.qstart):
				newhit.seq_comp[i] = {fill_gap: 1}
				#newhit.base[i] = '-'
		else:
			for i in range(hit1.qend+1, hit2.qstart):
				newhit.seq_comp[i] = {1: {fill_gap: 1}, 2: {fill_gap: 1}, 3:{fill_gap: 1}}
				#newhit.base[i] = '-'
		if newhit.count.get('too much') is None:
			if newhit.seqids[seqid][-1][2] == '+':
				if newhit.seqids[seqid][-1][0] - newhit.seqids[seqid][-2][1] <= 1:
					newhit.seqids[seqid][-2][1] = newhit.seqids[seqid][-1][1]
					newhit.seqids[seqid].pop()
			elif newhit.seqids[seqid][-2][0] - newhit.seqids[seqid][-1][1] <= 1:
				newhit.seqids[seqid][-2][0] = newhit.seqids[seqid][-1][0]
				newhit.seqids[seqid].pop()
	# debug
	#print(vars(hit1))
	#print(vars(hit2))
	#print(vars(newhit))
	#input()
	return newhit

""" # graph-based assembly, but the results seem not better

# assemble the unique paths and simplify the hits and paths
def sub_assemble(hit_link, hits, mol_type='dna', split_len=None):
	# debug
	#print(hit_link)
	used = {}
	paths = []
	path = []
	start_hits = []
	for i, jdict in sorted(hit_link.items()):
		js = list(jdict.keys())
		if used.get(i) is not None or len(js) == 0:
			continue
		if js[-1] < i:
			continue
		if js[0] > i:
			start_hits.append(i)
		if len(js) > 1:
			if js[-2] > i:
				continue
		path.append(i)
		j = js[-1]
		while len(hit_link[j]) == 2:
			jnext = list(hit_link[j].keys())
			if jnext[0] > j or jnext[1] < j:
				break
			path.append(j)
			j = jnext[1]
			if len(hit_link[j]) == 1 or (len(hit_link[j]) > 2 and list(hit_link[j].keys())[1] > j):
				path.append(j)
				break
		if len(path) > 2:
			for j in path:
				used[j] = 1
			paths.append(path.copy())
		path = []
	for path in paths:
		for i in path[1:]:
			hits[path[0]] = merge_hits(hits[path[0]], hits[i], mol_type, split_len)
		del hit_link[path[-1]][path[-2]]
		hit_link[path[0]] = deepcopy(hit_link[path[-1]])
	for path in paths:
		for i in path[1:]:
			hit_link[i] = {}
			hits[i] = None
	# debug
	#print(hit_link)
	#print(paths)
	#input()
	return hit_link, hits, start_hits

# get all the paths to assemble the hits through recursion
def get_path(hit_link, start_num, path=[], paths=[]):
	if not path and not paths:
		path.append(start_num)
	if hit_link.get(path[-1]):
		firsti = list(hit_link[path[-1]].keys())[0]
		for i in hit_link[path[-1]].keys():
			# avoid long redundant path
			ifrep = False
			if i > firsti:
				for j in range(len(paths)-1, -1, -1):
					if i in paths[j] and path[-1] in paths[j]:
						ifmid = False
						for k in hit_link[path[-1]].keys():
							if k < i and k in paths[j]:
								ifmid = True
								break
						if ifmid:
							ifrep = True
							break
			if not ifrep:
				path.append(i)
				paths = get_path(hit_link, start_num, path, paths)
				path.pop()
	else:
		paths.append(path)
		# debug
		#print('get_path', paths)
		#input()
	return deepcopy(paths)

# conduct reference-based assembly
def assemble_on_ref(hits, overlap_len=30, overlap_pident=98, mol_type='dna', split_len=None):
	# construct the valid links between the hits
	if mol_type != 'dna':
		overlap_len = overlap_len / 3
	hits.sort(key = lambda x : x.qend)
	hits.sort(key = lambda x : x.qstart)
	hit_link = {}
	for i in range(len(hits) - 1):
		if hits[i] is None:
			continue
		hit_link.setdefault(i, {})
		j = i + 1
		while j < len(hits) and hits[j] is not None and hits[i].qend - hits[j].qstart + 1 >= overlap_len:
			end_pos = min(hits[i].qend, hits[j].qend)
			# debug
			#print(i, j, hits[i].qstart, hits[i].qend, hits[j].qstart, hits[j].qend, end_pos - hits[j].qstart + 1)
			overlap = end_pos - hits[j].qstart + 1
			if overlap < overlap_len:
				j += 1
				continue
			nident = 0
			if mol_type == 'dna':
				for k in range(hits[j].qstart, end_pos+1):
					nident += calculate_ident(hits[i].seq_comp[k], hits[j].seq_comp[k])
			else:
				for k in range(hits[j].qstart, end_pos+1):
					for pos in [1, 2, 3]:
						nident += calculate_ident(hits[i].seq_comp[k][pos], hits[j].seq_comp[k][pos])
				nident = nident / 3
			pident = nident / overlap * 100
			# debug
			#print(pident)
			if pident >= overlap_pident:
				if pident == 100 and hits[i].qstart <= hits[j].qstart and hits[i].qend >= hits[j].qend:
					hits[i] = merge_hits(hits[i], hits[j], mol_type, split_len)
					hits[j] = None
					for num in hit_link.get(j, {}).keys():
						del hit_link[num][j]
					hit_link[j] = {}
				else:
					# remove the suboptimal links
					suboptimal = False
					for num in hit_link[i].keys():
						if num > i and (hit_link[i][num][0] > pident or (hit_link[i][num][0] == pident and hit_link[i][num][1] > overlap)):
							suboptimal = True
							break
					hit_link.setdefault(j, {})
					if not suboptimal:
						for num in hit_link[j].keys():
							if num < j and (hit_link[num][j][0] > pident or (hit_link[num][j][0] == pident and hit_link[num][j][1] > overlap)):
								suboptimal = True
								break
					if not suboptimal:
						#count1 = sum(hits[i].count.values())
						#count2 = sum(hits[j].count.values())
						#hit_link[i][j] = [pident, overlap, count1, count2]
						for num in deepcopy(hit_link[i]).keys():
							if num > i and (hit_link[i][num][0] < pident or (hit_link[i][num][0] == pident and hit_link[i][num][1] < overlap)):
								del hit_link[i][num]
								del hit_link[num][i]
						for num in deepcopy(hit_link[j]).keys():
							if num < j and (hit_link[num][j][0] < pident or (hit_link[num][j][0] == pident and hit_link[num][j][1] < overlap)):
								del hit_link[j][num]
								del hit_link[num][j]
						hit_link[i][j] = [pident, overlap]
						hit_link[j][i] = 1
						# remove the short redundant path
						for num in list(deepcopy(hit_link)[j].keys())[:-1]:
							if i > num and i in hit_link.get(num, {}).keys():
								del hit_link[num][j]
								del hit_link[j][num]
			j += 1

	# assemble the unique paths and simplify the hits and paths
	hit_link, hits, start_hits = sub_assemble(hit_link, hits, mol_type, split_len)
	# record the hits without start position and then remove the reversed links (j-to-i links above)
	for i, js in deepcopy(hit_link).items():
		for j in js.keys():
			if j < i:
				del hit_link[i][j]
	# debug
	#print(hit_link, start_hits)
	#input()

	# construct the path based on the links
	paths = []
	for i in start_hits:
		paths.extend(get_path(hit_link, i, [], []))
		# debug
		#print(paths)
		#input()
	# add the hits without links if exists
	used = {}
	for i in range(len(hits)):
		if hits[i] is None:
			continue
		used[i] = 0
		for path in paths:
			if i in path:
				used[i] += 1
		if used[i] == 0:
			paths.append([i])
			used[i] = 1
	if len(hits) > 0 and hits[0].count.get('too much') is not None:
		for i in range(len(hits)):
			if hits[i] is None:
				continue
			hits[i].count['too much'] = hits[i].count['too much'] / used[i]
	# debug
	#print(hit_link, start_hits, paths, used)
	#input()

	# assemble the hits into contigs based on the paths
	contigs = []
	for path in paths:
		contig = deepcopy(hits[path[0]])
		for i in path[1:]:
			contig = merge_hits(contig, hits[i], mol_type, split_len)
		contigs.append(contig)
	contigs.sort(key = lambda x : x.qend)
	contigs.sort(key = lambda x : x.qstart)
	return contigs
"""

# conduct reference-based assembly
def assemble_on_ref(hits, overlap_len=30, overlap_pident=98, mol_type='dna', split_len=None):
	# construct the valid links between the hits
	if mol_type != 'dna':
		overlap_len = overlap_len / 3
	hits.sort(key = lambda x : x.qend)
	hits.sort(key = lambda x : x.qstart)
	hit_link = []
	for i in range(len(hits) - 1):
		if hits[i] is None:
			continue
		j = i + 1
		while j < len(hits) and hits[j] is not None and hits[i].qend - hits[j].qstart + 1 >= overlap_len:
			end_pos = min(hits[i].qend, hits[j].qend)
			# debug
			#print(i, j, hits[i].qstart, hits[i].qend, hits[j].qstart, hits[j].qend, end_pos - hits[j].qstart + 1)
			overlap = end_pos - hits[j].qstart + 1
			if overlap < overlap_len:
				j += 1
				continue
			nident = 0
			if mol_type == 'dna':
				for k in range(hits[j].qstart, end_pos+1):
					nident += calculate_ident(hits[i].seq_comp[k], hits[j].seq_comp[k])
			else:
				for k in range(hits[j].qstart, end_pos+1):
					for pos in [1, 2, 3]:
						nident += calculate_ident(hits[i].seq_comp[k][pos], hits[j].seq_comp[k][pos])
				nident = nident / 3
			pident = nident / overlap * 100
			# debug
			#print(pident)
			if pident >= overlap_pident:
				if pident == 100 and hits[i].qend >= hits[j].qend:
					hits[i] = merge_hits(hits[i], hits[j], mol_type, split_len)
					hits[j] = None
				else:
					hit_link.append([i, j, pident, nident])
			j += 1
	to_remove = []
	for i in range(len(hit_link)):
		if hits[hit_link[i][0]] is None or hits[hit_link[i][1]] is None:
			to_remove.append(i)
	for i in sorted(to_remove, reverse=True):
		hit_link.pop(i)
	while hit_link:
		hit_link.sort(key = lambda x : x[3], reverse=True)
		hit_link.sort(key = lambda x : x[2], reverse=True)
		i = hit_link[0][0]
		j = hit_link[0][1]
		hits[i] = merge_hits(hits[i], hits[j], mol_type, split_len)
		hits[j] = None
		to_remove = []
		for k in range(len(hit_link)):
			if hit_link[k][0] in [i, j] or hit_link[k][1] in [i, j]:
				to_remove.append(k)
		for k in sorted(to_remove, reverse=True):
			hit_link.pop(k)
		for k in range(i):
			if hits[k] is None:
				continue
			end_pos = min(hits[i].qend, hits[k].qend)
			overlap = end_pos - hits[i].qstart + 1
			if overlap < overlap_len:
				continue
			nident = 0
			if mol_type == 'dna':
				for n in range(hits[i].qstart, end_pos+1):
					nident += calculate_ident(hits[i].seq_comp[n], hits[k].seq_comp[n])
			else:
				for n in range(hits[i].qstart, end_pos+1):
					for pos in [1, 2, 3]:
						nident += calculate_ident(hits[i].seq_comp[n][pos], hits[k].seq_comp[n][pos])
				nident = nident / 3
			pident = nident / overlap * 100
			# debug
			#print(pident)
			if pident >= overlap_pident:
				hit_link.append([k, i, pident, nident])
		for k in range(i+1, len(hits)):
			if hits[k] is None:
				continue
			end_pos = min(hits[i].qend, hits[k].qend)
			overlap = end_pos - hits[k].qstart + 1
			if overlap < overlap_len:
				continue
			nident = 0
			if mol_type == 'dna':
				for n in range(hits[k].qstart, end_pos+1):
					nident += calculate_ident(hits[i].seq_comp[n], hits[k].seq_comp[n])
			else:
				for n in range(hits[k].qstart, end_pos+1):
					for pos in [1, 2, 3]:
						nident += calculate_ident(hits[i].seq_comp[n][pos], hits[k].seq_comp[n][pos])
				nident = nident / 3
			pident = nident / overlap * 100
			# debug
			#print(pident)
			if pident >= overlap_pident:
				hit_link.append([i, k, pident, nident])
	contigs = []
	for hit in hits:
		if hit is not None:
			contigs.append(deepcopy(hit))
	contigs.sort(key = lambda x : x.qend)
	contigs.sort(key = lambda x : x.qstart)
	return contigs

# convert the assembled contig information into simple sequence information (the most common bases instead of SNPs)
def simplify_contigs(contigs, mol_type='dna', gencode=1, strict=False):
	for i in range(len(contigs)):
		for j in contigs[i].seq_comp.keys():
			bases = contigs[i].seq_comp[j]
			if mol_type == 'dna':
				if strict and len(bases) > 1:
					contigs[i].seq_comp[j] = {'n': 1}
					contigs[i].base[j] = 'n'
				else:
					base = most_common(bases)
					contigs[i].seq_comp[j] = {base: 1}
					contigs[i].base[j] = base
			else:
				if strict and len(bases[1]) > 1:
					base1 = 'n'
				else:
					base1 = most_common(bases[1])
				if strict and len(bases[2]) > 1:
					base2 = 'n'
				else:
					base2 = most_common(bases[2])
				if strict and len(bases[3]) > 1:
					base3 = 'n'
				else:
					base3 = most_common(bases[3])
				if base1 == '-' and base2 == '-' and base3 == '-':
					base = '-'
				else:
					try:
						base = translate(base1 + base2 + base3, table=gencode)
					except:
						base = 'X'
				contigs[i].seq_comp[j]= {1: {base1: 1}, 2: {base2: 1}, 3: {base3: 1}}
				contigs[i].base[j] = base
	return contigs

# remove the foreign sequences based on conservative scores
def remove_out_seq(contigs, site_composition, outgroup_score, outgroup_weight=1):
	to_remove = []
	last_out_score = [{}]
	for i in range(len(contigs)):
		contig = contigs[i]
		start_pos = min(contig.seq_comp.keys())
		end_pos = max(contig.seq_comp.keys())
		if last_out_score[0]:
			last_overlap = min(last_out_score[2], end_pos) - start_pos + 1
		else:
			last_overlap = 0
		score = 0
		out_score = {}
		for seqid, seqscore in outgroup_score.items():
			if seqscore.get(end_pos) is None:
				continue
			out_score[seqid] = 0
			for pos in contig.seq_comp.keys():
				if seqscore.get(pos) is None:
					del out_score[seqid]
					break
		# use the alternative outgroup score if all outgroups are unavailable
		if not out_score:
			out_score['PhyloAln_Alt'] = 0
		for pos, base in contig.base.items():
			score += site_composition[pos].get(base, 0)
		for seqid in out_score.keys():
			if last_overlap / (end_pos - start_pos + 1) < 0.6 or last_out_score[0].get(seqid) is None or last_overlap / (last_out_score[2] - last_out_score[1] + 1) < 0.6:
				for pos in contig.seq_comp.keys():
					out_score[seqid] += outgroup_score[seqid][pos]
			else:
				# use the existing scores of the calculated regions to save time
				out_score[seqid] = last_out_score[0][seqid]
				for pos in range(last_out_score[1], start_pos):
					out_score[seqid] -= outgroup_score[seqid][pos]
				for pos in range(end_pos + 1, last_out_score[2] + 1):
					out_score[seqid] -= outgroup_score[seqid][pos]
				for pos in range(last_out_score[2] + 1, end_pos + 1):
					out_score[seqid] += outgroup_score[seqid][pos]
		# debug
		#print(i, contig.seqids, score, out_score)
		#input()
		if score < min(out_score.values()) * outgroup_weight:
			to_remove.append(i)
		last_out_score = [out_score.copy(), start_pos, end_pos]
	for i in reversed(to_remove):
		contigs.pop(i)
	return contigs

# merged the splitted hits before predicting genes
def merge_split(hits, mol_type='dna'):
	contigs = []
	used = {}
	for i in range(len(hits)):
		if used.get(i) is not None:
			continue
		used[i] = 1
		seqid = list(hits[i].seqids.keys())[0]
		strand = hits[i].seqids[seqid][0][2]
		contig = deepcopy(hits[i])
		for j in range(i+1, len(hits)):
			if used.get(j) is not None or hits[j].seqids.get(seqid) is None:
				continue
			if hits[j].seqids[seqid][0][2] != strand:
				continue
			if contig.pos.get(hits[j].qstart) is None or hits[j].pos.get(contig.qend) is None:
				continue
			if mol_type != 'dna' and (hits[j].seqids[seqid][0][0] - contig.seqids[seqid][0][0]) % 3 != 0:
				continue
			if contig.pos[hits[j].qstart] == hits[j].pos[hits[j].qstart] and contig.pos[contig.qend] == hits[j].pos[contig.qend]:
				# debug
				#print(vars(contig))
				#print(vars(hits[j]))
				if hits[j].qend > contig.qend:
					for pos in range(contig.qend+1, hits[j].qend+1):
						contig.seq_comp[pos] = deepcopy(hits[j].seq_comp[pos])
						contig.base[pos] = hits[j].base[pos]
						contig.pos[pos] = deepcopy(hits[j].pos[pos])
					contig.qend = hits[j].qend
					contig.qpos.update(hits[j].qpos)
					if strand == '+':
						contig.seqids[seqid][0][1] = hits[j].seqids[seqid][0][1]
					else:
						contig.seqids[seqid][0][0] = hits[j].seqids[seqid][0][0]
				used[j] = 1
				# debug
				#print(vars(contig))
				#input()
		contigs.append(contig)
	contigs.sort(key = lambda x : x.qend)
	contigs.sort(key = lambda x : x.qstart)
	return contigs

# predict genes with introns from the (splitted) sequences
def predict_from_frag(contigs, intron_len, site_composition, mol_type='dna', gencode=1, final_seq='consensus', split_len=None, long_gaps=False):
	# reorganize the hits from different regions
	hits = []
	for contig in contigs:
		for seqid, SEs in contig.seqids.items():
			for SE in SEs:
				hit = deepcopy(contig)
				hit.seqids = {seqid: [SE.copy()]}
				hit.count = {seqid: 1}
				hit.map_pos(seqid, mol_type)
				hits.append(hit)
	hits.sort(key = lambda x : x.qend)
	hits.sort(key = lambda x : x.qstart)
	if split_len:
		hits = merge_split(hits, mol_type)

	# find the clostest regions and merge into one sequences
	contigs = []
	used = {}
	for i in range(len(hits)):
		if used.get(i) is not None:
			continue
		used[i] = 1
		seqid = list(hits[i].seqids.keys())[0]
		strand = hits[i].seqids[seqid][0][2]
		contig = deepcopy(hits[i])
		j = i
		while j < len(hits):
			start = hits[j].seqids[seqid][0][0]
			end = hits[j].seqids[seqid][0][1]
			closest = None
			min_dist = intron_len + 1
			for k in range(j+1, len(hits)):
				if used.get(k) is not None or hits[k].seqids.get(seqid) is None:
					continue
				if hits[k].seqids[seqid][0][2] != strand:
					continue
				if hits[j].qend - hits[k].qstart + 1 > 10:
					continue
				#if hits[j].qend - hits[k].qstart + 1 > 5:
				#	if hits[j].pos.get(hits[k].qstart) is None or hits[k].pos.get(hits[j].qend) is None:
				#		continue
				#	if mol_type != 'dna' and (hits[k].seqids[seqid][0][0] - hits[j].seqids[seqid][0][0]) % 3 != 0:
				#		continue
				#	if hits[j].pos[hits[k].qstart] == hits[k].pos[hits[k].qstart] and hits[j].pos[hits[j].qend] == hits[k].pos[hits[j].qend]:
				#		closest = k
				#		break
				#else:
				if long_gaps and mol_type != 'dna' and (hits[k].seqids[seqid][0][0] - hits[j].seqids[seqid][0][0]) % 3 != 0:
					continue
				if strand == '+':
					dist = hits[k].seqids[seqid][0][0] - end - 1
				else:
					dist = start - hits[k].seqids[seqid][0][1] - 1
				if dist < min_dist:
					if strand == '+' and hits[k].seqids[seqid][0][0] >= start or (strand == '-' and hits[k].seqids[seqid][0][1] <= end):
						closest = k
						min_dist = dist
			# debug
			#print(j, closest, min_dist)
			#input()
			if closest is None:
				break
			if min_dist < intron_len + 1:
				step = 1
				if mol_type != 'dna':
					step = step * 3
				# find the cut points and fix the overlap reference positions
				if hits[j].qend - hits[closest].qstart + 1 > 0:
					max_pos = None
					max_score = 0
					for k in range(hits[j].qend - hits[closest].qstart + 2):
						kscore = 0
						for pos in range(hits[closest].qstart, hits[j].qend - k + 1):
							kscore += site_composition[pos].get(hits[j].base[pos], 0)
						for pos in range(hits[j].qend - k + 1, hits[j].qend + 1):
							kscore += site_composition[pos].get(hits[closest].base[pos], 0)
						if kscore > max_score:
							max_pos = hits[j].qend - k + 1
							max_score = kscore
					qpos = max_pos - 1
					while (hits[j].pos.get(qpos) is None or hits[j].pos.get(qpos, {}) == {1: None, 2: None, 3: None}) and qpos > hits[j].qstart:
						qpos -= 1
					if mol_type == 'dna':
						if strand == '+':
							contig.seqids[seqid][-1][1] = hits[j].pos[qpos]
						else:
							contig.seqids[seqid][-1][0] = hits[j].pos[qpos]
					else:
						if strand == '+':
							contig.seqids[seqid][-1][1] = hits[j].pos[qpos][3]
						else:
							contig.seqids[seqid][-1][0] = hits[j].pos[qpos][3]
					for pos in range(qpos+1, hits[j].qend + 1):
						del contig.seq_comp[pos]
					contig.qend = qpos
					qpos = max_pos
					while (hits[closest].pos.get(qpos) is None or hits[closest].pos.get(qpos, {}) == {1: None, 2: None, 3:None}) and qpos < hits[closest].qend:
						qpos += 1
					if mol_type == 'dna':
						if strand == '+':
							hits[closest].seqids[seqid][0][0] = hits[closest].pos[qpos]
						else:
							hits[closest].seqids[seqid][0][1] = hits[closest].pos[qpos]
					else:
						if strand == '+':
							hits[closest].seqids[seqid][0][0] = hits[closest].pos[qpos][1]
						else:
							hits[closest].seqids[seqid][0][1] = hits[closest].pos[qpos][1]
					for pos in range(hits[closest].qstart, qpos):
						del hits[closest].seq_comp[pos]
					hits[closest].qstart = qpos
				# find the cut points and fix the overlap target positions
				kstart = hits[closest].seqids[seqid][0][0]
				kend = hits[closest].seqids[seqid][0][1]
				max_pos = None
				max_score = 0
				if mol_type != 'dna' and (start - kstart) %3 != 0:
					pass
				elif strand == '+' and kstart <= end:
					for k in range(end+1, kstart-1, -step):
						kscore = 0
						for pos in range(kstart, k, step):
							qpos = hits[j].qpos.get(pos)
							if qpos is not None:
								kscore += site_composition[qpos].get(hits[j].base[qpos], 0)
						for pos in range(k, end+1, step):
							qpos = hits[closest].qpos.get(pos)
							if qpos is not None:
								kscore += site_composition[qpos].get(hits[closest].base[qpos], 0)
						if kscore > max_score:
							max_pos = k
							max_score = kscore
					spos = max_pos - 1
					while hits[j].qpos.get(spos) is None and spos > start + 1:
						spos -= 1
					for pos in range(hits[j].qpos[spos]+1, hits[j].qend+1):
						del contig.seq_comp[pos]
					contig.qend = hits[j].qpos[spos]
					contig.seqids[seqid][-1][1] = spos
					spos = max_pos
					while hits[closest].qpos.get(spos) is None and spos < kend:
						spos += 1
					for pos in range(hits[closest].qstart, hits[closest].qpos[spos]):
						del hits[closest].seq_comp[pos]
					hits[closest].qstart = hits[closest].qpos[spos]
					hits[closest].seqids[seqid][0][0] = spos
				elif strand == '-' and kend >= start:
					for k in range(start, kend+2, step):
						kscore = 0
						for pos in range(k, kend+1, step):
							qpos = hits[j].qpos.get(pos)
							if qpos is not None:
								kscore += site_composition[qpos].get(hits[j].base[qpos], 0)
						for pos in range(start, k, step):
							qpos = hits[closest].qpos.get(pos)
							if qpos is not None:
								kscore += site_composition[qpos].get(hits[closest].base[qpos], 0)
						if kscore > max_score:
							max_pos = k
							max_score = kscore
					spos = max_pos
					while hits[j].qpos.get(spos) is None and spos < end:
						spos += 1
					for pos in range(hits[j].qpos[spos]+1, hits[j].qend+1):
						del contig.seq_comp[pos]
					contig.qend = hits[j].qpos[spos]
					contig.seqids[seqid][-1][0] = spos
					spos = max_pos - 1
					while hits[closest].qpos.get(spos) is None and spos > kstart:
						spos -= 1
					for pos in range(hits[closest].qstart, hits[closest].qpos[spos]):
						del hits[closest].seq_comp[pos]
					hits[closest].qstart = hits[closest].qpos[spos]
					hits[closest].seqids[seqid][0][1] = spos
			if long_gaps:
				contig = merge_hits(contig, hits[closest], mol_type, split_len, fill_gap='-')
			else:
				contig = merge_hits(contig, hits[closest], mol_type, split_len, fill_gap='N')
			used[closest] = 1
			j = closest
		del contig.pos
		del contig.qpos
		contigs.append(contig)
	return simplify_contigs(contigs, mol_type, gencode, final_seq == 'consensus_strict')

# assemble the reads and remove the foreign sequences of each species
def hit2aln(sp, hits, stat_num, mol_type='dna', gencode=1, no_assemble=False, overlap_len=30, overlap_pident=98, site_composition={}, outgroup_score=None, outgroup_weight=0.9, split_len=None, intron_len=None, final_seq='consensus'):
	if no_assemble:
		hits.sort(key = lambda x : x.qend)
		hits.sort(key = lambda x : x.qstart)
		contigs = hits
		del hits
	else:
		contigs = assemble_on_ref(hits, overlap_len, overlap_pident, mol_type, split_len)
		del hits
		""" # graph-based assembly, but the results seem not better
		# assemble the contigs until they can not be assembled
		while len(contigs) > 1:
			new_contigs = assemble_on_ref(deepcopy(contigs), overlap_len, overlap_pident, mol_type, split_len)
			if len(new_contigs) < len(contigs):
				contigs = deepcopy(new_contigs)
				del new_contigs
			else:
				del new_contigs
				break
		if len(contigs) > 1 and contigs[0].count.get('too much') is None:
			used = {}
			for contig in contigs:
				for seqid in contig.count.keys():
					used.setdefault(seqid, 0)
					used[seqid] += 1
			for contig in contigs:
				for seqid in contig.count.keys():
					contig.count[seqid] = 1 / used[seqid]
		"""
		contigs = simplify_contigs(contigs, mol_type, gencode)
		stat_num['assembled contigs'] = len(contigs)
		if outgroup_score:
			if len(contigs) > 0:
				contigs = remove_out_seq(contigs, site_composition, outgroup_score, outgroup_weight)
			stat_num['ingroup clean seqs'] = len(contigs)

	if no_assemble and not final_seq.startswith('consensus') and len(contigs) > 0 and contigs[0].count.get('too much') is None:
		if intron_len:
			contigs = predict_from_frag(contigs, intron_len, site_composition, mol_type, gencode, final_seq, split_len)
		else:
			# connect the fragments from the same sequences without introns (e.g., transcripts), between which the gaps are always shorter than 3000 bp / 1000 aa
			contigs = predict_from_frag(contigs, 3000, site_composition, mol_type, gencode, final_seq, split_len, long_gaps=True)
		stat_num['predicted gene seqs'] = len(contigs)

	# pre-process the hits to reduce the memory when not to assemble (too many hits)
	if no_assemble and final_seq != 'all' and len(contigs) > 0:
		if final_seq.startswith('consensus'):
			contigs[0] = simplify_contigs([contigs[0]], mol_type, gencode, final_seq == 'consensus_strict')[0]
		elif final_seq == 'length':
			contigs.sort(key = lambda x : len(x.seq_comp), reverse = True)
		elif final_seq == 'expression':
			contigs.sort(key = lambda x : sum(x.count.values()), reverse = True)
		contigs = [contigs[0]]
	return sp, contigs, stat_num

# detect and remove cross contamination of assembled sequences
def cross_decontam(alnid, spinfos, mol_type='dna', min_overlap=30, min_pident=98, min_expression=0.2, min_fold=5):
	total_count = {}
	for line in open('total_counts.tsv'):
		arr = line.rstrip().split("\t")
		total_count[arr[0]] = int(arr[1])
	to_remove = {}
	for i in range(len(spinfos)):
		spinfo = spinfos[i]
		to_remove[i] = []
		for j in range(len(spinfo[1])):
			spinfos[i][1][j].addRPKM(total_count[spinfo[0]])
	statout = open(os.path.join('stat_info', alnid + '.cross_info.tsv'), 'w')
	statout.write("species1\tquery_start1\tquery_end1\tspecies2\tquery_start2\tquery_end2\toverlap_length\tpident\tpseudoRPKM1\tpseudoRPKM2\texpression_fold\tseq1\tseq2\n")
	for i in range(len(spinfos) - 1):
		sp1, contigs1, stat_num1 = spinfos[i]
		for j in range(i + 1, len(spinfos)):
			sp2, contigs2, stat_num2 = spinfos[j]
			for m in range(len(contigs1)):
				for n in range(len(contigs2)):
					if contigs1[m].qstart > contigs2[n].qstart:
						hit1 = contigs2[n]
						hit2 = contigs1[m]
					else:
						hit1 = contigs1[m]
						hit2 = contigs2[n]
					nident = 0
					end_pos = min(hit1.qend, hit2.qend)
					overlap_len = end_pos - hit2.qstart + 1
					if mol_type == 'dna':
						if overlap_len < min_overlap:
							continue
						for pos in range(hit2.qstart, end_pos+1):
							if hit1.seq_comp.get(pos) is None or hit2.seq_comp.get(pos) is None:
								continue
							if hit1.base[pos] == hit2.base[pos]:
								nident += 1
					else:
						overlap_len = overlap_len * 3
						if overlap_len < min_overlap:
							continue
						for pos in range(hit2.qstart, end_pos+1):
							if hit1.seq_comp.get(pos) is None or hit2.seq_comp.get(pos) is None:
								continue
							for k in [1, 2, 3]:
								if hit1.seq_comp[pos][k] == hit2.seq_comp[pos][k]:
									nident += 1
					pident = nident / overlap_len * 100
					if pident < min_pident:
						continue
					exp_fold = contigs1[m].pseudoRPKM / contigs2[n].pseudoRPKM
					seq1 = ''
					for pos in range(contigs1[m].qstart, contigs1[m].qend+1):
						if mol_type == 'dna':
							seq1 += contigs1[m].base.get(pos, 'n')
						else:
							for k in [1, 2, 3]:
								seq1 += list(contigs1[m].seq_comp.get(pos, {k: {'n': 1}})[k].keys())[0]
					seq2 = ''
					for pos in range(contigs2[n].qstart, contigs2[n].qend+1):
						if mol_type == 'dna':
							seq2 += contigs2[n].base.get(pos, 'n')
						else:
							for k in [1, 2, 3]:
								seq2 += list(contigs2[n].seq_comp.get(pos, {k: {'n': 1}})[k].keys())[0]
					statout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sp1, contigs1[m].qstart, contigs1[m].qend, sp2, contigs2[n].qstart, contigs2[n].qend, overlap_len, pident, contigs1[m].pseudoRPKM, contigs2[n].pseudoRPKM, exp_fold, seq1, seq2))
					# debug
					#if pident < min_pident:
					#	continue
					if exp_fold > min_fold and contigs1[m].pseudoRPKM >= min_expression and n not in to_remove[j]:
						to_remove[j].append(n)
					elif exp_fold < 1 / min_fold and contigs2[n].pseudoRPKM >= min_expression and m not in to_remove[i]:
						to_remove[i].append(m)
	statout.close()
	# debug
	#print(to_remove)
	#input()
	for i, remove_indexes in to_remove.items():
		for m in sorted(remove_indexes, reverse=True):
			spinfos[i][1].pop(m)
		spinfos[i][2]['cross clean seqs'] = len(spinfos[i][1])
	return spinfos

# output the final result alignments and the stat information
def output(alnid, spinfos, aln_len, ref=None, mol_type='dna', gencode=1, no_assemble=False, sep='.', nuclfill='N', protfill='X', final_seq='consensus', keep_seqid=False):
	fseqs = {}
	fprotseqs = {}
	# add the reference sequences
	if ref:
		if mol_type != 'dna':
			if mol_type == 'codon':
				codonref = lib.read_fastx(os.path.join('map', alnid + '.codon.ref.fas'), 'fasta')
				fseqs.update(codonref)
			fprotseqs.update(ref)
		else:
			fseqs.update(ref)

	# print the stat information and generate the final sequences
	statout = open(os.path.join('stat_info', alnid + '.stat_info.tsv'), 'w')
	statout.write("species\t" + "\t".join(spinfos[0][2].keys()) + "\tfinal seqs\n")
	asout = open(os.path.join('stat_info', alnid + '.assemble_info.tsv'), 'w')
	asout.write("sequence ID\tspecies\tread count\tvalid length\ttarget regions\n")
	for spinfo in spinfos:
		sp, contigs, stat_num = spinfo
		seqs = []
		protseqs = []
		seqids = []
		counts = []
		lens = []
		if len(contigs) == 0:
			pass
		elif final_seq.startswith('consensus'):
			# the sequences consisting of more hits have more weights
			if len(contigs) > 1:
				contigs.sort(key = lambda x : sum(x.count.values()), reverse = True)
				contig = deepcopy(contigs[0])
				if contig.count.get('too much') is None:
					for seqid in contig.count.keys():
						contig.count[seqid] = 1
				for hit in contigs[1:]:
					contig = merge_hits(contig, hit, mol_type)
				contigs[0] = simplify_contigs([contig], mol_type, gencode, final_seq == 'consensus_strict')[0]
			if mol_type == 'dna':
				seq = list(nuclfill * aln_len)
				protseq = None
			else:
				seq = list(nuclfill * (3*aln_len))
				protseq = list(protfill * aln_len)
			seq, protseq = contigs[0].print_seq(seq, protseq)
			seqs.append(''.join(seq).replace('n', nuclfill))
			if mol_type != 'dna':
				protseqs.append(''.join(protseq))
			seqids = [contigs[0].seqids]
			counts = [sum(contigs[0].count.values())]
			lens = [len(contigs[0].seq_comp) - sum(1 for value in contigs[0].base.values() if value == '-')]
		elif final_seq == 'all':
			# output all assembled sequences without consensus or selection
			for contig in contigs:
				if mol_type == 'dna':
					seq = list(nuclfill * aln_len)
					protseq = None
				else:
					seq = list(nuclfill * (3*aln_len))
					protseq = list(protfill * aln_len)
				seq, protseq = contig.print_seq(seq, protseq)
				seqs.append(''.join(seq))
				if mol_type != 'dna':
					protseqs.append(''.join(protseq))
				seqids.append(contig.seqids)
				counts.append(sum(contig.count.values()))
				lens.append(len(contig.seq_comp) - sum(1 for value in contig.base.values() if value == '-'))
		else:
			if final_seq == 'expression':
				contigs.sort(key = lambda x : sum(x.count.values()), reverse = True)
			else:
				# for longest sequence
				contigs.sort(key = lambda x : len(x.seq_comp), reverse = True)
			if mol_type == 'dna':
				seq = list(nuclfill * aln_len)
				protseq = None
			else:
				seq = list(nuclfill * (3*aln_len))
				protseq = list(protfill * aln_len)
			seq, protseq = contigs[0].print_seq(seq, protseq)
			seqs.append(''.join(seq))
			if mol_type != 'dna':
				protseqs.append(''.join(protseq))
			seqids = [contigs[0].seqids]
			counts = [sum(contigs[0].count.values())]
			lens = [len(contigs[0].seq_comp) - sum(1 for value in contigs[0].base.values() if value == '-')]
		statout.write("{}\t{}\t{}\n".format(sp, "\t".join(map(str, stat_num.values())), len(seqs)))
		for i in range(len(seqs)):
			if keep_seqid:
				if no_assemble and not final_seq.startswith('consensus'):
					seqid = ' '.join(seqids[i].keys()).split('_PhAlTag_')[0]
					if fseqs.get(seqid) is not None:
						seqid0 = seqid
						j = 1
						seqid = f"{seqid0}.alt{j}"
						while fseqs.get(seqid) is not None:
							j += 1
							seqid = f"{seqid0}.alt{j}"
				else:
					seqid = sep.join([sp, alnid, str(i+1)]) + ' ' + ' '.join(seqids[i].keys()).replace('_PhAlTag_' + sp, '')
			else:
				seqid = sep.join([sp, alnid, str(i+1)])
			fseqs[seqid] = seqs[i]
			if mol_type != 'dna':
				fprotseqs[seqid] = protseqs[i]
			targets = []
			for rawid, SEs in seqids[i].items():
				for SE in SEs:
					targets.append(f"{rawid},{SE[2]},{SE[0]}-{SE[1]}")
			asout.write(f"{seqid.split()[0]}\t{sp}\t{counts[i]}\t{lens[i]}\t{';'.join(targets).replace('_PhAlTag_' + sp, '')}\n")
	statout.close()
	asout.close()

	# output the final sequences
	seqout = open(os.path.join('nt_out', alnid + '.fa'), 'w')
	for seqid, seqstr in fseqs.items():
		seqout.write(">{}\n{}\n".format(seqid, seqstr))
	seqout.close()
	if mol_type != 'dna':
		seqout = open(os.path.join('aa_out', alnid + '.fa'), 'w')
		for seqid, seqstr in fprotseqs.items():
			seqout.write(">{}\n{}\n".format(seqid, seqstr))
		seqout.close()

def main(args):
	# read the arguments
	alnid = sys.argv[1]
	cpu = int(sys.argv[3])
	# enter the directory
	os.chdir(sys.argv[2])
	paras = {}
	for line in open('parameters.config'):
		para_name, para = line.strip().split(":\t")
		if para == '[]':
			para = []
		elif para.startswith("['") and para.endswith("']"):
			para = para.strip("[]'").split(', ')
		elif para.startswith("'") and para.endswith("'"):
			para = para.strip("'")
		paras[para_name] = para
	if paras['mol_type'].startswith('dna'):
		mol_type = 'dna'
	else:
		mol_type = paras['mol_type']
	gencode = int(paras['gencode'])
	if paras['split_len'] == 'None':
		split_len = None
	else:
		split_len = int(paras['split_len'])
	trim_pos = int(paras['trim_pos'])
	extend_pos = int(paras['extend_pos'])
	pos_freq = float(paras['pos_freq'])
	if paras['no_assemble'] == 'True':
		no_assemble = True
	else:
		no_assemble = False
	overlap_len = int(paras['overlap_len'])
	overlap_pident = float(paras['overlap_pident'])
	if paras['no_out_filter'] == 'True':
		no_out_filter = True
	else:
		no_out_filter = False
	outgroup = paras['outgroup']
	ingroup = paras['ingroup']
	sep = paras['sep']
	outgroup_weight = float(paras['outgroup_weight'])
	if paras['intron_len'] == 'None':
		intron_len = None
	else:
		intron_len = int(paras['intron_len'])
	if paras['no_cross_species'] == 'True':
		no_cross_species = True
	else:
		no_cross_species = False
	cross_overlap_len = int(paras['cross_overlap_len'])
	cross_overlap_pident = float(paras['cross_overlap_pident'])
	min_exp = float(paras['min_exp'])
	min_exp_fold = float(paras['min_exp_fold'])
	unknow_symbol = paras['unknow_symbol']
	unknow_prot = paras['unknow_prot']
	final_seq = paras['final_seq']
	if paras['no_ref'] == 'True':
		no_ref = True
	else:
		no_ref = False
	if paras['keep_seqid'] == 'True':
		keep_seqid = True
		info_max_seqs = float('inf')
	else:
		keep_seqid = False
		info_max_seqs = int(paras['info_max_seqs'])
	species = paras['species']
	merge = {}
	# read the target sequences and avoid consuming too much memory for reading merged redundances
	for line in open(os.path.join('map', alnid + '.hit_info.tsv')):
		arr = line.rstrip().split("\t")
		merge[arr[1]] = {}
	for line in open('merged_redundance.txt'):
		seqids = line.rstrip().split(' ')
		if merge.get(seqids[0]) is None:
			continue
		for sp in species:
			merge[seqids[0]][sp] = []
		for seqid in seqids:
			sp = seqid.split('_PhAlTag_')[-1].split('_fastx')[0]
			merge[seqids[0]][sp].append(seqid)
	refseqs = lib.read_fastx(os.path.join('map', alnid + '.ref.fas'), 'fasta')
	# calculate the site composition of the alignments and the conservative score of the outgroup sequence(s)
	if mol_type == 'dna':
		site_composition, outgroup_score = get_ref_score(refseqs, outgroup, ingroup, sep, unknow_symbol)
	else:
		site_composition, outgroup_score = get_ref_score(refseqs, outgroup, ingroup, sep, unknow_prot)
	if no_out_filter:
		outgroup_score = None

	# read the hit information
	print("\nReading the hits...")
	hits = {}
	stat_nums = {}
	for sp in species:
		hits[sp] = []
		stat_nums[sp] = {'raw hits': 0}
		if no_assemble and not no_out_filter:
			stat_nums[sp]['ingroup clean seqs'] = 0
	seqs = lib.read_fastx(os.path.join('map', alnid + '.targets.fas'), 'fasta_db')
	for line in open(os.path.join('map', alnid + '.hit_info.tsv')):
		arr = line.rstrip().split("\t")
		if len(merge.get(arr[1], {})) == 0:
			sp = arr[1].split('_PhAlTag_')[-1].split('_fastx')[0]
			stat_nums[sp]['raw hits'] +=1
		else:
			for sp in species:
				if merge[arr[1]][sp]:
					stat_nums[sp]['raw hits'] += len(merge[arr[1]][sp])
	for line in open(os.path.join('map', alnid + '.hit_info.tsv')):
		arr = line.rstrip().split("\t")
		newarr, start_trim, end_trim = pretreat_hits(arr, site_composition, outgroup_score, no_assemble, outgroup_weight, trim_pos, pos_freq)
		if not newarr:
			continue
		if len(merge.get(arr[1], {})) == 0:
			sp = arr[1].split('_PhAlTag_')[-1].split('_fastx')[0]
			if no_assemble and not no_out_filter:
				stat_nums[sp]['ingroup clean seqs'] += 1
			if no_assemble and final_seq.startswith('consensus') and len(hits[sp]) > 0:
				hits[sp][0].consensus(newarr[2:], [arr[1]], seqs[arr[1]].upper(), site_composition, start_trim, end_trim, extend_pos, pos_freq, mol_type, gencode, split_len)
				# debug
				#print(sp, vars(hits[sp][0]))
				#input()
			else:
				hits[sp].append(read_hit(newarr[2:], [arr[1]], seqs[arr[1]].upper(), site_composition, start_trim, end_trim, extend_pos, pos_freq, gencode, split_len, stat_nums[sp]['raw hits'] > info_max_seqs))
		else:
			for sp in species:
				if merge[arr[1]][sp]:
					if no_assemble and not no_out_filter:
						stat_nums[sp]['ingroup clean seqs'] += len(merge[arr[1]][sp])
					if no_assemble and final_seq.startswith('consensus') and len(hits[sp]) > 0:
						hits[sp][0].consensus(newarr[2:], merge[arr[1]][sp], seqs[arr[1]].upper(), site_composition, start_trim, end_trim, extend_pos, pos_freq, mol_type, gencode, split_len)
					else:
						hits[sp].append(read_hit(newarr[2:], merge[arr[1]][sp], seqs[arr[1]].upper(), site_composition, start_trim, end_trim, extend_pos, pos_freq, gencode, split_len, stat_nums[sp]['raw hits'] > info_max_seqs))
	del merge
	del seqs

	# assemble the reads and remove the foreign sequences
	print("\nAssembling the reads and removing putative foreign sequences of each species...")
	sys.setrecursionlimit(100000)
	args_list = []
	kwds = {'mol_type': mol_type, 'gencode': gencode, 'no_assemble': no_assemble, 'overlap_len': overlap_len, 'overlap_pident': overlap_pident, 'site_composition': site_composition, 'outgroup_score': outgroup_score, 'outgroup_weight': outgroup_weight, 'split_len': split_len, 'intron_len': intron_len, 'final_seq': final_seq}
	for sp in species:
		args_list.append((sp, hits[sp], stat_nums[sp]))
	spinfos = lib.run_mp(hit2aln, args_list, cpu, kwds=kwds)

	# remove cross contamination
	if not no_cross_species and len(spinfos) > 1:
		print("\nRemoving cross contamination...")
		spinfos = cross_decontam(alnid, spinfos, mol_type, cross_overlap_len, cross_overlap_pident, min_exp, min_exp_fold)

	# print the final result alignments and the stat information
	print("\nPrinting the results...")
	aln_len = len(list(refseqs.values())[0])
	if no_ref:
		refseqs = {}
	output(alnid, spinfos, aln_len, refseqs, mol_type, gencode, no_assemble, sep, unknow_symbol, unknow_prot, final_seq, keep_seqid)
	print(f"\nComplished alignment of {alnid}")

if __name__ == "__main__":
	main(sys.argv[1:])


