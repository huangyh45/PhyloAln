#!/usr/bin/env python3

import sys
import os
from copy import deepcopy
from Bio.Seq import reverse_complement, translate
from Bio import AlignIO
try:
    import library as lib
except ImportError:
    import lib.library as lib

# class to process the information of the read or sequence hits
class read_hit:

	# inherit the key parameters from the HSPfrag objects
	def __init__(self, hmmresult):
		self.hit_id = hmmresult.hit_id
		self.hit_start = hmmresult.hit_start + 1
		self.hit_end = hmmresult.hit_end
		self.hit_seq = str(hmmresult.hit.seq).upper()
		self.query_start = hmmresult.query_start + 1
		self.query_end = hmmresult.query_end
		self.query_seq = str(hmmresult.query.seq).upper()

	# map the sequence of the target region of target read without assembly
	def map_read_seqstr(self, seqid, seqstr, rep, moltype='dna', gencode=1):
		self.moltype = moltype
		self.gencode = gencode
		self.read_count = rep
		self.seq_comp = {}
		pos = 0
		j = 0
		jpos = 0
		if moltype.startswith('dna'):
			if seqid.endswith('_rev'):
				seqstr = reverse_complement(seqstr)
			for i in range(self.hit_start-1, self.hit_end):
				while self.hit_seq[pos] in ['.', '-']:
					self.seq_comp[self.query_start + j] = {'-': rep}
					j += 1
					pos += 1
					jpos += 1
				if self.query_seq[jpos] not in ['.', '-']:
					self.seq_comp[self.query_start + j] = {seqstr[i]: rep}
					j += 1
				pos += 1
				jpos += 1
		else:
			# for protein or codon-to-nucleotide alignments (codon alignments)
			if seqid.endswith('rev'):
				seqstr = reverse_complement(seqstr)
			i = int(seqid.split('_')[-1].replace('pos', '').replace('rev', ''))
			codon_pos = i - 1
			for protpos in range(self.hit_start-1, self.hit_end):
				i = protpos * 3 + codon_pos
				while self.hit_seq[pos] in ['.', '-']:
					self.seq_comp[self.query_start + j] = {1: {'-': rep}, 2: {'-': rep}, 3: {'-': rep}}
					j += 1
					pos += 1
					jpos += 1
				if self.query_seq[jpos] not in ['.', '-']:
					self.seq_comp[self.query_start + j] = {1: {seqstr[i]: rep}, 2: {seqstr[i+1]: rep}, 3: {seqstr[i+2]: rep}}
					j += 1
				pos += 1
				jpos += 1

	# print the sequence from the hit information
	def print_seq(self, seq, protseq=None):
		for pos, base in self.seq_comp.items():
			if protseq:
				for i in [1, 2, 3]:
					seq[3*pos+i-4] = base['codon'][i-1]
				protseq[pos-1] = base['base']
			else:
				seq[pos-1] = base
		return seq, protseq

	# calculate pseudo RPKM of the hit using read counts
	def addRPKM(self, total_count):
		self.pseudoRPKM = 1000 * 1000000 * self.read_count / total_count / len(self.seq_comp)

# read the alignments with stockholm format for HMM construction (conservative sites for output alignments)
def read_sto(stofile):
	seqs = {}
	alignment = AlignIO.read(stofile, "stockholm")
	cons_pos = alignment._per_col_annotations['reference_annotation']
	for seq in alignment:
		seqstr = str(seq.seq)
		seqs[seq.id] = ''
		for i in range(len(seqstr)):
			if cons_pos[i] == 'x':
				seqs[seq.id] += seqstr[i]
	return seqs

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

# calculate overlap length and percent identity between two hits for cross decontamination
def compare_hits(hitA, hitB, min_overlap=30, min_pident=98):
	if hitA.query_start > hitB.query_start:
		hit1 = deepcopy(hitB)
		hit2 = deepcopy(hitA)
	else:
		hit1 = deepcopy(hitA)
		hit2 = deepcopy(hitB)
	nident = 0
	end_pos = min(hit1.query_end, hit2.query_end)
	if hit1.moltype.startswith('dna'):
		overlap_len = end_pos - hit2.query_start + 1
		if overlap_len < min_overlap:
			return 0, 0
		for i in range(hit2.query_start, end_pos+1):
			if hit1.seq_comp[i] == hit2.seq_comp[i]:
				nident += 1
	else:
		# for codon alignments
		overlap_len = (end_pos - hit2.query_start + 1) * 3
		if overlap_len < min_overlap:
			return 0, 0
		for i in range(hit2.query_start, end_pos+1):
			for k in [0, 1, 2]:
				if hit1.seq_comp[i]['codon'][k] == hit2.seq_comp[i]['codon'][k]:
					nident += 1
	if overlap_len == 0:
		pident = 0
	else:
		pident = nident / overlap_len * 100
	if pident < min_pident:
		return 0, 0
	else:
		return overlap_len, pident

# class to assemble the reads and process the reads mapped to the alignment
class Assembler:

	# construct the site composition of the alignments and calculate the conservative score of the outgroup sequence
	def __init__(self, group_name, alnfile, species, outgroup=None, sep='.'):
		self.species = species
		self.name = group_name
		self.alnfile = alnfile
		stofile = os.path.join('ref_hmm', group_name + '.sto')
		self.aln = read_sto(stofile)
		# process the gap symbols
		for seqid, seqstr in self.aln.items():
			self.aln[seqid] = seqstr.replace('.', '-')
			self.aln[seqid] = seqstr.replace('~', '-')
		for seqid in self.aln.keys():
			if outgroup is None:
				self.outgroup = seqid
				self.outgroup_seq = self.aln.pop(seqid)
				break
			if seqid == outgroup or seqid.startswith(outgroup + sep):
				self.outgroup = seqid
				self.outgroup_seq = self.aln.pop(seqid)
				break
		self.site_composition = {}
		self.outgroup_score = {}
		j = 0
		for i in range(len(self.outgroup_seq)):
			self.site_composition[i+1] = {}
			for seqstr in self.aln.values():
				base = seqstr[i]
				if self.site_composition[i+1].get(base) is None:
					self.site_composition[i+1][base] = 1
				else:
					self.site_composition[i+1][base] += 1
			self.outgroup_score[i+1] = self.site_composition[i+1].get(self.outgroup_seq[i], 0)

	# map the new hits with or without assembly to the original mapped hits
	def map_binning(self, hit, seqid, seqstr, read_count, min_overlap=30, min_pident=98):
		for n in range(len(self.hits)):
			refhit = self.hits[n]
			newhit = deepcopy(refhit)
			end_pos = min(refhit.query_end, hit.query_end)
			pos = 0
			jpos = 0
			if self.moltype.startswith('dna'):
				overlap_len = end_pos - hit.query_start + 1
				if overlap_len < min_overlap:
					continue
				if seqid.endswith('_rev'):
					seqstr = reverse_complement(seqstr)
				nident = 0
				j = hit.hit_start - 1
				for i in range(hit.query_start, end_pos+1):
					while hit.query_seq[jpos] in ['.', '-']:
						j += 1
						pos += 1
						jpos += 1
					if hit.hit_seq[pos] in ['.', '-']:
						nident += calculate_ident(refhit.seq_comp[i], {'-': read_count})
						if newhit.seq_comp[i].get('-') is None:
							newhit.seq_comp[i]['-'] = read_count
						else:
							newhit.seq_comp[i]['-'] += read_count
					else:
						nident += calculate_ident(refhit.seq_comp[i], {seqstr[j]: read_count})
						if newhit.seq_comp[i].get(seqstr[j]) is None:
							newhit.seq_comp[i][seqstr[j]] = read_count
						else:
							newhit.seq_comp[i][seqstr[j]] += read_count
						j += 1
					pos += 1
					jpos += 1
				if nident / overlap_len * 100 < min_pident:
					continue
				i = end_pos + 1
				while i < hit.query_end + 1:
					if hit.hit_seq[pos] in ['.', '-']:
						newhit.seq_comp[i] = {'-': read_count}
						i += 1
					elif hit.query_seq[jpos] in ['.', '-']:
						j += 1
					else:
						newhit.seq_comp[i] = {seqstr[j]: read_count}
						i += 1
						j += 1
					pos += 1
					jpos += 1
				if end_pos < hit.query_end:
					newhit.query_end = hit.query_end
				newhit.read_count += read_count
				self.hits[n] = newhit
				#debug
				#print("overlap length:", str(overlap_len), "pident:", str(nident / overlap_len * 100), "\nMerging\n", vars(refhit), "\nand\n", vars(hit), "seqid: ", seqid, "seq:", seqstr, "\ninto\n", vars(newhit))                        
				#input()
				hit = None
				break
			else:
				# for codon alignments
				overlap_len = (end_pos - hit.query_start + 1) * 3
				if overlap_len < min_overlap:
					continue
				if seqid.endswith('rev'):
					seqstr = reverse_complement(seqstr)
				i = int(seqid.split('_')[-1].replace('pos', '').replace('rev', ''))
				codon_pos = i - 1
				nident = 0
				j = hit.hit_start * 3 + codon_pos - 3
				for i in range(hit.query_start, end_pos+1):
					while hit.query_seq[jpos] in ['.', '-']:
						j += 3
						pos += 1
						jpos += 1
					if hit.hit_seq[pos] in ['.', '-']:
						for k in [1, 2, 3]:
							nident += calculate_ident(refhit.seq_comp[i][k], {'-': read_count})
							if newhit.seq_comp[i][k].get('-') is None:
								newhit.seq_comp[i][k]['-'] = read_count
							else:
								newhit.seq_comp[i][k]['-'] += read_count
					else:
						for k in [1, 2, 3]:
							nident += calculate_ident(refhit.seq_comp[i][k], {seqstr[j+k-1]: read_count})
							if newhit.seq_comp[i][k].get(seqstr[j+k-1]) is None:
								newhit.seq_comp[i][k][seqstr[j+k-1]] = read_count
							else:
								newhit.seq_comp[i][k][seqstr[j+k-1]] += read_count
						j += 3
					pos += 1
					jpos += 1
				if nident / overlap_len * 100 < min_pident:
					continue
				i = end_pos + 1
				while i < hit.query_end + 1:
					if hit.hit_seq[pos] in ['.', '-']:
						newhit.seq_comp[i] = {1: {'-': read_count}, 2: {'-': read_count}, 3: {'-': read_count}}
						i += 1
					elif hit.query_seq[jpos] in ['.', '-']:
						j += 3
					else:
						newhit.seq_comp[i] = {1: {seqstr[j]: read_count}, 2: {seqstr[j+1]: read_count}, 3: {seqstr[j+2]: read_count}}
						i += 1
						j += 3
					pos += 1
					jpos += 1
				if end_pos < hit.query_end:
					newhit.query_end = hit.query_end
				newhit.read_count += read_count
				self.hits[n] = newhit
				#debug
				#print("overlap length:", str(overlap_len), "pident:", str(nident / overlap_len * 100), "\nMerging\n", vars(refhit), "\nand\n", vars(hit), "seqid: ", seqid, "seq:", seqstr, "\ninto\n", vars(newhit))
				#input()
				hit = None
				break

		# if not assembled to any mapped hits, mapped it simply to the alignments
		if hit is not None:
			hit.map_read_seqstr(seqid, seqstr, read_count, self.moltype, self.gencode)
			self.hits.append(hit)

	# map and assemble the hits
	def map_read_seq(self, hmmresults, fasta, moltype='dna', gencode=1, no_assemble=False, min_overlap=30, min_pident=98):
		self.moltype = moltype
		self.gencode = gencode
		self.stat_num = {'raw hits': len(hmmresults)}
		self.hits = []
		targets = []
		readids = []
		hmmhits = []
		hmmresults.sort(key = lambda x : x.query_end)
		hmmresults.sort(key = lambda x : x.query_start)
		for hmmresult in hmmresults:
			if moltype.startswith('dna'):
				hitid = hmmresult.hit_id
				if hitid.endswith('_rev'):
					hitid = hitid.replace('_rev', '')
				targets.append(hmmresult.hit_id)
				readids.append(hitid)
				hmmhits.append(read_hit(hmmresult))
			else:
				targets.append(hmmresult.hit_id)
				readids.append('_'.join(hmmresult.hit_id.split('_')[:-1]))
				hmmhits.append(read_hit(hmmresult))
		reads = lib.read_fastx(fasta, 'fasta', select_list=readids)
		
		# cluster the identical hits first
		if not no_assemble:
			reps = {}
			for i in range(len(targets)):
				identifier = '_'.join([reads[readids[i]].upper(), str(hmmhits[i].query_start), str(hmmhits[i].query_end)])
				if reps.get(identifier) is None:
					reps[identifier] = []
				reps[identifier].append(readids[i])

		# map and assemble each hit
		for i in range(len(targets)):
			hit = hmmhits[i]
			#debug
			#print(vars(hit))
			#input()
			seqstr = reads[readids[i]].upper()
			if no_assemble:
				hit.map_read_seqstr(targets[i], seqstr, 1, moltype, gencode)
				self.hits.append(hit)
			else:
				identifier = '_'.join([seqstr, str(hit.query_start), str(hit.query_end)])
				if reps[identifier][-1] != 'mapped': 
					self.map_binning(hit, targets[i], seqstr, len(reps[identifier]), min_overlap, min_pident)
					reps[identifier].append('mapped')

		# record the hit number
		if not no_assemble:
			self.stat_num['assembled contigs'] = len(self.hits)

	# convert the assembled contig information into simple sequence information (the most common bases instead of SNPs)
	def simplify_hit_info(self):
		for i in range(len(self.hits)):
			for j in self.hits[i].seq_comp.keys():
				bases = self.hits[i].seq_comp[j]
				if self.hits[i].moltype.startswith('dna'):
					self.hits[i].seq_comp[j] = most_common(bases)
				else:
					base1 = most_common(bases[1])
					base2 = most_common(bases[2])
					base3 = most_common(bases[3])
					if base1 == '-' and base2 == '-' and base3 == '-':
						base = '-'
					else:
						try:
							base = translate(base1 + base2 + base3, table=self.hits[i].gencode)
						except:
							base = 'X'
					self.hits[i].seq_comp[j] = {'base': base, 'codon': [base1, base2, base3]}

	# remove the foreign sequences based on conservative scores
	def remove_out_hit(self, weight=1):
		to_remove = []
		for i in range(len(self.hits)):
			hit = self.hits[i]
			score = 0
			outgroup_score = 0
			for pos, seq_comp in hit.seq_comp.items():
				if hit.moltype.startswith('dna'):
					base = seq_comp
				else:
					base = seq_comp['base']
				score += self.site_composition[pos].get(base, 0)
				outgroup_score += self.outgroup_score[pos]
			if score < outgroup_score * weight:
				to_remove.append(i)
		for i in reversed(to_remove):
			self.hits.pop(i)
		self.stat_num['ingroup clean seqs'] = len(self.hits)

	# output the final sequences
	def output_sequence(self, mode='consensus', nuclfill='N', protfill='X'):
		seqs = []
		protseqs = []
		if len(self.hits) == 0:
			return seqs, protseqs
		if mode.startswith('consensus'):
			seq_info = {}
			# the sequences consisting of more hits have more weights
			self.hits.sort(key = lambda x : x.read_count, reverse = True)
			for hit in self.hits:
				for pos, base in hit.seq_comp.items():
					if hit.moltype.startswith('dna'):
						if seq_info.get(pos) is None:
							seq_info[pos] = {}
						if seq_info[pos].get(base) is None:
							seq_info[pos][base] = 1
						else:
							seq_info[pos][base] += 1
					else:
						if seq_info.get(pos) is None:
                                                        seq_info[pos] = {1: {}, 2: {}, 3:{}}
						for k in [1, 2, 3]:
							if seq_info[pos][k].get(base['codon'][k-1]) is None:
								seq_info[pos][k][base['codon'][k-1]] = 1
							else:
								seq_info[pos][k][base['codon'][k-1]] += 1
			if self.moltype.startswith('dna'):
				seq = list(nuclfill * len(self.outgroup_seq))
			else:
				seq = list(nuclfill * (3*len(self.outgroup_seq)))
				protseq = list(protfill * len(self.outgroup_seq))
			for pos, bases in seq_info.items():
				if self.moltype.startswith('dna'):
					if mode == 'consensus' or len(set(bases)) == 1:
						seq[pos-1] = most_common(bases)
				else:
					for i in [1, 2, 3]:
						if mode == 'consensus' or len(set(bases[i])) == 1:
							seq[3*pos+i-4] = most_common(bases[i])
					codon = ''.join(seq[3*pos-3:3*pos])
					if codon == '---':
						protseq[pos-1] = '-'
					else:
						try:
							protseq[pos-1] = translate(codon, table=self.gencode)
						except:
							protseq[pos-1] = 'X'
			seqs.append(''.join(seq))
			if not self.moltype.startswith('dna'):
				protseqs.append(''.join(protseq))
		elif mode == 'all':
			# output all assembled sequences without consensus or selection
			for hit in self.hits:
				if self.moltype.startswith('dna'):
					seq = list(nuclfill * len(self.outgroup_seq))
					protseq = None
				else:
					seq = list(nuclfill * (3*len(self.outgroup_seq)))
					protseq = list(protfill * len(self.outgroup_seq))
				seq, protseq = hit.print_seq(seq, protseq)
				seqs.append(''.join(seq))
				if not self.moltype.startswith('dna'):
					protseqs.append(''.join(protseq))
		else:
			if mode == 'expression':
				self.hits.sort(key = lambda x : x.read_count, reverse = True)
			else:
				# for longest sequence
				self.hits.sort(key = lambda x : len(x.seq_comp), reverse = True)
			if self.moltype.startswith('dna'):
				seq = list(nuclfill * len(self.outgroup_seq))
				protseq = None
			else:
				seq = list(nuclfill * (3*len(self.outgroup_seq)))
				protseq = list(protfill * len(self.outgroup_seq))
			seq, protseq = self.hits[0].print_seq(seq, protseq)
			seqs.append(''.join(seq))
			if not self.moltype.startswith('dna'):
				protseqs.append(''.join(protseq))
		return seqs, protseqs

# a single task for assembly and foreign sequence removal
def generate_assembly(group_name, alnfile, species, hmmresults, np=8, moltype='dna', gencode=1, no_assemble=False, overlap_len=30, overlap_pident=98, no_out_filter=False, outgroup=None, sep='.', outgroup_weight=1, final_seq='consensus'):
	assembler = Assembler(group_name, alnfile, species, outgroup=outgroup, sep=sep)
	assembler.map_read_seq(hmmresults, os.path.join('map_' + species, group_name + '.targets.fa'), moltype=moltype, gencode=gencode, no_assemble=no_assemble, min_overlap=overlap_len, min_pident=overlap_pident)
	assembler.simplify_hit_info()
	if not no_out_filter:
		assembler.remove_out_hit(outgroup_weight)
	# pre-process the hits to reduce the memory when not to assemble (too many hits)
	if no_assemble and len(assembler.hits) > 0:
		if final_seq.startswith('consensus'):
			seq_info = {}
			for hit in assembler.hits:
				for pos, base in hit.seq_comp.items():
					if hit.moltype.startswith('dna'):
						if seq_info.get(pos) is None:
							seq_info[pos] = {}
						if seq_info[pos].get(base) is None:
							seq_info[pos][base] = 1
						else:
							seq_info[pos][base] += 1
					else:
						if seq_info.get(pos) is None:
                                                        seq_info[pos] = {1: {}, 2: {}, 3:{}}
						for k in [1, 2, 3]:
							if seq_info[pos][k].get(base['codon'][k-1]) is None:
								seq_info[pos][k][base['codon'][k-1]] = 1
							else:
								seq_info[pos][k][base['codon'][k-1]] += 1
			start_pos = 0
			end_pos = 0
			for i in range(1, len(assembler.outgroup_seq) + 1):
				if seq_info.get(i) is None:
					continue
				if assembler.moltype.startswith('dna'):
					if final_seq == 'consensus' or len(seq_info[i]) == 1:
						assembler.hits[0].seq_comp[i] = most_common(seq_info[i])
					else:
						assembler.hits[0].seq_comp[i] = None
				else:
					if final_seq == 'consensus' or len(seq_info[i][1]) == 1:
						base1 = most_common(seq_info[i][1])
					else:
						base1 = 'N'
					if final_seq == 'consensus' or len(seq_info[i][2]) == 1:
						base2 = most_common(seq_info[i][2])
					else:
						base2 = 'N'
					if final_seq == 'consensus' or len(seq_info[i][3]) == 1:
						base3 = most_common(seq_info[i][3])
					else:
						base3 = 'N'
					if base1 == '-' and base2 == '-' and base3 == '-':
						base = '-'
					else:
						try:
							base = translate(base1 + base2 + base3, table=self.hits[i].gencode)
						except:
							base = 'X'
					assembler.hits[0].seq_comp[i] = {'base': base, 'codon': [base1, base2, base3]}
				if not start_pos:
					start_pos = i
				end_pos = i
			assembler.hits[0].query_start = start_pos
			assembler.hits[0].query_end = end_pos
		elif final_seq == 'length':
			assembler.hits.sort(key = lambda x : len(x.seq_comp), reverse = True)
		if final_seq != 'all':
			assembler.hits = [assembler.hits[0]]
	return assembler

# assemble the hits and remove putative foreign sequences in multiple processes
def generate_assembly_mp(alns, species, all_hmmres, np=8, moltype='dna', gencode=1, no_assemble=False, overlap_len=30, overlap_pident=98, no_out_filter=False, outgroup=None, sep='.', outgroup_weight=1, final_seq='consensus'):
	print("\nAssembling the reads of {} and removing putative foreign sequences...".format(species))
	assemblers = {}
	args_list = []
	kwds = {'moltype': moltype, 'gencode': gencode, 'no_assemble': no_assemble, 'overlap_len': overlap_len, 'overlap_pident': overlap_pident, 'no_out_filter': no_out_filter, 'outgroup': outgroup, 'sep': sep, 'outgroup_weight': outgroup_weight, 'final_seq': final_seq}
	for group_name, alnfile in alns.items():
		args_list.append((group_name, alnfile, species, all_hmmres[group_name]))
	asbs = lib.run_mp(generate_assembly, args_list, np, kwds=kwds)
	for assembler in asbs:
		assemblers[assembler.name] = assembler
	return assemblers

# detect and remove cross contamination of assembled sequences
def cross_decontam(assemblers, total_reads, min_overlap=30, min_pident=98, min_expression=0.2, min_fold=2):
	to_remove = {}
	for sp1 in assemblers.keys():
		to_remove[sp1] = []
		for i in range(len(assemblers[sp1].hits)):
			assemblers[sp1].hits[i].addRPKM(total_reads[sp1])
	statout = open(os.path.join('stat_info', list(assemblers.values())[0].name + '.cross_info.tsv'), 'w')
	statout.write("species1\tquery_start1\tquery_end1\tspecies2\tquery_start2\tquery_end2\toverlap_length\tpident\tpseudoRPKM1\tpseudoRPKM2\texpression_fold\tseq1\tseq2\n")
	for sp1, assembler1 in assemblers.items():
		for sp2, assembler2 in assemblers.items():
			if sp1 == sp2:
				continue
			for i in range(len(assembler1.hits)):
				for j in range(len(assembler2.hits)):
					overlap_len, pident = compare_hits(assembler1.hits[i], assembler2.hits[j], min_overlap, min_pident)
					if overlap_len:
						exp_fold = assembler1.hits[i].pseudoRPKM / assembler2.hits[j].pseudoRPKM
						seq1 = ''
						for pos in range(assembler1.hits[i].query_start, assembler1.hits[i].query_end+1):
							if assembler1.moltype.startswith('dna'):
								seq1 += assembler1.hits[i].seq_comp[pos]
							else:
								for k in [1, 2, 3]:
									seq1 += assembler1.hits[i].seq_comp[pos]['codon'][k-1]
						seq2 = ''
						for pos in range(assembler2.hits[j].query_start, assembler2.hits[j].query_end+1):
							if assembler2.moltype.startswith('dna'):
								seq2 += assembler2.hits[j].seq_comp[pos]
							else:
								for k in [1, 2, 3]:
									seq2 += assembler2.hits[j].seq_comp[pos]['codon'][k-1]
						statout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(assembler1.species, assembler1.hits[i].query_start, assembler1.hits[i].query_end, assembler2.species, assembler2.hits[j].query_start, assembler2.hits[j].query_end, overlap_len, pident, assembler1.hits[i].pseudoRPKM, assembler2.hits[j].pseudoRPKM, exp_fold, seq1, seq2))
						if exp_fold > min_fold and assembler1.hits[i].pseudoRPKM >= min_expression and j not in to_remove[sp2]:
							to_remove[sp2].append(j)
						elif exp_fold < 1 / min_fold and assembler2.hits[j].pseudoRPKM >= min_expression and i not in to_remove[sp1]:
							to_remove[sp1].append(i)
	statout.close()
	for sp, remove_indexes in to_remove.items():
		for i in sorted(remove_indexes, reverse=True):
			assemblers[sp].hits.pop(i)
		assemblers[sp].stat_num['cross clean seqs'] = len(assemblers[sp].hits)
	return assemblers

# generate conservative reference codon alignments (corresponding to protein alignments for HMM construction)
def generate_codon_ref(assembler):
	codon_seqs = lib.read_fastx(assembler.alnfile, 'fasta')
	prot_seqs = assembler.aln
	prot_seqs[assembler.outgroup] = assembler.outgroup_seq
	seqs = {assembler.outgroup: ''}
	for seqid, protstr in prot_seqs.items():
		seqs[seqid] = ''
		i = 0
		j = 0
		while j < len(protstr) and i + 3 <= len(codon_seqs[seqid]):
			protbase = protstr[j]
			codon = codon_seqs[seqid][i:i+3]
			#debug
			#print(seqid, 'codon', i, codon, 'protein', j, protbase)
			#input()
			if codon == '---':
				if protbase == '-':
					seqs[seqid] += codon
					j += 1
				i += 3
			elif translate(codon, table=assembler.gencode) == protbase:
				seqs[seqid] += codon
				j += 1
				i += 3
			elif protbase == '-':
				seqs[seqid] += '---'
				j += 1
			else:
				i += 3
		if j < len(protstr):
			print("\nError: fail to match the codon columns to the protein alignment HMM: {}-{}-{}!".format(assembler.name, seqid, j))
			return None
	#old version
	#i = 0
	#for j in range(len(assembler.outgroup_seq)):
	#	ifmatch = False
	#	while not ifmatch:
	#		if i + 3 > len(codon_seqs[assembler.outgroup]):
	#			break
	#		ifmatch = True
	#		sites = {}
	#		for seqid, protstr in prot_seqs.items():
	#			protbase = protstr[j]
	#			codon = codon_seqs[seqid][i:i+3]
	#			#debug
	#			print(seqid, 'codon', i, codon, 'protein', j, protbase)
	#			input()
	#			if codon == '---':
	#				if protbase != '-':
	#					ifmatch = False
	#					break
	#				else:
	#					sites[seqid] = codon
	#			elif translate(codon, table=assembler.gencode) != protbase:
	#				ifmatch = False
	#				break
	#			else:
	#				sites[seqid] = codon
	#		i += 3
	#	if not ifmatch:
	#		print("\nError: fail to match the codon columns to the protein alignment HMM: {}-{}!".format(assembler.name, j+1))
	#		return None
	#	for seqid in seqs.keys():
	#		seqs[seqid] += sites[seqid]
	return seqs

# a single task for cross contamination and final sequence output
def cross_and_output(group_name, assemblers, sps, total_reads, moltype='dna', gencode=1, no_assemble=False, no_cross_species=False, min_overlap=30, min_pident=98, min_exp=0.2, min_fold=2, unknow='unknow', final_seq='consensus', no_ref=False, sep='.'):
	if not no_cross_species and len(sps) > 1 and not no_assemble:
		assemblers = cross_decontam(assemblers, total_reads, min_overlap=min_overlap, min_pident=min_pident, min_expression=min_exp, min_fold=min_fold)
	if unknow == 'unknow':
		nuclfill = 'N'
		protfill = 'X'
	else:
		nuclfill = unknow
		protfill = unknow
	seqs = {}
	protseqs = {}
	if not no_ref:
		# prepare the reference alignments
		if moltype.startswith('dna'):
			seqs[assemblers[sps[0]].outgroup] = assemblers[sps[0]].outgroup_seq
			seqs.update(assemblers[sps[0]].aln)
		else:
			protseqs[assemblers[sps[0]].outgroup] = assemblers[sps[0]].outgroup_seq
			protseqs.update(assemblers[sps[0]].aln)
			if moltype == 'codon':
				seqs = generate_codon_ref(assemblers[sps[0]])
				if seqs is None:
					return 1, group_name

	# print the stat information
	statout = open(os.path.join('stat_info', group_name + '.stat_info.tsv'), 'w')
	statout.write("species\t" + "\t".join(assemblers[sps[0]].stat_num.keys()) + "\tfinal seqs\n")
	for sp in sps:
		seq_list, protseq_list = assemblers[sp].output_sequence(mode=final_seq, nuclfill=nuclfill, protfill=protfill)
		statout.write("{}\t{}\t{}\n".format(sp, "\t".join(map(str, assemblers[sp].stat_num.values())), len(seq_list)))
		for i in range(len(seq_list)):
			seqid = sep.join([sp, group_name, str(i+1)])
			seqs[seqid] = seq_list[i]
			if not moltype.startswith('dna'):
				protseqs[seqid] = protseq_list[i]
	statout.close()

	# cross decontamination for consensus without assembly
	if not no_cross_species and len(sps) > 1 and no_assemble and final_seq.startswith('consensus'):
		crossout = open(os.path.join('stat_info', group_name + '.cross_info.tsv'), 'w')
		crossout.write("seqid1\tseqid2\toverlap_length\tpident\tpseudoRPKM1\tpseudoRPKM2\texpression_fold\tseq1\tseq2\n")
		for sp1 in sps:
			for sp2 in sps:
				if sp1 == sp2:
					continue
				seqid1 = None
				seqid2 = None
				for seqid in seqs.keys():
					if seqid.startswith(sp1 + sep):
						seqid1 = seqid
					elif seqid.startswith(sp2 + sep):
						seqid2 = seqid
				if seqid1 is None or seqid2 is None:
					continue
				overlap_len = 0
				nident = 0
				length1 = 0
				length2 = 0
				for i in range(len(seqs[seqid1])):
					if seqs[seqid1][i] != nuclfill:
						length1 += 1
					if seqs[seqid2][i] != nuclfill:
						length2 += 1
					if seqs[seqid1][i] == nuclfill or seqs[seqid2][i] == nuclfill:
						continue
					overlap_len += 1
					if seqs[seqid1][i] == seqs[seqid2][i]:
						nident += 1
				if overlap_len < min_overlap:
					continue
				pident = 0
				if overlap_len > 0:
					pident = nident / overlap_len * 100
				if pident < min_pident:
					continue
				pseudoRPKM1 = 1000 * 1000000 * assemblers[sp1].stat_num['raw hits'] / total_reads[sp1] / length1
				pseudoRPKM2 = 1000 * 1000000 * assemblers[sp2].stat_num['raw hits'] / total_reads[sp2] / length2
				exp_fold = pseudoRPKM1 / pseudoRPKM2
				crossout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seqid1, seqid2, overlap_len, pident, pseudoRPKM1, pseudoRPKM2, exp_fold, seqs[seqid1], seqs[seqid2]))
				if exp_fold > min_fold and pseudoRPKM1 > min_exp:
					seqs.pop(seqid2)
					if not assemblers[sp2].moltype.startswith('dna'):
						protseqs.pop(seqid2)
				elif exp_fold < 1 / min_fold and pseudoRPKM2 > min_exp:
					seqs.pop(seqid1)
					if not assemblers[sp1].moltype.startswith('dna'):
						protseqs.pop(seqid1)
		crossout.close()

	# output the final sequences
	seqout = open(os.path.join('nt_out', group_name + '.fa'), 'w')
	for seqid, seqstr in seqs.items():
		seqout.write(">{}\n{}\n".format(seqid, seqstr))
	seqout.close()
	if not moltype.startswith('dna'):
		seqout = open(os.path.join('aa_out', group_name + '.fa'), 'w')
		for seqid, seqstr in protseqs.items():
			seqout.write(">{}\n{}\n".format(seqid, seqstr))
		seqout.close()
	elif moltype == 'dna_codon':
		lib.trans_seq(os.path.join('nt_out', group_name + '.fa'), os.path.join('aa_out', group_name + '.fa'), gencode, dna_codon_unknow=protfill)
	return 0, group_name

# cross contamination and final sequence output in multiple processes
def cross_and_output_mp(groups, sps, assemblers, total_reads, np = 8, moltype='dna', gencode=1, no_assemble=False, no_cross_species=False, min_overlap=30, min_pident=98, min_exp=0.2, min_fold=2, unknow='unknow', final_seq='consensus', no_ref=False, sep='.'):
	print("\nRemoving cross contamination and printing the final sequences...")
	if not os.path.isdir('stat_info'):
		os.mkdir('stat_info')
	if not os.path.isdir('nt_out'):
		os.mkdir('nt_out')
	if moltype != 'dna' and not os.path.isdir('aa_out'):
		os.mkdir('aa_out')
	if not no_cross_species and len(sps) > 1 and no_assemble:
		if final_seq.startswith('consensus'):
			print("\nWarning: not to assemble and {}, cross decontamination will be conducted after consensus!".format(final_seq))
		else:
			print("\nWarning: not to assemble and {}, cross decontamination will be disabled due to no read counts!".format(final_seq))
	args_list = []
	kwds = {'sps': sps, 'total_reads': total_reads, 'moltype': moltype, 'gencode': gencode, 'no_assemble': no_assemble, 'no_cross_species': no_cross_species, 'min_overlap': min_overlap, 'min_pident': min_pident, 'min_exp': min_exp, 'min_fold': min_fold, 'unknow': unknow, 'final_seq': final_seq, 'no_ref': no_ref, 'sep': sep}
	for group_name in groups:
		args_list.append((group_name, assemblers[group_name]))
	iferrors = lib.run_mp(cross_and_output, args_list, np, kwds=kwds)
	# if errors occur, print the message
	errors = []
	for iferror in iferrors:
		if iferror[0] == 1:
			errors.append(iferror[1])
	if errors:
		print("\nError: fail to output due to invalid reference codon alignments: {}".format(', '.join(errors)))
		sys.exit(1)

# concatenate the short alignments if split the reference alignements
def concatenate(dirname, alns, fill='-'):
	seqs = {}
	total_len = 0
	for aln in alns.keys():
		aln_seqs = lib.read_fastx(os.path.join(dirname, aln + '.fa'), 'fasta')
		for seqid, seqstr in aln_seqs.items():
			if seqid.endswith('.' + aln + '.1'):
				seqid = seqid.replace('.' + aln + '.1', '')
			if seqs.get(seqid) is None:
				seqs[seqid] = fill * total_len + seqstr
			else:
				seqs[seqid] += seqstr
			aln_len = len(seqstr)
		for seqid in seqs.keys():
			if aln_seqs.get(seqid) is None and aln_seqs.get(seqid + '.' + aln + '.1') is None:
				seqs[seqid] += fill * aln_len
		total_len += aln_len
	outfile = open(os.path.join(dirname, 'aln.concatenated.fa'), 'w')
	for seqid, seqstr in seqs.items():
		outfile.write(">{}\n{}\n".format(seqid, seqstr))
	outfile.close()
