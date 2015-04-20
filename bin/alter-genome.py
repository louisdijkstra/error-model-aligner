#!/usr/bin/env python

"""
Copyright (C) 2015 Louis Dijkstra

This file is part of error-model-aligner

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import random

__author__ = "Louis Dijkstra"

usage = """%prog <variants.vcf> <ref.fasta> <ref.chromosome-lengths> 	

	<variants.vcf> 	VCF file that contains the variants that must be 
			incorporated in the reference. File is generated
			by "create-random-vcf.py". NOTE: File must be sorted!
	<ref.fasta> 	the reference. 
	<ref.chromosome-lengths> File with two columns. First column contains 
			the labels of the chromosomes in the fasta file. The other
			contains their lengths in bp. 

Creates a altered genome from a given reference and a list of variants. 
The altered genome is printed to standard output. 
"""

basepairs = ['A', 'C', 'G', 'T']

def returnVariantLength(info):
	info = info.split(';')
	for i, item in enumerate(info):
		if item.split('=')[0] == 'SVLEN':
			return int(item.split('=')[1])

def generateInsertionSequence(variant_length):
	ins_seq = []
	for i in range(variant_length):
		ins_seq.append(random.choice(basepairs))
	return ins_seq

def createSNP(ref_allele):
	valid_basepairs = basepairs[:]
	if ref_allele != 'N':
		valid_basepairs.remove(ref_allele)
	return random.choice(valid_basepairs)

def alterChromosome(index, vcf_records, chromosome, seq, n_char_per_line=50):

	ref, alt = "", ""
	diff, variant_length = 0, 0 
	for i in range(index, len(vcf_records)):
		line = vcf_records[i]
		if line.startswith("#"):
			continue
		values = line.split()
		if values[0].strip() != chromosome:
			index = i 
			break 
		ref = values[3].strip()
		alt = values[4].strip()
		pos = int(values[1].strip())
		if alt == '<DEL>': 
			variant_length = returnVariantLength(values[7])
			seq = seq[:pos+1+diff] + seq[pos+variant_length+1+diff:] # remove 
			diff -= variant_length
		elif alt == '<INS>':
			variant_length = returnVariantLength(values[7])
			seq = seq[:pos+1+diff] + generateInsertionSequence(variant_length) + seq[pos+1+diff:] # add sequence
			diff += variant_length
		else: # SNP
			seq[pos+diff] = createSNP(seq[pos+diff])
			
	print(">%s dna:chromosome length:%d"%(chromosome, len(seq)))
	for i, c in enumerate(seq):
		print(c, end = '')
		if (i+1) % n_char_per_line == 0:
			print('') # new line
	if len(seq) % n_char_per_line != 0:
		print('')

	return index

def updateIndex(index, vcf_records, chromosome):
	for i in range(index, len(vcf_records)):
		line = vcf_records[i]
		if line.startswith("#"):
			continue
		values = line.split()
		if values[0].strip() != chromosome:
			return i 
	return len(vcf_records)

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-w", action="store", dest="width", default=50, type=int,
				  		help="Number of base pairs per line. (Default=50)") 
	parser.add_option("-x", action="store", dest="chromosome", default=None, 
				  		help="Processes only this chromosome. (Default=process all)") 
	(options, args) = parser.parse_args()
	
	if (len(args)!=3):
		parser.print_help()
		return 1

	vcf_file = open(args[0], 'r')
	ref_file = open(args[1], 'r')
	meta_ref_file = open(args[2], 'r')

	chromosomes = []	
	for line in meta_ref_file:
		chromosomes.append(line.split()[0].strip()) 

	if options.chromosome != None:
		if options.chromosome in chromosomes:
			chromosomes = [options.chromosome]
		else:
			print('ERROR: chromosome %s is absent in the given fasta file. Check option -x.'%options.chromosome) 
			return 1

	vcf_records = vcf_file.readlines()
	index = 0 

	# allocate memory
	header 	= ref_file.next()
	chromosome = header[1:].split()[0].strip()
	seq 	= []

	process_this_chromosome = (chromosome in chromosomes)

	# process every chromosome seperately 
	for line in ref_file:
		if line[0] == '>': # new chromosome starts! 
			# process the read chromosome
			if process_this_chromosome: 
				index = alterChromosome(index, vcf_records, chromosome, seq, n_char_per_line=options.width)	
			else: 
				index = updateIndex(index, vcf_records, chromosome)		

			# new chromosome
			header = line 
			chromosome = header[1:].split()[0].strip()
			process_this_chromosome = (chromosome in chromosomes)
			seq = []
		else:
			if process_this_chromosome: 
				seq += line.strip()
	if process_this_chromosome: 
		alterChromosome(index, vcf_records, chromosome, seq, n_char_per_line=options.width)	 

if __name__ == '__main__':
	sys.exit(main())
