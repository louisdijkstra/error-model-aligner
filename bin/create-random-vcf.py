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

usage = """%prog <length-file>

	<length-file> 	File that contains the length of the chromosomes of 
			the reference used. 

Creates a VCF-file with a predefined number of SNPs and indels.  
"""

# global variables for chromosome data
chrom_label, chrom_length, chrom_prob = [], [], [] 
deletions, insertions, snps = [], [], [] # contain the variants already placed 

def randomPosition():
	"""Returns a random position from a random chromosome. Probability of 
	   being selected depends on the length of the chromosome."""
	rand = random.random()

	low = 0.0
	for i, prob in enumerate(chrom_prob):
		if low <= rand <= low + prob:
			return i, random.randint(0, chrom_length[i] - 1) 
		low += prob 


def possibleDeletion(chrom_index, position, length):
	"""Returns True when the deletion can be placed there, otherwise False"""
	if position + length >= chrom_length[chrom_index]:
		return False

	for deletion in deletions[chrom_index]:
		if position + length < deletion[0]:
			continue
		elif position > deletion[1]:
			continue
		return False
	return True

def possibleInsertion(chrom_index, position):
	"""Returns True when the insertion can be placed there, otherwise False"""
	for deletion in deletions[chrom_index]:
		if deletion[0] <= position <= deletion[1]:
			return False
	return position not in insertions[chrom_index]

def possibleSNP(chrom_index, position):
	"""Returns True when the SNP can be placed there, otherwise False"""
	for deletion in deletions[chrom_index]:
		if deletion[0] <= position <= deletion[1]:
			return False
	return (position not in insertions[chrom_index]) and (position not in snps[chrom_index])

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-d", action="store", dest="n_deletions", default=100, type=int,
				  		help="Number of created deletions per length. (Default=100)") 
	parser.add_option("-i", action="store", dest="n_insertions", default=100, type=int,
				  		help="Number of created insertions per length. (Default=100)")
	parser.add_option("-k", action="store", dest="min_length", default=1, type=int,
				  		help="Minimal indel length. (Default=1)")
	parser.add_option("-l", action="store", dest="max_length", default=1000, type=int,
				  		help="Maximal indel length. (Default=1000)")
	parser.add_option("-s", action="store", dest="n_snps", default=100000, type=int,
				  		help="Number of SNPs. (Default=100000)")
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	# read in the length file
	lengthfile = open(args[0], 'r') 
	for line in lengthfile:
		values = line.split()
		chrom_label.append(values[0].strip())
		chrom_length.append(int(values[1].strip()))
		deletions.append([])
		insertions.append([])
		snps.append([])
		if chrom_label[-1] == 'MT':
			break 

	total_length = float(sum(chrom_length))
	for length in chrom_length:
		chrom_prob.append(float(length) / total_length)

	# print header
	print('##fileformat=VCFv4.1')
	print('##source=create-random-vcf.py')
	print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
	
	# add deletions
	for length in range(options.max_length, options.min_length-1, -1):
		for i in range(options.n_deletions):			
			chrom_index, position = randomPosition()
			while not possibleDeletion(chrom_index, position, length):
				chrom_index, position = randomPosition()
			deletions[chrom_index].append([position, position + length - 1])
			print("%s\t%d\t.\t.\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=%d"
					%(chrom_label[chrom_index], position, length)) 

	# add insertions
	for length in range(options.max_length, options.min_length-1, -1):
		for i in range(options.n_insertions):			
			chrom_index, position = randomPosition()
			while not possibleInsertion(chrom_index, position):
				chrom_index, position = randomPosition()
			insertions[chrom_index].append(position)
			print("%s\t%d\t.\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=%d"
					%(chrom_label[chrom_index], position, length)) 

	# add SNPs
	for i in range(options.n_snps):
		chrom_index, position = randomPosition()
		while not possibleSNP(chrom_index, position):
				chrom_index, position = randomPosition()
		snps[chrom_index].append(position)
		print("%s\t%d\t.\tN\tN\t.\tPASS\t."
					%(chrom_label[chrom_index], position)) 
	
if __name__ == '__main__':
	sys.exit(main())

