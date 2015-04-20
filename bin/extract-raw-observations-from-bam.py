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
import vcf
import pysam
import numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')

from Indel import *
from Alignments import *
from BAMProcessor import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <bam-file> <read-length> <result-dir>

	<bam-file>	BAM file (sorted & indexed)
	<read-length>	length of a single read
	<result-dir>	directory in which the results are stored
		
The output is stored in the directory <result-dir>. The output is of
the following form: 

	histogram-data.{insert-sizes|length-deletion-splits|length-insertion-splits|meta}

The extension denotes the type of data that is contained in the file:
	
	.insert-sizes - the insert sizes 
	.length-deletion-splits - the length of the deletion splits found
	.length-insertion-splits - the length of the insertion splits found

Every of the aboves files is organized in two column (tab-seperated): 

	x_1	c_1
	x_2	c_2
	...	...
	x_n	c_n

where x_1 is the minimal value found and x_n is the maximum value found. (Note: x_{i+1} = x_i + 1). 
c_i is the count for x_i. 

In case there are no such observations, the file remains empty. 

The file extension .meta contains the following values: 

	<# alignments> <# alignments with deletion splits> <# alignments with insertion splits>
"""

class CountData:

	def __init__(self):
		self.n = 0 # number of data points added
		self.min_value = float('Inf')
		self.max_value = float('-Inf')
		self.count = defaultdict(int) # count data is stored here	

	def add(self, value, count):
		self.n += count 
		if value < self.min_value: 
			self.min_value = value 
		if value > self.max_value:
			self.max_value = value 
		self.count[value] += count

	def print(self, output_stream):
		if self.n == 0:
			return 0
		for i in range(self.min_value, self.max_value + 1):
			print("%d\t%d"%(i, self.count[i]), file=output_stream) 


def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Prints regularly how many alignments have been processed.")
	(options, args) = parser.parse_args()

	if (len(args)!=3):
		parser.print_help()
		return 1

	bam_reader 	= pysam.Samfile(args[0], "rb")
	read_length 	= int(args[1])
	result_dir 	= args[2] 
	if result_dir[-1] != '/':
		result_dir += '/'
	
	insert_sizes 		= CountData() # insert size histogram
	length_deletion_splits 	= CountData() # length of deletion splits histogram
	length_insertion_splits = CountData() # length of insertion splits histogram
	# number of alignments found:
	n_align 		= 0 
	n_align_del_split 	= 0 
	n_align_ins_split 	= 0 

	for align in bam_reader.fetch():
		if align.isize == 0: # alignment is unmapped
			continue
		if align.isize - 2 * read_length > 0:
			insert_sizes.add(align.isize - 2 * read_length, 1) 
		n_align += 1
		if options.verbose and n_align % 100000 == 0: 
			print('Having processed %d alignments'%n_align, file=sys.stderr)
		insertion_split_present = False 
		deletion_split_present = False 
		i = align.pos
		for (cigar_type, cigar_length) in align.cigar: # walk through the cigar string
			if cigar_type == 1: # insertion
				insertion_split_present = True
				length_insertion_splits.add(cigar_length, 1)	
			elif cigar_type == 2: # deletion
				deletion_split_present = True
				length_deletion_splits.add(cigar_length, 1)	
			i += cigar_length 
		if deletion_split_present: 
			n_align_del_split += 1 
		if insertion_split_present:
			n_align_ins_split += 1 	
		

	bam_reader.close() 
	
	# print results to file
	insert_sizes.print(open(result_dir + 'histogram-data.insert-sizes', 'w'))
	length_deletion_splits.print(open(result_dir + 'histogram-data.length-deletion-splits', 'w'))
	length_insertion_splits.print(open(result_dir + 'histogram-data.length-insertion-splits', 'w')) 
	print("%d\t%d\t%d"%(n_align, n_align_del_split, n_align_ins_split), file=(open(result_dir + 'histogram-data.meta', 'w'))) 

if __name__ == '__main__':
	sys.exit(main())

