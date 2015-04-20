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
import math

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))[:-3] + 'python')

from Indel import *
from Alignments import *
from BAMProcessor import *

__author__ = "Louis Dijkstra"

usage = """%prog [options] <vcf-file> <bam-file> 

	<vcf-file> 	VCF file containing all the variants 
				that underlie the BAM file
	<bam-file>	BAM file containing the alignments
				NOTE: File needs to be sorted 
				and indexed. 
	
Estimates the 'epsilon_p' parameter of the model; the probability 
that there an alignment contains no split while an indel is present.  

The output is as follows: 

	l_1 	e_1	c_1
	l_2	e_2	c_2
	.	.	.	
	.	.	.	
	.	.	.	
	l_n	e_n	c_n

where the first column denotes the length of the indel, the second the estimates
and the third the confidence interval (e_i +/- c_i). 

The same is repeated for the insertions. 
"""

def processDeletion (vcf_record, bam_reader, search_range = 5000): 
	"""Processes a deletion (vcf record) given a BAM reader."""
	deletion = Deletion(vcf_record)

	# allocate data 
	n		= 0 # number of overlapping alingments
	n_split		= 0 # number of overlapping alignments with a split 

	# fetch the alignments in the vicinity of the deletion
	for align in bam_reader.fetch(vcf_record.CHROM, max(0, deletion.start - 1 - search_range), deletion.end + 1 + search_range):
		if align.isize == 0: # alignment is unmapped
			continue

		start = align.pos
		end = align.pos + align.alen - 1 	

		# determine whether the alignment is an overlapping one
		if start < min(deletion.centerpoints) and end > max(deletion.centerpoints):
			# process overlapping alignment here!
			n += 1
			split_present = False
			i = start
			for (cigar_type, cigar_length) in align.cigar: # walk through the cigar string
				if cigar_type == 1: # insertion
					split_present = True	
				elif cigar_type == 2: # deletion
					split_present = True	
				i += cigar_length 
			if split_present: 
				n_split += 1 
	
	return n, n - n_split

def processInsertion (vcf_record, bam_reader, search_range = 5000): 
	"""Processes a deletion (vcf record) given a BAM reader."""
	insertion = Insertion(vcf_record)

	# allocate data 
	n		= 0 # number of overlapping alingments
	n_split		= 0 # number of overlapping alignments with a split 

	# fetch the alignments in the vicinity of the deletion
	for align in bam_reader.fetch(vcf_record.CHROM, max(0, insertion.position - search_range), insertion.position + 1 + search_range):
		if align.isize == 0: # alignment is unmapped
			continue

		start = align.pos
		end = align.pos + align.alen - 1 	

		# determine whether the alignment is an overlapping one
		if start < insertion.position and end  > insertion.position:
			# process overlapping alignment here!
			n += 1
			split_present = False
			i = start
			for (cigar_type, cigar_length) in align.cigar: # walk through the cigar string
				if cigar_type == 1: # insertion
					split_present = True	
				elif cigar_type == 2: # deletion
					split_present = True	
				i += cigar_length 
			if split_present: 
				n_split += 1 
	
	return n, n - n_split
	
def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
						help="Range to search for potentially relevant reads (Default = 5000 bp)")
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Prints regularly how many variants have been processed.")
	(options, args) = parser.parse_args()

	if (len(args)!=2):
		parser.print_help()
		return 1

	vcf_reader = vcf.Reader(open(args[0]))
	bam_reader = pysam.Samfile(args[1], "rb")

	n_obs_del 	= defaultdict(float) # number of observations per length 
	no_split_del 	= defaultdict(float) # number of observations with no split per length
	n_obs_ins 	= defaultdict(float) # number of observations per length 
	no_split_ins 	= defaultdict(float) # number of observations with no split per length
	
	n = 0 
	for vcf_record in vcf_reader:  
		n += 1 
		if options.verbose and n % 1000 == 0: 
			print("Processed %d variants"%n)  

		if isDeletion(vcf_record):
			length = returnIndelLength(vcf_record)
			n_obs, no_split = processDeletion(vcf_record, bam_reader, search_range=options.search_range)
			n_obs_del[length] += n_obs
			no_split_del[length] += no_split
		elif isInsertion(vcf_record):
			length = returnIndelLength(vcf_record)
			n_obs, no_split = processInsertion(vcf_record, bam_reader, search_range=options.search_range)
			n_obs_ins[length] += n_obs
			no_split_ins[length] += no_split

	bam_reader.close() 
	
	#print results to screen
	print("*** RESULTS ***")
	print("\nDeletions\n")
	print("length\testimate\tCI\tN")
	print("------\t--------\t--\t-")
	for length in range(1,1001):
		if n_obs_del[length] == 0.0: 
			print("%d\t--------\t--------\t0"%length)
		else:
			estimate = no_split_del[length] / n_obs_del[length]
			ci = 1.96 * (estimate * (1.0 - estimate)) / math.sqrt(n_obs_del[length])
			print("%d\t%lf\t%lf\t%d"%(length, estimate, ci, n_obs_del[length])) 
	

	print("\nInsertions\n")
	print("length\testimate\tCI")
	print("------\t--------\t--")
	for length in range(1,1001):
		if n_obs_ins[length] == 0.0: 
			print("%d\t--------\t--------\t0"%length)
		else:
			estimate = no_split_ins[length] / n_obs_ins[length]
			ci = 1.96 * (estimate * (1.0 - estimate)) / math.sqrt(n_obs_ins[length])
			print("%d\t%lf\t%lf\t%d"%(length, estimate, ci, n_obs_ins[length])) 

if __name__ == '__main__':
	sys.exit(main())

