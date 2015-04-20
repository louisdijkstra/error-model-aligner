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

usage = """%prog [options] <vcf-file> <bam-file> 

	<vcf-file> 	VCF file containing all the variants 
				that underlie the BAM file
	<bam-file>	BAM file containing the alignments
				NOTE: File needs to be sorted 
				and indexed. 

Outputs for every indel in the given VCF file the following lines 
of output: 

line1: {+/-} <indel-length>

'+' in case of an insertion and '-' in case of a deletion. 

line2: <#IS alignments> <#O alignments> <#O alignments with deletions>
 		<#O alignments with insertions>

#IS alignments - number of insert size (IS) alignments
#O alignments - number of overlapping alignments (potential split-reads)
#O alignments with deletions - overlapping alignments with a deletion split
#O alignments with insertions - overlapping alignments with a insertion split

line3: insert sizes
line4: length deletion splits
line5: distance between centerpoints of the indel and the deletion splits
line6: length insertion splits
line7: distance between centerpoints of the indel and the insertion splits
"""

def returnChromosome (chromosome):
	"""Returns an integer denoting the chromosome."""
	if len(chromosome) > 3 and chromosome[:3] == 'chr':
		chromosome = chromosome[3:]
	if chromosome == 'X' or chromosome == 'x':
		return 23
	if chromosome == 'Y' or chromosome == 'y':
		return 24	
	if chromosome == 'MT' or chromosome == 'M':
		return 25
	return int(chromosome)


def processDeletion (vcf_record, bam_reader, search_range = 5000): 
	"""Processes a deletion (vcf record) given a BAM reader."""
	deletion = Deletion(vcf_record)

	# allocate data 
	n_is		= 0 # number of insert size observations
	n_o 		= 0 # number of overlapping alingments
	n_o_del		= 0 # number of overlapping alignments with deletion splits
	n_o_ins		= 0 # number of overlapping alignments with insertion splits
	insert_sizes 	= [] 
	length_del	= [] # length of the deletion splits
	length_ins 	= [] # length of the insertion splits
	dist_del	= [] # distance between the centerpoints of the deletion split 
	dist_ins 	= [] # distance between the centerpoints of the insertion split

	alignment_dict = defaultdict(list)
	# fetch the alignments in the vicinity of the deletion
	for align in bam_reader.fetch(vcf_record.CHROM, max(0, deletion.start - 1 - search_range), deletion.end + 1 + search_range):
		if align.isize == 0: # alignment is unmapped
			continue

		start = align.pos
		end = align.pos + align.alen - 1 	

		# determine whether the alignment is an overlapping one
		if start < min(deletion.centerpoints) and end > max(deletion.centerpoints):
			# process overlapping alignment here!
			n_o += 1
			insertion_split_present = False 
			deletion_split_present = False 
			i = start
			for (cigar_type, cigar_length) in align.cigar: # walk through the cigar string
				if cigar_type == 1: # insertion
					insertion_split_present = True
					length_ins.append(cigar_length)	
					dist_ins.append(returnMinimumDifference([i], deletion.centerpoints))		
				elif cigar_type == 2: # deletion
					deletion_split_present = True
					length_del.append(cigar_length)
					centerpoints_split = [i + cigar_length / 2]
					if cigar_length % 2 == 0:
						centerpoints_split = [i + cigar_length / 2 - 1, i + cigar_length / 2]
					dist_del.append(returnMinimumDifference(centerpoints_split, deletion.centerpoints))

				i += cigar_length 
			if deletion_split_present: 
				n_o_del += 1 
			if insertion_split_present:
				n_o_ins += 1 	
		else:
			alignment_dict[align.qname].append(align) 
	
	# walk through the paired-end reads
	for qname, alignments in alignment_dict.iteritems():
		if len(alignments) == 2: # paired-end read 	
			# process insert size alignment here!
			align_l, align_r = alignments[0], alignments[1]
			if align_r.pos + align_r.alen - 1 < align_l.pos:
				temp = align_r
				align_r = align_l
				align_l = temp
			interval_segm_start = align_l.pos + align_l.alen
			interval_segm_end  = align_r.pos - 1
			if interval_segm_start <= min(deletion.centerpoints) and interval_segm_end >= max(deletion.centerpoints):
				insert_sizes.append(interval_segm_end - interval_segm_start + 1)
				n_is += 1

	# print data to output
	print("- %d"%deletion.length)
	print("%d %d %d %d"%(n_is, n_o, n_o_del, n_o_ins))
	for i, insert_size in enumerate(insert_sizes):
		print("%d "%insert_size, end = '')
	print()
	for i, l in enumerate(length_del):
		print("%d "%l, end = '')
	print()
	for i, d in enumerate(dist_del):
		print("%d "%d, end = '')
	print()
	for i, l in enumerate(length_ins):
		print("%d "%l, end = '')
	print()
	for i, d in enumerate(dist_ins):
		print("%d "%d, end = '')
	print()


def processInsertion (vcf_record, bam_reader, search_range = 5000): 
	"""Processes an insertion (vcf record) given a BAM reader."""
	insertion = Insertion(vcf_record)

	# allocate data 
	n_is		= 0 # number of insert size observations
	n_o 		= 0 # number of overlapping alingments
	n_o_del		= 0 # number of overlapping alignments with deletion splits
	n_o_ins		= 0 # number of overlapping alignments with insertion splits
	insert_sizes 	= [] 
	length_del	= [] # length of the deletion splits
	length_ins 	= [] # length of the insertion splits
	dist_del	= [] # distance between the centerpoints of the deletion split 
	dist_ins 	= [] # distance between the centerpoints of the insertion split

	alignment_dict = defaultdict(list)
	# fetch the alignments in the vicinity of the deletion
	for align in bam_reader.fetch(vcf_record.CHROM, max(0, insertion.position - search_range), insertion.position + 1 + search_range):
		if align.isize == 0: # alignment is unmapped
			continue
		
		start = align.pos
		end = align.pos + align.alen - 1 

		# determine whether the alignment is an overlapping one
		if start < insertion.position and end  > insertion.position:
			# process overlapping alignment here!
			n_o += 1
			insertion_split_present = False 
			deletion_split_present = False 
			i = align.pos
			for (cigar_type, cigar_length) in align.cigar: # walk through the cigar string
				if cigar_type == 1: # insertion
					insertion_split_present = True
					length_ins.append(cigar_length)	
					dist_ins.append(abs(i - insertion.position - 1))		
				elif cigar_type == 2: # deletion
					deletion_split_present = True
					length_del.append(cigar_length)
					centerpoints_split = [i + cigar_length / 2]
					if cigar_length % 2 == 0:
						centerpoints_split = [i + cigar_length / 2 - 1, i + cigar_length / 2]
					dist_del.append(returnMinimumDifference([insertion.position], centerpoints_split))

				i += cigar_length 
			if deletion_split_present: 
				n_o_del += 1 
			if insertion_split_present:
				n_o_ins += 1 	
		else:
			alignment_dict[align.qname].append(align) 
	
	# walk through the paired-end reads
	for qname, alignments in alignment_dict.iteritems():
		if len(alignments) == 2: # paired-end read 	
			# process insert size alignment here!
			align_l, align_r = alignments[0], alignments[1]
			if align_r.pos + align_r.alen - 1 < align_l.pos:
				temp = align_r
				align_r = align_l
				align_l = temp
			interval_segm_start = align_l.pos + align_l.alen
			interval_segm_end  = align_r.pos - 1
			if interval_segm_start <= insertion.position and interval_segm_end >= insertion.position:
				insert_sizes.append(interval_segm_end - interval_segm_start + 1)
				n_is += 1

	# print data to output
	print("+ %d"%insertion.length)
	print("%d %d %d %d"%(n_is, n_o, n_o_del, n_o_ins))
	for i, insert_size in enumerate(insert_sizes):
		print("%d "%insert_size, end = '')
	print()
	for i, l in enumerate(length_del):
		print("%d "%l, end = '')
	print()
	for i, d in enumerate(dist_del):
		print("%d "%d, end = '')
	print()
	for i, l in enumerate(length_ins):
		print("%d "%l, end = '')
	print()
	for i, d in enumerate(dist_ins):
		print("%d "%d, end = '')
	print()


def liesInInterval(x, minimum, maximum):
	if minimum != None:
		if x < minimum: 
			return False
	if maximum != None:
		if x > maximum:
			return False
	return True
		

def main():
	parser = OptionParser(usage=usage)
	
	parser.add_option("--deletions-only", action="store_true", dest="deletions_only", default=False, 
				  		help="Only deletions are processed.")
	parser.add_option("--insertions-only", action="store_true", dest="insertions_only", default=False, 
				  		help="Only insertions are processed.")
	parser.add_option("-k", action="store", dest="min_length", default=None, type=int,
				  		help="Minimal length of an indel to be considered. (Default = no minimum)")
	parser.add_option("-l", action="store", dest="max_length", default=None, type=int,
				  		help="Maximal length of an indel to be considered. (Default = no maximum)")
	parser.add_option("-r", action="store", dest="search_range", default=5000, type=int, 
						help="Range to search for potentially relevant reads (Default = 5000 bp)")
	(options, args) = parser.parse_args()

	if (len(args)!=2):
		parser.print_help()
		return 1

	vcf_reader = vcf.Reader(open(args[0]))
	bam_reader = pysam.Samfile(args[1], "rb")
	
	if options.deletions_only: 
		for vcf_record in vcf_reader:  
			if isDeletion(vcf_record):
				if liesInInterval(returnIndelLength(vcf_record), options.min_length, options.max_length): 
					processDeletion(vcf_record, bam_reader, search_range=options.search_range)

	if options.insertions_only: 
		for vcf_record in vcf_reader:  
			if isInsertion(vcf_record):
				if liesInInterval(returnIndelLength(vcf_record), options.min_length, options.max_length): 
					processInsertion(vcf_record, bam_reader, search_range=options.search_range)

	if not options.deletions_only and not options.insertions_only:
		for vcf_record in vcf_reader:  
			if isDeletion(vcf_record):
				if liesInInterval(returnIndelLength(vcf_record), options.min_length, options.max_length): 
					processDeletion(vcf_record, bam_reader, search_range=options.search_range)
			elif isInsertion(vcf_record):
				if liesInInterval(returnIndelLength(vcf_record), options.min_length, options.max_length): 
					processInsertion(vcf_record, bam_reader, search_range=options.search_range)

	bam_reader.close() 

if __name__ == '__main__':
	sys.exit(main())

