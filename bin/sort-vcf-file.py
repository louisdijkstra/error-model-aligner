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
import vcf
import sys
import operator

__author__ = "Louis Dijkstra"

usage = """%prog [options] <vcf-file> <output-vcf-file>

Sorts all records in a VCF file by chromosome and then by position. 
 
	<vcf-file> 		original VCF file
	<output-vcf-file>	name where the sorted VCF records will be stored	
"""

def main():

	parser = OptionParser(usage=usage)
	parser.add_option("-r", action="store", dest="reference", default=None,  
						help="Sorts the VCF records by chromosome as sorted in the .chromosome-lengths file (Default = only autosomes and the two sex chromosomes: 1,2,...,22,X,Y)")
	parser.add_option("-x", action="store", dest="chromosome", default=None, 
						help="Processes only this chromosome.")
	(options, args) = parser.parse_args()
	
	if (len(args)!=2):
		parser.print_help()
		return 1

	list_of_chromosome_labels = [] # list of the labels (contains the order information)
	list_of_records = dict() # list of VCF records per chromosome label
	if options.reference != None: 
		for line in open(options.reference, 'r'):
			chromosome_label = line.split()[0].strip() 
			list_of_chromosome_labels.append(chromosome_label)
			list_of_records[chromosome_label] = [] # initialize dictionary with empty list
		if options.chromosome != None:
			if not options.chromosome in list_of_chromosome_labels:
				print("ERROR: chromosome %s is not in the reference (see option -x)"%options.chromosome)
				return 1 
			else:
				list_of_chromosome_labels = [options.chromosome]
				list_of_records[options.chromosome] = [] 
	else: # default 
		if options.chromosome != None: # only process this chromosome
			list_of_chromosome_labels = [options.chromosome]
			list_of_records[options.chromosome] = [] 
		else: 
			list_of_chromosome_labels = map(str, range(1,23)) + ['X', 'Y']	
			for chromosome_label in list_of_chromosome_labels: 
				list_of_records[chromosome_label] = [] 

	vcf_filename 		= os.path.abspath(args[0])
	vcf_output_filename	= os.path.abspath(args[1])

	vcf_reader 		= vcf.Reader(open(vcf_filename))
	vcf_writer 		= vcf.Writer(open(vcf_output_filename, 'w'), vcf_reader)

	list_of_neglected_chromosomes = set()

	# Walk through all records in the VCF file
	for vcf_record in vcf_reader: 
		if not vcf_record.CHROM in list_of_chromosome_labels: # not the chromosome we are currently interested in.		
			list_of_neglected_chromosomes.add(vcf_record.CHROM)
			continue
		list_of_records[vcf_record.CHROM].append(vcf_record)
		
	# sort the chromosomes by position
	for chromosome_label in list_of_chromosome_labels: 
		list_of_records[chromosome_label].sort(key=operator.attrgetter('POS')) 

	# write them to the new vcf file:
	for chromosome_label in list_of_chromosome_labels: 
		for vcf_record in list_of_records[chromosome_label]:
			vcf_writer.write_record(vcf_record)
		
	# print chromosomes ignored to output
	print("VCF file %s is sorted and stored @ %s"%(args[0], args[1]))
	
	if len(list_of_neglected_chromosomes) == 0: 
		print("No chromosomes were ignored.")
	else: 
		print("The following chromosomes were IGNORED:")
		for label in list_of_neglected_chromosomes:
			print("\t%s"%label)

	vcf_writer.close()

if __name__ == '__main__':
	sys.exit(main())

