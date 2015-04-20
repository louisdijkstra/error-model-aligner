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

usage = """%prog [options] <bam-file> <read-length> 

Outputs all the reads that have insertion splits above the read length. 
"""

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Prints regularly how many alignments have been processed.")
	(options, args) = parser.parse_args()

	if (len(args)!=2):
		parser.print_help()
		return 1

	bam_reader 	= pysam.Samfile(args[0], "rb")
	read_length 	= int(args[1])

	n_align = 0 

	for align in bam_reader.fetch():
		if align.isize == 0: # alignment is unmapped
			continue
		n_align += 1
		if options.verbose and n_align % 100000 == 0: 
			print('Having processed %d alignments'%n_align, file=sys.stderr)
		for (cigar_type, cigar_length) in align.cigar: # walk through the cigar string
			if cigar_type == 1: # insertion
				if cigar_length > read_length: 
					print(align) 		
	bam_reader.close() 
	
if __name__ == '__main__':
	sys.exit(main())

