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

usage = """%prog <ref.fasta>

	<ref.fasta> 	the reference. 

Outputs the length of each chromosome in the given fasta file
"""

def main():
	parser = OptionParser(usage=usage)
	(options, args) = parser.parse_args()
	
	if (len(args)!=1):
		parser.print_help()
		return 1

	ref_file = open(args[0], 'r')
	chromosome_label =  ref_file.next()[1:].split()[0].strip()	
	n_bps = 0 

	for line in ref_file:
		if line[0] == '>':
			print("%s\t%d"%(chromosome_label, n_bps)) 
			chromosome_label = line[1:].split()[0].strip()	
			n_bps = 0 
		else:
			n_bps += len(line.strip()) 
		
	print("%s\t%d"%(chromosome_label, n_bps)) 		
		 

if __name__ == '__main__':
	sys.exit(main())
