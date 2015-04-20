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

usage = """%prog <#chrom> <#bps>

	<#chrom> 	Number of chromosomes
	<#bps> 		Number of base pairs per chromosome

Creates a randomly generated fasta file with a given number of chromosomes
with a fixed number of base pairs. 
"""

basepairs = ['A', 'C', 'G', 'T']	

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-w", action="store", dest="width", default=50, type=int,
				  		help="Number of base pairs per line. (Default=50)") 
	(options, args) = parser.parse_args()
	
	if (len(args)!=2):
		parser.print_help()
		return 1

	n_chromosomes 	= int(args[0])
	n_bps 		= int(args[1])
	
	for c in range(n_chromosomes):
		print(">%d length:%d"%(c+1, n_bps))
		for bp in range(1, n_bps+1):
			print(random.choice(basepairs), end = '')
			if bp % options.width == 0: 
				print() # newline 

	if n_bps % options.width != 0:
		print()
if __name__ == '__main__':
	sys.exit(main())
