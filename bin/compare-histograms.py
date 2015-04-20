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
import matplotlib.pyplot as plt

__author__ = "Louis Dijkstra"

usage = """%prog [options] <histogram-file1> <histogram-file2>

	<histogram-file1> Contains the histogram data. First
				column are the labels. Second
				column contains the counts.
			  	Data is displayed in BLUE. 
	<histogram-file1> Contains the histogram data. 
				Data is displayed in RED.
"""

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-a", action="store", dest="alpha", default=0.5, type=float, 
					help="Opacity. (Default=0.5)")
	parser.add_option("-w", action="store", dest="width", default=0.8, type=float, 
					help="Bar width. (Default=0.8)")
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	x_values1 	= []
	counts1		= []

	for line in open(args[0], 'r'):
		values = map(int, line.split()) 
		x_values1.append(values[0])
		counts1.append(values[1])

	x_values2 	= []
	counts2		= []

	for line in open(args[1], 'r'):
		values = map(int, line.split()) 
		x_values2.append(values[0])
		counts2.append(values[1])

	plt.bar(x_values1, counts1, align='center', width=options.width, alpha=options.alpha, color='blue')
	plt.bar(x_values2, counts2, align='center', width=options.width, alpha=options.alpha, color='red')
	plt.show()
	

if __name__ == '__main__':
	sys.exit(main())
