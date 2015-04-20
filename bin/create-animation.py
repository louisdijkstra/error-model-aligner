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
from matplotlib import animation
from matplotlib import colors
import numpy as np

__author__ = "Louis Dijkstra"

usage = """%prog [options] <histogram-dir> <variant-type> <extension> <x-label>

	<histogram-dir>	The folder with the histogram data as generated by 
				"generate-histogram-files.py"
	<variant-type> 	Either "deletion" or "insertion"
	<extension>	One of the extensions given to the histogram data files,
			e.g., insert-sizes, n-obs
	<x-label>	Label for the x-axis. Spaces are allowed. 
	
Creates an animation of the histogram data in the files. Starts with length 1 and goes 
up to length 1000 bp. 
"""

NULL_COLOR = "#e6550d"
NON_NULL_COLOR = "#2b8cbe"

class HistogramData:
	
	def __init__(self, histogram_filename, normalize=False):
		inputfile 	= open(histogram_filename, 'r')
		self.x_values 	= []
		self.count	= []
		self.max_count 	= float('-Inf')
		for line in inputfile:
			values = map(int, line.split())
			self.x_values.append(values[0])
			self.count.append(values[1])
			if self.max_count < values[1]:
				self.max_count = values[1]
		if normalize: 
			total = float(sum(self.count))
			if total != 0: 
				for i in range(len(self.count)):
					self.count[i] /= total  
				self.max_count /= total 


	def isEmpty(self):
		return self.count == [] 

	def returnMinimumX(self):
		if self.x_values == []:
			return None
		return self.x_values[0]

	def returnMaximumX(self):
		if self.x_values == []:
			return None
		return self.x_values[-1]

	def extendValues(self, min_x, max_x):
		if min_x != self.x_values[0]:
			diff = self.x_values[0] - min_x 
			self.x_values = range(min_x, self.x_values[0]) + self.x_values
			self.count = [0] * diff + self.count
		if max_x != self.x_values[-1]:
			diff = max_x - self.x_values[-1]
			self.count += [0] * diff
			self.x_values += range(self.x_values[-1] + 1, max_x + 1) 

	def setToZero (self, min_x, max_x):
		self.x_values = range(min_x, max_x + 1)
		self.count = [0] * (max_x - min_x + 1)
		
parser = OptionParser(usage=usage)
parser.add_option("--normalize", action="store_true", dest="normalization", default=False,  
				help="Normalizes the data.")
parser.add_option("-k", action="store", dest="min_x", default=None, type=int,  
				help="Minimum x-value (Default: no minimum)")
parser.add_option("-l", action="store", dest="max_x", default=None, type=int,   
				help="Maximum x-value (Default: no maximum)")
parser.add_option("-n", action="store", dest="null_case", default=None,   
				help="File containing the data of the null-case. When used, normalization is turned on automatically. (Default: none used)")
parser.add_option("-o", action="store", dest="outputfilename", default=None,  
				help="Animation is stored under this name (mp4)")
parser.add_option("-w", action="store", dest="width", default=1.0, type=float, 
				help="Bar width. (Default=1.0)")
(options, args) = parser.parse_args()

if (len(args)<4):
	parser.print_help()
	exit(EXIT_FAILURE)

result_dir 	= args[0]
variant_type 	= args[1]
extension 	= args[2]
x_label		= ""
for i in range(3, len(args) - 1):
	x_label += args[i] + ' '
x_label += args[-1]  

if result_dir[-1] != '/':
	result_dir += '/'

if extension[0] != '.':
	extension = '.' + extension 

print("Reading data...")

null_case = None
histograms = []
overall_maximum_count = float('-Inf')
overall_minimum_x = float('Inf')
overall_maximum_x = float('-Inf')

if options.null_case != None:
	options.normalization = True
	null_case = HistogramData(options.null_case, normalize=True) 
	overall_maximum_count = null_case.max_count
	overall_minimum_x = null_case.x_values[0]
	overall_maximum_x = null_case.x_values[-1] 

for i in range(1,1001):
	histogram_filename = result_dir + variant_type + '.length' + str(i) + extension
	histogram = HistogramData(histogram_filename, normalize=options.normalization)
	histograms.append(histogram)
	
	if not histogram.isEmpty():
		if overall_minimum_x > histogram.returnMinimumX():
			overall_minimum_x = histogram.returnMinimumX()
		if overall_maximum_x < histogram.returnMaximumX():
			overall_maximum_x = histogram.returnMaximumX()
		if overall_maximum_count < histogram.max_count: 
			overall_maximum_count = histogram.max_count
print("DONE reading data...")

print("\nMinimal x-value: %d\nMaximal x-value: %d\nMaximal count: %d\n"%(overall_minimum_x, overall_maximum_x, overall_maximum_count))

print("Updating data...")
if options.null_case != None: 
	null_case.extendValues(overall_minimum_x, overall_maximum_x)

for histogram in histograms:
	if histogram.isEmpty():
		histogram.setToZero(overall_minimum_x, overall_maximum_x)
	else: 
		histogram.extendValues(overall_minimum_x, overall_maximum_x)
print("DONE updating data...")

print("Setting up the figure/axes etc...")
fig = plt.figure()

if options.min_x != None: 
	overall_minimum_x = options.min_x
if options.max_x != None:
	overall_maximum_x = options.max_x 
ax = plt.axes(xlim=(overall_minimum_x,overall_maximum_x), ylim=(0,overall_maximum_count))

plt.title(variant_type + " length: 1")
plt.xlabel(x_label)
if options.normalization: 
	plt.ylabel("probability")
else:
	plt.ylabel("frequency")

if options.null_case != None:
	plt.plot(null_case.x_values, null_case.count, color=NULL_COLOR, linewidth=2.0, label="null") 

rects = plt.bar(histograms[0].x_values, histograms[0].count, align='center', width=options.width, color=NON_NULL_COLOR, edgecolor=NON_NULL_COLOR, label=(variant_type) + ' present')
plt.legend() 
print("DONE setting up the figure/axes etc...")

n_hist = 1

def updatefig(*args):
	global n_hist
	print("Plotting histogram #%d"%(n_hist))
	global histograms, n_hist, plt
	plt.title(variant_type + " length: " + str(n_hist+1))
	count = histograms[n_hist].count 
	n_hist += 1
	for rect, h in zip(rects, count):
		rect.set_height(h)
	return rects,	

anim = animation.FuncAnimation(fig, updatefig, frames=1000, interval=1)

if options.outputfilename != None:
	anim.save(options.outputfilename + '.mp4', fps=10, extra_args=['-vcodec', 'libx264']) 
else: 
	plt.show()
