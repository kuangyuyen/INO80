#!/usr/bin/env python
# encoding: utf-8
usage = """
Usage:
This script will calculate the pearson correlation for all against all
It produces three output files
1. the text file where it contains all the correlation coefficiency values
2. the dendrogram.pdf
3. Heatmap.pdf for the correlation coefficiency values

Example: python correlation_coeff_v2.py -i input -o output
***

Created by Kuangyu Yen on 2013-02-25.
Copyright (c) 2013 __PughLab@PSU__. All rights reserved.
"""

import sys, getopt, csv, re
import numpy as np
import matplotlib.pyplot as plt
import pylab
from scipy.stats.stats import pearsonr
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list


#
# define functions
#
def calculate_coeff(infile, output):
	readinput, data, r, outhandle = csv.reader(open(infile, "rU"), delimiter="\t"), [], [], {}
	title = readinput.next()[1:]
	
	### transpose the data 
	for row in readinput:
		data.append(map(float, row[1:]))
	
	transpose_data = map(list, zip(*data))

	### calculate correlation coefficiency
	i = 0 
	while i < len(transpose_data):
		coeff = []
		j = 0
		while j < len(transpose_data):
			coeff.append(pearsonr(transpose_data[i], transpose_data[j])[0])
			j = j + 1
		r.append(coeff)
		outhandle[title[i]] = coeff
		
		i = i + 1
	
	print "Done calculating correlation coefficiency!"
	
	### rearrange correlation_coefficiency data according to dendrogram 
	### it will generate output text file as well as dendrogram.png
	link = linkage(r, method="average", metric="euclidean")
	dendrogram_list = leaves_list(link)
	 
	dendrogram_title = []
	for row in dendrogram_list:
		dendrogram_title.append(title[int(row)])
	
	out, png_array = csv.writer(open(output, "w"), delimiter="\t"), []
	out.writerow(["gene"]+dendrogram_title)

	for line in dendrogram_title:
		tmp, tmp_array = [], []
		tmp.append(line)
		i = 0
		while i < len(dendrogram_title):
			tmp.append(outhandle[line][dendrogram_list[i]])
			tmp_array.append(outhandle[line][dendrogram_list[i]])
			i = i + 1
		png_array.append(tmp_array)
		out.writerow(tmp)			
	
	out_dendrogram_png = output[:-4] +"_dendrogram.pdf"
	dendrogram(link, labels=title)
	pylab.savefig(out_dendrogram_png)
	
	print "Done generating dendrogram figure"
	return dendrogram_title, np.array(png_array)
	
	
def generate_heatmap_figure(label, data_array, output):
	fig = plt.figure(figsize=(4, 8), dpi=300)
	
	ax = fig.add_subplot(111, axisbg='white')

	cax = ax.imshow(data_array, vmin=0.0, vmax=1.0, cmap=plt.cm.jet, interpolation="nearest")
	ax.set_title("Heatmap for correlation coefficient", fontsize=10)
	ax.set_xticklabels(label, fontsize=8)
	ax.set_yticklabels(label, fontsize=8)

	cbar = fig.colorbar(cax, ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
	cbar.ax.set_yticklabels(["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"])


	outfile = output[:-4] + ".pdf"
	plt.savefig(outfile, dpi=300)

	print "Done generating heatmap figure!"
	plt.close()
		

if __name__ == '__main__':
	if len(sys.argv) < 2 or not sys.argv[1].startswith("-"): sys.exit(usage)

	# get arguments
	optlist, alist = getopt.getopt(sys.argv[1:], 'hi:o:')
	for opt in optlist:
		if opt[0] == "-h": sys.exit(usage)
		elif opt[0] == "-i": infile = opt[1]
		elif opt[0] == "-o": output  = opt[1]
		
	label, data_array = calculate_coeff(infile, output)
	label.insert(0, "")
	generate_heatmap_figure(label, data_array, output)

