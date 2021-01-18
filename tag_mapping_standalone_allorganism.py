#!/usr/bin/env python
# encoding: utf-8
"""
python tag_mapping_standalone_allorganism.py -s suffix -r ref -m model -u up_lim -d dn_lim

This script is modified from Zhenhai's ye_plusone_mapping.py 
It maps tags into any coor file (the coor file has to be in the same format as Tsukiyama et. al. defined plus one dyads.)
It maps tags with consideration of orientation of your coor


Usage:
Put all experimental data files (genetrack shifted tag index file in format desceibed below) into a folder; 
Make sure they have same extentioin (.txt or .tab).


***
suffix:		.txt or .tab
ref:		the coor file for mapping
model: choose chr information for different model organisms: yeast, human, mouse, fly
up_lim:		-300 if you wanna map to 300bp upstream; 100 if you wanna map from downstream 100bp
dn_lim: 	downstream limitation.

example:	python tag_mapping_standalone_allorganism.py -s .tab -r ref.txt -m model -u -300 -d 500
This command will map your tags to upstream 300 bp to downstream 500 bp from cooridates. (strand information of features has been taken care)

***

Created by Zhenhai Zhang on 2009-11-11.
Modified by Kuangyu Yen on 2011-11-03, 2012-04-09
Copyright (c) 2009 The Pennsylvania State Univ.. All rights reserved.
"""

import sys, os, csv, re
import random as rand
import operator
from genome_info import *
from itertools import ifilter
from socket import gethostname
from numpy import mean, array, zeros, ones
from math import sqrt, pi, exp
from time import time
from socket import gethostname
from datetime import datetime

sep="\t"


#
# define functions
#
def processParas(para_list, **keywords):
	# process the parameter information, all parameter start with "-paraname"
	# return a dictionary (paraName, paraValue)
	# remove the first parameter which is the program name
	para_list = para_list[1 :]
	kwgs, values = para_list[ :: 2], para_list[1 :: 2]
	if len(kwgs) != len(values):
		print "number of keywords and values does not equal"
		sys.exit(0)
	
	kwgs = map(lambda x : keywords[x[1 :]], kwgs)
	values = map(evalValues, values)
	return dict(zip(kwgs,values))
	
def evalValues(v):
	# Evaluate strings and return a value corresponding to its real type (int, float, list, tuple)
	try:	return eval(v)
	except:	return v

def getParas(my_dict, *args):
		if len(args) == 1:	return my_dict[args[0]]
		else:	return (my_dict[arg] for arg in args)

def getFiles(my_folder = os.getcwd(), suffix = ""):
	"""
	getting all files under current directory;
	if suffix is not empty, only files with given suffix will be return
	"""
	tmp_result, result = [], []
	def checkSuffix(f):
		return f.endswith(suffix) and os.path.isfile(f)

	for dirpath, dirname, filename in os.walk(my_folder):		#os.walk generate 3-tuple (dirpath, dirnames, filenames)
		tmp_result = filename

	#print "tmp_files:", tmp_result

	if suffix != "":
		result = [x for x in ifilter(checkSuffix, tmp_result)]
		return result
	else:
		return tmp_result
	print "getFile1"
		
def getFiles(suffix):
	all_files = os.listdir(os.getcwd())
	return [x for x in all_files if x.endswith(suffix)]
	print "getFile2"

def getFiles(folder, suffix):
	all_files = os.listdir(os.getcwd() + "/" + folder)
	return [x for x in all_files if x.endswith(suffix)]
	print "getFile3"


def swapExt(s):
	if s == ".txt":
		return ".tab"
	if s == ".tab":
		return ".txt"

def int2chrom(i,  prefix = "chr"):	
	if i < 10:
		return prefix + "0" + str(i)
	else:
		return prefix + str(i)


def populateTagOnChrom(f, chrom, chr_len_organism):
	print "populating reads from %s on %s" %(f, chrom)
	# f is a genetrack input file contain lines in format: chrom	index	forward	reverse
	# chrom is chromosome information in format "chr01" or "chr16"
	result = [0] * (chr_len_organism[chrom] + 1)
	def checkChrom(aList):
		return aList[0] == chrom

	reader = csv.reader(open(f, "rU"), delimiter = sep)
	reader.next()		#skip title line
	total_tags = 0
	for aLine in ifilter(checkChrom, reader):
		index, total = int(aLine[1]), int(aLine[2]) + int(aLine[3])
		try:
			result[index] += total
			total_tags += total
		except:
			print index, chr_len_organism[chrom], len(result)
			print aLine

	print "finished %d reads in total" %total_tags
	return result


def generate_any_coors(chrom):
	reader = csv.reader(open(ref, "rU"), delimiter = sep)
	reader.next()		# skip title line
	
	def checkChrom(aList):
	  return aList[3] == chrom

	for row in ifilter(checkChrom, reader):
		gene, strand, coor = row[4], row[5], int(row[-1])
		yield gene, strand, coor


def get_current_list(count_list, coor, strand, chrom_len):
	global up_lim, dn_lim

	start, end, left, right = 0, 0, [], []
	if strand in "WwFf+":
		start, end = coor + up_lim, coor + dn_lim
	else:
		start, end = coor - dn_lim + 1, coor - up_lim + 1

	if start < 0:
		start, left = 0, [0] * abs(start)

	if end > chrom_len:
		#end, right = chrom_len + 1, [0.0] * (end - chrom_len)
		end, right = chrom_len, [0] * (end - chrom_len)

	result = left + count_list[start : end] + right

	if strand in "CcRr-":
		result.reverse()

	return result



if __name__ == "__main__":
	
	# get arguments
	if len(sys.argv) <  11:
		print __doc__
		sys.exit(0)


	dict_args = processParas(sys.argv, s="suffix", r="ref", m="model", u="up_lim", d="dn_lim")
	suffix, up_lim, dn_lim = getParas(dict_args, "suffix", "up_lim", "dn_lim")
	
	ref = sys.argv[4]
	model = sys.argv[6]
	
	if model == "yeast":
	  chr_len_organism = chr_len_yeast
	elif model == "human":
	  chr_len_organism = chr_len_hg18 
	elif model == "mouse":
	  chr_len_organism = chr_len_mm9
	elif model == "fly":
	  chr_len_organism = chr_len_dm3
	else:
	  sys.exit(0)
	
	print chr_len_organism.items()
	
	title = range(up_lim, dn_lim)
	title.insert(0, "gene")
	
	infiles = getFiles("", suffix)
	for infile in infiles:
		outfile = infile[ : - len(suffix)] + "_u" + str(up_lim)  + "_d" + str(dn_lim) + swapExt(suffix)
		writer = csv.writer(open(outfile, "w"), delimiter = sep)
		writer.writerow(title)
		for chrom, chrom_len in chr_len_organism.items():
			print "mapping reads on chromosome: %s to features" %chrom,
			chrom_count_list = populateTagOnChrom(infile, chrom, chr_len_organism)
			
			for gene, strand, coor in generate_any_coors(chrom):
				curr_list = get_current_list(chrom_count_list, coor, strand, chrom_len)
				writer.writerow([gene] + curr_list)
			
			print "done!"
		print "done processing file %s" %infile

	
	
