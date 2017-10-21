#!/usr/bin/env python
"""\n This program is to count how many genes are supported by
mRNA evidence, with a cutoff using FPKM (e.g. 0.5). The input 
file is the output file from featureCounts. """


import argparse


def fpkm_cal(dataframe, lib_size, cutoff):
	all, exp = [], []
	for x in xrange(len(dataframe)):
		geneid = (dataframe[x][0]).split('mRNA')[0]
		all.append(geneid)
		fpkm = float(dataframe[x][2])*1e9/(lib_size*dataframe[x][1])
		if fpkm >= cutoff:
			exp.append(geneid)
	all = list(set(all))
	exp = list(set(exp))
	return len(all), len(exp)


def main():
	parser = argparse.ArgumentParser(prog="parse_featurecounts.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input_file', help="output file from featureCounts", required=True)
	parser.add_argument('-c', '--fpkm_cutoff', help="the minimum FPKM for gene expression", default=0.5)
	args = parser.parse_args()
    
	infile = open(args.input_file, "r")
	cutoff = float(args.fpkm_cutoff)

	last_item = ''
	length = 0
	readsc = 0
	lib_size = 0
	dataframe = []

	for line in infile:
		if 'scf' in line:
			line = line.rstrip()
			mRNA = (line.split('\t')[0]).split(':')[0]
			if mRNA == last_item or last_item == '':
				length += int(line.split('\t')[5])
				readsc += int(line.split('\t')[6])
			else:
				dataframe.append([last_item, length, readsc])
				length = int(line.split('\t')[5])
				readsc = int(line.split('\t')[6])		
			last_item = mRNA	

	dataframe.append([last_item, length, readsc])
	
	for x in xrange(len(dataframe)):
		lib_size += dataframe[x][2]
	
	all, exp = fpkm_cal(dataframe, lib_size, cutoff)
	ratio = float(exp)/all*100
	print "total genes: %d, expressed genes: %d (%.2f%%)" % (exp, all, ratio)


if __name__ == "__main__":
	main()
