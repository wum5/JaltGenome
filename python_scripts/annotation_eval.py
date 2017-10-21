#!/usr/bin/env python
"""\nThis program is to count how many annotated genes are supported by protein evidence
against curated database (TAIR10, Sprot, Treml) by using BLAST or Interproscan \n"""


import os, argparse, re


def blast_support(infile, cutoff):
	genes = []
	for line in infile:
		line = line.rstrip()
		gene = line.split('\t')[0]
		eval = float(line.split('\t')[10])
		if eval < cutoff:
			genes.append(gene)
	genes = list(set(genes))	
	return len(genes)


def interpro_support(infile, cutoff):
	genes = []
	for line in infile:
		line = line.rstrip()
		gene = line.split('\t')[0]
		try:	
			eval = float(line.split('\t')[8])
		except:
			eval = 1.0 
		if eval < cutoff:
			genes.append(gene)
	genes = list(set(genes))	
	return len(genes)
		
	
def main():
	parser = argparse.ArgumentParser(prog="annotation_eval.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input_file', help="the input BLAST/interproscan output", required=True)
	parser.add_argument('-f', '--input_format', help="input format (blast/interproscan)", required=True)
	parser.add_argument('-e', '--min_eval', help="the evalue cutoff", default=1e-5)
	parser.add_argument('-n', '--gene_number', help="the total number of genes", required=True)
	args = parser.parse_args()
    
	infile = open(args.input_file, "r")
	format = args.input_format
	cutoff = args.min_eval
	total_number = int(args.gene_number)
	
	if format == "blast":
		support_number = blast_support(infile, cutoff)
	elif format == "interproscan":	
		support_number = interpro_support(infile, cutoff)
	else:
		raise InputError("Input format should be blast or interproscan")
	
	ratio = float(support_number)/total_number*100	
	print "number of genes supported: %d out of %d (%.2f%%)" % (support_number, total_number, ratio)
	


if __name__ == "__main__":
	main()	
	