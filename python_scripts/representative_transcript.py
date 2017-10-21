#!/usr/bin/env python
"""\n This program is to extract the gene if alternative isoforms exist. 
The output file will include only the nonredundant transcripts.\n"""


from Bio import SeqIO 
import argparse
import numpy as np


def longest_isoform(infile):
	"""
	Return a dictionary listing each gene with the number of amino acids from the longest isoform
	"""
	protein_dict = {}
	for record in infile:
		name = record.id.split('mRNA')[0]
		if name not in protein_dict:
			protein_dict[name]=len(str(record.seq))
		elif len(str(record.seq)) > protein_dict[name]:
			protein_dict[name]=len(str(record.seq))
	return protein_dict
	

def remove_redundant(infile, protein_dict, outfile):
	"""
	remove redundant transcripts from a same gene based on the length of protein sequences
	"""
	update_list = []
	outfile = open(outfile, "w")
	
	for record in infile:
		name = record.id.split('mRNA')[0]
		size = len(str(record.seq))
		if protein_dict[name] == size:
			outfile.write(">"+record.id+'\n')
			outfile.write(str(record.seq)+'\n')
			protein_dict[name] = -1
			update_list.append(record.id)
			
	outfile.close()
	return update_list

	
	
def main():
	parser = argparse.ArgumentParser(prog="annotation_by_AED.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--fasta', required=True)	
	parser.add_argument('--gff', required=True)
	parser.add_argument('--outfile', required=True)
	args = parser.parse_args()	
	
	handle = open(args.fasta, "r")
	infile = SeqIO.parse(handle, "fasta")
	outfile = args.outfile
	protein_dict = longest_isoform(infile)
	infile.close()
	
	handle = open(args.fasta, "r")
	infile = SeqIO.parse(handle, "fasta")
	updata_list = remove_redundant(infile, protein_dict, outfile)
	infile.close()

	
if __name__ == "__main__":
	main()
	
	