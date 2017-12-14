#!/usr/bin/env python
"""\n This program is to extract the gene if alternative isoforms exist. 
The output file will include only the nonredundant transcripts.\n"""


from Bio import SeqIO 
import argparse
import numpy as np


def longest_isoform(infile):
	"""
	Return a dictionary listing each gene with the number of amino acids (AA) or bases 
	(DNA) from the longest isoform
	"""
	protein_dict = {}
	for record in SeqIO.parse(infile, "fasta"):
		if 'mRNA' in record.id:
			name = record.id.split('mRNA')[0]
		elif '.' in record.id:
			name = record.id.split('.')[0]
		if name not in protein_dict:
			protein_dict[name]=len(str(record.seq))
		elif len(str(record.seq)) > protein_dict[name]:
			protein_dict[name]=len(str(record.seq))
			
	return protein_dict
	

def remove_redundant(infile, type, seq_dict):
	"""
	remove redundant transcripts from a same gene based on the length of sequences
	"""
	update_list = []
	if type == "cds":
		outfile = open("representative_cds.fa", "w")
	elif type == "protein":
		outfile = open("representative_proteins.fa", "w")
	
	for record in SeqIO.parse(infile, "fasta"):
		if 'mRNA' in record.id:
			name = record.id.split('mRNA')[0]
		elif '.' in record.id:
			name = record.id.split('.')[0]
			
		size = len(str(record.seq))
		if seq_dict[name] == size:
			outfile.write(">"+record.id+'\n')
			outfile.write(str(record.seq)+'\n')
			seq_dict[name] = -1
			update_list.append(record.id)
			
	outfile.close()
	return update_list

	
	
def main():
	parser = argparse.ArgumentParser(prog="representative_transcripts.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--fasta', required=True)	
	parser.add_argument('--type', choices=("protein", "cds"), required=True)	
	args = parser.parse_args()	
	
	infile = args.fasta
	type = args.type
	seq_dict = longest_isoform(infile)
	
	updata_list = remove_redundant(infile, type, seq_dict)
	#print updata_list
	
	
if __name__ == "__main__":
	main()
	
	