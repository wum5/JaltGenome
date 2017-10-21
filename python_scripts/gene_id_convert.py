#!/usr/bin/env python
"""This script is to convert the gene ids generated from MAKER pipeline to better
readable gene ids in genome assembly, such as Scaff0000g00101.1. """


import argparse
from Bio import SeqIO


def parse_maker_id(old_id):
	elements = old_id.split('-')
	for key in elements:
		if 'scf' in key:
			new_id = key[:3]+key[-5:]
			return new_id


def sort_geneid(old_order):
	new_order = [i[0] for i in sorted(enumerate(old_order), key=lambda x:x[1])]
	gene_suffix = [x*10+10 for x in new_order]
	return gene_suffix


def order_index(gene_suffix):
	if 0 < gene_suffix < 100:
		return '000'+str(gene_suffix)
	elif 100 <= gene_suffix < 1000:
		return '00'+str(gene_suffix)
	elif 1000 <= gene_suffix < 10000:
		return '0'+str(gene_suffix)
	elif 10000 <= gene_suffix < 100000:
		return str(gene_suffix)


def create_checklist(infile, gfffile, species, version = '.1'):
	outfile = open('ids_checklist.txt', "w")
	last_scaff = ''
	for record in SeqIO.parse(infile, "fasta"):
		with open(gfffile, "r") as handle:
			for line in handle:
				line = line.rstrip()
				if record.id in line:
					start = int(line.split('\t')[3])
					scaffid = parse_maker_id(record.id)
					if scaffid == last_scaff:
						gene_oldid.append(record.id)
						gene_start.append(start)
						last_scaff = scaffid
					else:
						if last_scaff != '':
							gene_suffix = sort_geneid(gene_start)
							for x in xrange(len(gene_suffix)):
								gene_newid = species+last_scaff+'g'+order_index(gene_suffix[x])+version
								outfile.write("%s\t%s\n" % (gene_oldid[x], gene_newid))
						gene_oldid = [record.id]
						gene_start = [start]
						last_scaff = scaffid
	outfile.close()

	

def main():
	parser = argparse.ArgumentParser(prog="gene_id_convert.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--fasta_file', required=True)
	parser.add_argument('-g', '--gff_file', required=True)
	args = parser.parse_args()

	infile = args.fasta_file
	gfffile = args.gff_file
	
	species = 'Jalsin_'
	create_checklist(infile, gfffile, species, version = '.1')

	
	

if __name__ == "__main__":
	main()
	
	