#!/usr/bin/env python
"""This script is to convert the gene ids generated from MAKER pipeline to more
human readable gene IDs in genome assembly according to the relative position of 
genes on the corresponding scaffolds, such as Scaff0010g00010.1. (the first gene
on the scaffold Scaff0010)"""


import argparse, re
from Bio import SeqIO


def parse_scaffold(old_id):
	elements = old_id.split('-')
	for key in elements:
		if 'scf' in key:
			new_id = key[:3]+key[-5:]
			return new_id


def sort_geneid(old_order):
	indices = sorted(xrange(len(old_order)), key=lambda k: old_order[k])
	print indices
	new_order = [0 for x in xrange(len(indices))]
	n = 1
	for x in xrange(len(new_order)):
		new_order[indices[x]] += n
		n += 1	
	gene_suffix = [x*10 for x in new_order]
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
		

def convert_ids(infile, gfffile, species, version = '.1'):
	outfile = open('ids_checklist.txt', "w")
	updated_gff = open('updated_geneids.gff3', "w")
	updated_seq = open('updated_geneids.fasta', "w")
	last_scaff = ''
	
	for record in SeqIO.parse(infile, "fasta"):
		scaffid = parse_scaffold(record.id)
	
		with open(gfffile, "r") as handle:
			for line in handle:
				line = line.rstrip()
				
				pattern = (record.id).split('-mRNA-')[0]
				
				if (pattern+';') in line or (pattern+'-') in line:
					if line.split('\t')[2] == "gene":
						start = int(line.split('\t')[3])
					
						if scaffid == last_scaff:
							gene_oldid.append(record.id)
							gene_seq.append(record.seq)
							gene_start.append(start)
							last_scaff = scaffid
							gff_lines.append(line)
											
						else:
							if last_scaff != '':
								gene_suffix = sort_geneid(gene_start)
							
								for x in xrange(len(gene_suffix)):
									gene_newid = species+last_scaff+'g'+order_index(gene_suffix[x])+version
									outfile.write("%s\t%s\n" % (gene_oldid[x], gene_newid))
								
									updated_seq.write('>'+gene_newid+'\n')
									updated_seq.write(str(gene_seq[x])+'\n')
								
									replaced = re.sub((gene_oldid[x]).split('-mRNA-')[0], gene_newid, gff_lines[x])
									replaced2 = re.sub('-mRNA-', '.', replaced)
									updated_gff.write(replaced2+'\n')
							
								
							gene_seq = [record.seq]
							gene_oldid = [record.id]
							gene_start = [start]
							last_scaff = scaffid
							gff_lines = [line]
							
					else:
						gff_lines[-1] = gff_lines[-1]+'\n'+line
						
	
	gene_suffix = sort_geneid(gene_start)
							
	for x in xrange(len(gene_suffix)):
		gene_newid = species+last_scaff+'g'+order_index(gene_suffix[x])+version
		outfile.write("%s\t%s\n" % (gene_oldid[x], gene_newid))
								
		updated_seq.write(gene_newid+'\n')
		updated_seq.write(str(gene_seq[x])+'\n')
								
		replaced = re.sub((gene_oldid[x]).split('-mRNA-')[0], gene_newid, gff_lines[x])
		replaced2 = re.sub('-mRNA-', '.', replaced)
		updated_gff.write(replaced2+'\n')		

	outfile.close()
	updated_gff.close()

	

def main():
	parser = argparse.ArgumentParser(prog="convert_geneid.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--fasta_file', required=True)
	parser.add_argument('-g', '--gff_file', required=True)
	args = parser.parse_args()

	infile = args.fasta_file
	gfffile = args.gff_file
	
	species = 'Jalsin_'
	convert_ids(infile, gfffile, species, version = '.1')

	
	

if __name__ == "__main__":
	main()
	
	