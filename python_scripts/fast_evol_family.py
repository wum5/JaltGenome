#!/usr/bin/env python
"""This script is to pull out the genes within the significantly expanded or contracted
gene families identified by CAFE. For example, python ../../Scripts/fast_evol_family.py 
-i contracted_families_jalt.txt --c updated_MCLout --g pop_genes.txt 
-s sinuosa -o contracted_genes.txt"""

import argparse

def lineage_specific_genes(cluster, species, genefile, outfile, curr_clst):
	genes = cluster.split('\t')
	#print genes
	for x in xrange(len(genes)):
		if species in genes[x]:
			#print genes[x]
			geneid = genes[x].split('|')[1]
			with open(genefile, "r") as handle:
				for line in handle:
					line = line.rstrip()
					if geneid in line:
						outfile.write(curr_clst+'\t'+line+'\n')


def main():
	parser = argparse.ArgumentParser(prog="cafe_input.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--infile', '-i',required=True, \
	help=("the expanded/contracted gene lists estimated by CAFE"))
	parser.add_argument('--outfile', '-o', required=True)
	parser.add_argument('--clusterfile', '-c', required=True, \
	help=("the cluster file used for input of CAFE"))
	parser.add_argument('--genefile', '-g', required=True, \
	help=("the list of all genes in target species with functional annotations"))
	parser.add_argument('--species', '-s', help=("the target species"), required=True)
	args = parser.parse_args()
	
	infile = open(args.infile,"r")
	outfile = open(args.outfile, "w")
	clusterfile = open(args.clusterfile,"r")
	genefile = args.genefile
	species = args.species

	gene_families = infile.readline().split(',')
	#print(len(gene_families))
	#print gene_families

	num = 0
	for line in clusterfile:
		line = line.rstrip()
		num += 1
		curr = 'cluster_'+str(num)
		#print curr
		for x in xrange(len(gene_families)):
			#print key.split('[')[0]
			if curr == gene_families[x].split('[')[0]:
				lineage_specific_genes(line,species,genefile,outfile,curr) 	

	infile.close()
	outfile.close()
	clusterfile.close()

if __name__ == "__main__":
	main()
	
	