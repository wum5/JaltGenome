#!/usr/bin/env python
"""This script is to firstly filter gene trees (the concatenated set of the RAxML bipartition trees) 
by comparing the average bootstrap values against a cutoff threshold. Then, the script is to 
investigate the unsolved internode in my investiated speices tree,
python parse_internode.py -i genetrees.tre -b 80 -s C.annuumm J.sinuosa S.lycopersicum S.tuberosum
here species id are ordered as A,B,C1,C2. 
"""

import argparse, re
import numpy as np
from Bio import Phylo
from cStringIO import StringIO


def bootstraps_extract(tree):
	rawlist = re.findall('\[\d+\]',tree)
	outlist = [int(x[1:][:-1]) for x in rawlist]
	return outlist


def filterTree(bipartition_trees, cutoff):
	outfile = open("trees_bootstrap_"+str(cutoff), "w")
	with open(bipartition_trees, "r") as trees:
		for line in trees:
			tree = line.rstrip()
			bootstraps = bootstraps_extract(tree)
			if np.mean(bootstraps) >= cutoff:
			#if min(bootstraps) >= cutoff:
				outfile.write(tree+'\n')
	outfile.close()
				
				
def parse_tre(infile, species):
	"""
	in this case, we are looking for the frequency of internodes support 1. (A,(B,(C1,C2)));
	((A,(C1,C2)),B); or ((A,B),(C1,C2)).
	"""
	pattern0 = set([species[2],species[3]])
	pattern1 = set([species[1],species[2],species[3]])
	pattern2 = set([species[0],species[2],species[3]])
	pattern3 = set([species[0],species[1]])
	
	total, case1, case2, case3 = 0,0,0,0
	
	for line in infile:
		treedata = line.rstrip()
		handle = StringIO(treedata)
		tree = Phylo.read(handle, "newick")
		signature = []
		total += 1
        
		for node in tree.get_nonterminals():
			subnodes = set([x.name for x in node.get_terminals()])
			if subnodes in [pattern0, pattern1, pattern2, pattern3]:
				signature.append(subnodes)
					
		#print signature		
		if pattern1 in signature and pattern0 in signature:
			case1 +=1
		elif pattern2 in signature and pattern0 in signature:
			case2 +=1
		elif pattern3 in signature and pattern0 in signature:
			case3 +=1
	
	total = case1+case2+case3
	ratio1 = float(case1)/total*100
	ratio2 = float(case2)/total*100
	ratio3 = float(case3)/total*100		
	return total, ratio1, ratio2, ratio3
	
					
				
def main():
	parser = argparse.ArgumentParser(prog="parse_internode.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--infile', \
	help=("the concatenated set of the RAxML bipartition trees"), required=True)
	parser.add_argument('-b', '--bootstrap_cutoff', \
	help=("a designated bootstrap cutoff"), default=80)
	parser.add_argument('-s', '--species', nargs='*', required=True)
	args = parser.parse_args()
	
	infile = args.infile
	cutoff = int(args.bootstrap_cutoff)
	species = list(args.species)
	
	filterTree(infile, cutoff)
	infile = open("trees_bootstrap_"+str(cutoff), "r")
	total, case1, case2, case3 = parse_tre(infile, species)
	print "%d, %.2f%%, %.2f%%, %.2f%%" % (total, case1, case2, case3)
	infile.close()
	
	

if __name__ == "__main__":
	main()


