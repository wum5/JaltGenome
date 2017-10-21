#!/usr/bin/env python
"""
This script is to investiagate the unsolved internode in my investiated speices tree,
python parse_internode.py -i genetrees.tre -s C.annuumm J.sinuosa S.lycopersicum S.tuberosum
here species id are ordered as A,B,C1,C2. 
"""


from Bio import Phylo
from cStringIO import StringIO
import argparse


def parse_tre(infile, species):
	"""
	in this case, we are looking for the frequency of internodes support 1. (A,B,(C1,C2));
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
	parser.add_argument('-i', '--infile', required=True)
	parser.add_argument('-s', '--species', nargs='*', required=True)
	args = parser.parse_args()

	infile = open(args.infile, "r")
	species = list(args.species)
	#print species

	total, case1, case2, case3 = parse_tre(infile,species)
	print "%d, %.2f%%, %.2f%%, %.2f%%" % (total, case1, case2, case3)
	infile.close()


if __name__ == "__main__":
	main()
	
	
