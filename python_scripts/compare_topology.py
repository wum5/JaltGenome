from Bio import Phylo
from cStringIO import StringIO
import sys, os


def get_nodes(tree):
	
	node_list = []	
	for node in tree.get_nonterminals():
		subnodes = set([x.name for x in node.get_terminals()])
		node_list.append(subnodes)
	
	return node_list
	
	
def comparison(tree1,tree2):
	
	nodelist1 = get_nodes(tree1)
	nodelist2 = get_nodes(tree2)	
	same_topology = True
		
	for x in xrange(len(nodelist1)):
		if nodelist1[x] not in nodelist2:
			same_topology = False
			break
	
	return same_topology
		
		

if __name__ == "__main__":
	
	if len(sys.argv) != 3:
		print "usage: python compare_topology.py geneTrees speciesTree"
		sys.exit()
	
	inDIR = sys.argv[1]		
	speciestree = open(sys.argv[2], "r")	
	outfile = open("gene_tree_stat.txt", "w")
	handle = StringIO(speciestree.readline().rstrip())
	speciesTree = Phylo.read(handle, "newick")
	
	agree_num = 0
	diff_num = 0
	gene_tre_list = []
	
	for file in os.listdir(inDIR):
		with open (inDIR+'/'+file, "r") as infile:
			line = infile.readline()
			line = line.rstrip()
			handle = StringIO(line)
			tree = Phylo.read(handle, "newick")
		
			if comparison(tree, speciesTree):
				agree_num += 1
				print file.split('RAxML_bestTree.')[1]
		
			for x in xrange(len(gene_tre_list)):
				if comparison(tree, speciesTree):
					diff_num -= 1
					break
		
		diff_num += 1
		gene_tre_list.append(tree)
		
	
	print agree_num, diff_num	
	outfile.write("number of genes trees agree with species tree: %d\n" % (agree_num))
	outfile.write("number of different gene tree topologies: %d\n" % (diff_num))
	
	speciestree.close()
	outfile.close()
	