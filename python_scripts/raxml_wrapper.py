#!/usr/bin/env python
"""
Change the raxml command (RAXML_CMD) depend on which flavor of raxml and where the
executable is in your machine

Input: a dir of cleaned alignments in fasta format and end with ".fas"
Output: tree names are clusterID.raxml
"""

ALIGNMENT_FILE_ENDING = ".fas"
RAXML_CMD = "raxmlHPC-PTHREADS"
#RAXML_CMD = "raxmlHPC-PTHREADS-SSE3-icc"
import os, argparse


if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog="SlindingWindows", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--inDIR', help=("the input directory"), required=True)
	parser.add_argument('--outDIR', help=("the output directory"), required=True)
	parser.add_argument('--outgroup', help=("the outgroup species"), default=None)
	parser.add_argument('--cpus', help=("number of cpus for parrallel"), default=2)
	parser.add_argument('--type', help=("seqeuce type (DNA/aa)"), choices=['DNA', 'aa'])
	parser.add_argument('--bootstrap', help=("bootstrap"), action="store_true")
	args = parser.parse_args()
	
		
	inDIR = args.inDIR+"/"
	outDIR = args.outDIR+"/"
	num_cores = str(args.cpus)
	outgroup = args.outgroup
	bootstrap = args.bootstrap
	
	if args.type == "aa": model = "PROTCATWAG"
	elif args.type == "DNA": model = "GTRCAT" 	
	else:
		print "Input data type: DNA or aa"
		sys.exit()
	
	os.chdir(inDIR)
	cmd="mkdir temp"
	os.system(cmd)
	
	l = len(ALIGNMENT_FILE_ENDING)
	for j in os.listdir(inDIR):
		if j[-l:] != ALIGNMENT_FILE_ENDING:
			continue
		clusterID = j.split(".")[0]
		print clusterID
		if outgroup:
			cmd = RAXML_CMD+" -T "+num_cores+" -p 12345 -s "+j+" -n "+j+" -m "+model+" -o "+outgroup
		else:
			cmd = RAXML_CMD+" -T "+num_cores+" -p 12345 -s "+j+" -n "+j+" -m "+model
		if bootstrap:
			cmd += " -f a -x 12345 -N 100"
		else:
			pass
		#print cmd
		try:
			os.system(cmd)
		except:
			raise Exception("Error in phylogeny reconstruction")
	
	cmd="mv *bestTree* "+outDIR
	os.system(cmd)

	if bootstrap:
		cmd="mv *bootstrap* "+outDIR
		os.system(cmd)		

	cmd="mv RAxML* temp"
	os.system(cmd)
	
