#!/usr/bin/env python
"""
This script is to run branch-site+MNM model on multiple sequence alignments 
"""

from scipy import stats
import argparse, os


def parse_LNL(BS_out):
	
	with open(BS_out) as infile:
		for line in infile:
			line = line.rstrip()
			if '-' in line and '{' in line:
				LNL = float('-'+(line.split(',')[0]).split('-')[1])
	
	return LNL




def main():
	parser = argparse.ArgumentParser(prog="paml_BS_MNM.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--inDIR', required=True)
	parser.add_argument('-d', '--MNMDIR', required=True)
	parser.add_argument('-o', '--output', required=True)
	args = parser.parse_args()

	outfile = open(args.output, "w")
	outfile.write("geneID\t2difL\tpval\n")
	inDIR = args.inDIR
	MNMDIR = args.MNMDIR
	
	for i in os.listdir(inDIR):
		if not i.startswith("Jalsin"): continue
		cmd = "cp "+inDIR+"/"+i+" "+MNMDIR+"/tempfile"
		os.system(cmd)
		cmd = "HYPHYMPI "+MNMDIR+"/BranchSites_delta_null.bf"
		os.system(cmd)
		cmd = "HYPHYMPI "+MNMDIR+"/BranchSites_delta_alt.bf"
		os.system(cmd)
		
		LNL_null = parse_LNL(MNMDIR+"/tempfile.null.out")
		LNL_alt = parse_LNL(MNMDIR+"/tempfile.alt.out")
		chisq = 2*(LNL_alt-LNL_null)
		pval = stats.chi2.pdf(chisq, 1)
		
		outfile.write("%s\t%.3f\t%.6f\n" % (i, chisq, pval))
		
	outfile.close()
		 


if __name__ == "__main__":
	main()
	