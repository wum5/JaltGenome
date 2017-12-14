#!/usr/bin/env python
"""
This script is to run PAML branch-site model on multiple sequence alignments if the 
gene trees are consistent with the species tree
"""

from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML.chi2 import cdf_chi2
import argparse, os


def codeml_file():
	       
	bsfile="""
	noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
	verbose = 1   * 1: detailed output, 0: concise output
	runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
	* 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
	
	seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
	CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
	clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
	aaDist = 0   * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
	model = 2   * models for codons:
	* 0:one, 1:b, 2:2 or more dN/dS ratios for branches
	NSsites = 2   * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
	* 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal
	icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
	Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
	fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
	kappa = 2   * initial or fixed kappa
	fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
	omega = 1   * initial or fixed omega, for codons or codon-based AAs
	getSE = 0       * 0: don't want them, 1: want S.E.s of estimates
	RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
	Small_Diff = .45e-6  * Default value.
	cleandata = 1       * remove sites with ambiguity data (1:yes, 0:no)?
	fix_blength = 0       * 0: ignore, -1: random, 1: initial, 2: fixed
	"""
  	
  	nullfile="""
  	noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
  	verbose = 1   * 1: detailed output, 0: concise output
  	runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
  	* 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
  	
  	seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
  	CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
  	clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
  	aaDist = 0   * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
  	model = 2   * models for codons:
  	* 0:one, 1:b, 2:2 or more dN/dS ratios for branches
  	NSsites = 2   * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
  	* 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal
  	icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
  	Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
  	fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
  	kappa = 2   * initial or fixed kappa
  	fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate
  	omega = 1   * initial or fixed omega, for codons or codon-based AAs
  	getSE = 0       * 0: don't want them, 1: want S.E.s of estimates
  	RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  	Small_Diff = .45e-6  * Default value.
  	cleandata = 1       * remove sites with ambiguity data (1:yes, 0:no)?
  	fix_blength = 0       * 0: ignore, -1: random, 1: initial, 2: fixed
  	"""
		
 	with open("alt_codeml.ctl", "w") as outfile:
 		outfile.write(bsfile)
 
	with open("null_codeml.ctl", "w") as outfile:
 		outfile.write(nullfile)	


def run_paml(infile, treefile):
	cml = codeml.Codeml(alignment = infile, tree = treefile, \
	out_file = "results.out", working_dir = "./temp")
	cml.read_ctl_file("alt_codeml.ctl")
	#cml.print_options()
	fore = cml.run()
	fore_dat = ((fore.get('NSsites')).get(2))
	fore_lnL = fore_dat.get('lnL')
	fore_omegas = (fore_dat.get('parameters')).get('site classes')
	fore_omega = ((fore_omegas.get(2)).get('branch types')).get('foreground')
	cml.read_ctl_file("null_codeml.ctl")
	#cml.print_options()
	null = cml.run()
	null_lnL = ((null.get('NSsites')).get(2)).get('lnL')

	df = 1
	chisq = 2*(fore_lnL-null_lnL)
	try:	
		pval = cdf_chi2(df, chisq)
	except:
		pval = 1.000
	return fore_omega, chisq, pval


def main():
	parser = argparse.ArgumentParser(prog="paml_BS.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--inDIR', required=True)
	parser.add_argument('-o', '--output', required=True)
	parser.add_argument('-t', '--species_tree', required=True)
	args = parser.parse_args()

	codeml_file()
	outfile = open(args.output, "w")
	outfile.write("geneID\tdnds\t2difL\tpval\n")
	
	for i in os.listdir(args.inDIR):
		if not i.startswith("Jalsin"): continue
		print i
		fore_omega, chisq, pval = run_paml(args.inDIR+'/'+i, args.species_tree)
		outfile.write("%s\t%.3f\t%.3f\t%f\n" % (i,fore_omega,chisq,pval))
		 


if __name__ == "__main__":
	main()
	
	