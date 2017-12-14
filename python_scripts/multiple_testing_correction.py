#!/usr/bin/env python
"""
This script is to calculate the adjust p-value using multiple testing method based on
a list of p values given
"""

import pandas as pd 
import argparse


def multiple_testing_correction(pvalues, correction_type="FDR"):
    """
    Consistent with R - print
    correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05,
                                          0.069, 0.07, 0.071, 0.09, 0.1])	
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    sample_size = pvalues.shape[0]
    qvalues = empty(sample_size)
    if correction_type == "Bonferroni":
        # Bonferroni correction
        qvalues = sample_size * pvalues
    elif correction_type == "Bonferroni-Holm":
        # Bonferroni-Holm correction
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            qvalues[i] = (sample_size-rank) * pvalue
    elif correction_type == "FDR":
        # Benjamini-Hochberg, AKA - FDR test
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = sample_size - i
            pvalue, index = vals
            new_values.append((sample_size/rank) * pvalue)
        for i in range(0, int(sample_size)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            qvalues[index] = new_values[i]
    return qvalues
    
 
def function_assign(funct_file, gene_list):
	'''
	function file should have two columns:
	Solyc00g005040.3.1	Potassium channel
	Solyc00g005096.1.1	RWP-RK domain
	'''
	geneid, functL, output = [], [], []
	with open(funct_file, "r") as infile:
		for line in infile:
			if 'Jal' not in line: continue
			line = line.rstrip()
			geneid.append(line.split('\t')[0])
			try:
				functL.append(line.split('\t')[1])
			except:
				functL.append('unknown')
	
	for item in gene_list:
		for x in xrange(len(geneid)):
			if item.split('.')[0] in geneid[x]:
				output.append(functL[x])
	
	return output			

   
    
def main():
	'''
	here is the function to parse the p-value from my custom paml script
	'''
	parser = argparse.ArgumentParser(prog="multiple_testing_correction.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--infile', required=True)
	parser.add_argument('-o', '--outfile', required=True)
	parser.add_argument('-f', '--function_file', required=True)
	args = parser.parse_args()
	
	df = pd.read_csv(args.infile, sep='\t')
	df['fdr'] = multiple_testing_correction(df['pval'])
	adj_df = df.sort_values(by=['fdr'])
	
	adj_df['function'] = function_assign(args.function_file, adj_df['geneID'])
	adj_df.to_csv(args.outfile, sep='\t', index=False)
	

if __name__ == "__main__":
	main()


    