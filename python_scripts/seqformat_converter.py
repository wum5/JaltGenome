#!/usr/bin/env python
"""
This script is to convert the sequence alignment format between .fas/.phy/.nex
"""

from Bio import AlignIO
from Bio import Alphabet
import argparse, os


def seq_format_converter(inDIR, outDIR, inputFormat, outputFormat):
	for i in os.listdir(inDIR):
		if not i.endswith(inputFormat): continue	
		clusterID = i.split(inputFormat)[0]
		if inputFormat == '.fas' and outputFormat == '.phy':
			AlignIO.convert(inDIR+i, "fasta", outDIR+clusterID+".phy", "phylip-sequential", alphabet=Alphabet.generic_dna)
		elif inputFormat == '.phy' and outputFormat == '.fas':
			AlignIO.convert(inDIR+i, "phylip-sequential", outDIR+clusterID+".fas", "fasta", alphabet=Alphabet.generic_dna)
		elif inputFormat == '.phy' and outputFormat == '.nex':
			AlignIO.convert(inDIR+i, "phylip-sequential", outDIR+clusterID+".nex", "nexus", alphabet=Alphabet.generic_dna)	
		elif inputFormat == '.fas' and outputFormat == '.nex':
			AlignIO.convert(inDIR+i, "fasta", outDIR+clusterID+".nex", "nexus", alphabet=Alphabet.generic_dna)	


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--inDIR', required=True)
	parser.add_argument('-o', '--outDIR', required=True)
	parser.add_argument('-f1', '--inputFormat', required=True)
	parser.add_argument('-f2', '--outputFormat', required=True)
	args = parser.parse_args()

	inDIR = args.inDIR+"/"
	outDIR = args.outDIR+"/"
	inputFormat = args.inputFormat
	outputFormat = args.outputFormat
	
	seq_format_converter(inDIR, outDIR, inputFormat, outputFormat)


if __name__ == "__main__":
	main()
