#!/usr/bin/env python
"""
## this module is to treat post-aligned coding sequences, to
1) check stop codon within sequence
2) remove sites with indels, missing values or with low coverages
3) sequence length lower than 200
"""


import os, sys, argparse
from Bio import SeqIO


class OrfAln:
	
	def __init__(self, fileFASTA, seqname):
		self.fileFASTA = [f for f in sorted(fileFASTA, key=lambda x : x.id)]
		self.number = len(self.fileFASTA)
		for record in self.fileFASTA:	
			if record.id == seqname:
				self.number -= -1
	
	def StopCodons(self):
		for record in self.fileFASTA:
			stop = 0
			for i in xrange(0, len(record.seq)-3, 3):
				codon = str(record.seq)[i:i+3]
				if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
					raise Exception("stop codon in sequence")
	
	def site_depth(self, seqname):
		codon_cov = []
		num = 0 
		nucleotides = ['A', 'T', 'C', 'G']
		for record in self.fileFASTA:
			if num == 0:
				codon_cov =  [0 for x in xrange(0, len(record.seq)-3, 3)]
			for x in xrange(len(codon_cov)):
				if record.seq[x*3] in nucleotides and \
				record.seq[x*3+1] in nucleotides and \
				record.seq[x*3+2] in nucleotides and record.id != seqname:
					codon_cov[x] += 1
				else:
					pass
			num += 1 
		return codon_cov

	def seq_filter(self, seqname, length_cutoff, depth_cutoff):
		if depth_cutoff == 0:
			depth_cutoff = self.number
		self.StopCodons()
		seq_dict = {}
		seqLen = 1000000
		codon_cov = self.site_depth(seqname)
		for record in self.fileFASTA:
			if record.id == seqname:
				continue
			seq_dict[record.id] = ''
			for x in xrange(len(codon_cov)):
				if codon_cov[x] >= depth_cutoff:
					seq_dict[record.id] += str(record.seq[x*3:x*3+3])

		for key in seq_dict:
			currLen = seq_dict[key].count('A')+seq_dict[key].count('C')\
			+seq_dict[key].count('G')+seq_dict[key].count('T')
			seqLen = min(currLen, seqLen)		
		if seqLen > length_cutoff:
			return seq_dict, seqLen
		else:
			return None, None

			
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_folder', required=True)
	parser.add_argument('-o', '--output_folder', required=True)
	parser.add_argument('-s', '--sequence_delete', default='')
	parser.add_argument('-l', '--length_cutoff', default=200)
	parser.add_argument('-d', '--depth_cutoff', default=0)
	args = parser.parse_args()
	
	inDIR = args.input_folder+"/"
	outDIR = args.output_folder+"/"
	seqname = args.sequence_delete
	lencut = int(args.length_cutoff)
	depcut = int(args.depth_cutoff)
		
	for i in os.listdir(inDIR):
		if i[-4:] != ".fas": continue
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")		

		myseq = OrfAln(infile, seqname)
		seqdict, seqlen = myseq.seq_filter(seqname, lencut, depcut)

		if seqdict == None:
			continue
		print i
		
		aln_stat = open(outDIR+"alignment_stat.txt", "a")
		aln_stat.write("%s\t%d\n" % (i, seqlen))
		outfile = open(outDIR+i,"w")
		
		for key in seqdict:
			outfile.write('>'+key+'\n')
			outfile.write(seqdict[key]+'\n')

		handle.close()
		outfile.close()


			
