#!/usr/bin/env python
"""
This script is to mask pooly aligned regions using 15-bp sliding 
window. For example to mask regions with average pairwise divergent 
sites more than 5 (excluding A.thaliatha): 
python SlidingWindows.py -i from_prank/ -o trimed_align/ -e A.thaliatha -d 5 
"""


import os, argparse
from Bio import SeqIO


nucleotides = ['A', 'T', 'C', 'G']
dif = lambda seq1, seq2: sum(seq1[i] != seq2[i] for i in \
xrange(len(seq1)) if seq1[i] in nucleotides and seq2[i] in nucleotides)

def window_check(windows, **kw):
	maxEdit = []
	for x in xrange(len(windows)):
		num, allEdit = 0, 0.0
		for y in xrange(len(windows)):
			allEdit += dif(windows[x], windows[y])
			num += 1
		maxEdit.append(allEdit/num)
	return max(maxEdit)

def remove_by_mask(seq_id, seq, mask_cutoff=0.1):
	mask = False
	masked_region = seq.count('N')
	masked_ratio = round(float(masked_region)/len(seq),4)
	if masked_ratio > mask_cutoff:
		mask = True
		print seq_id+'\t'+str(masked_ratio*100)+'%\t(deleted)'
		return mask
	else:
		print seq_id+'\t'+str(masked_ratio*100)+'%'


def mask_process(inDIR, outDIR, outgroup, num_div, win_size):
	for i in os.listdir(inDIR):
		if i[-4:] != ".fas": continue
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")
		windows = []
		num, snp = 0, 0
		masked_region = win_size*'N'
	
		for record in infile:
			if outgroup not in str(record.id):
				if num == 0:
					for x in xrange(len(record.seq)-win_size):
						windows.append([str(record.seq[x:x+win_size])])
				else:
					for x in xrange(len(record.seq)-win_size):
						windows[x].append(str(record.seq[x:x+win_size]))
				num += 1 
											
		mask_win = []
		for x in xrange(len(windows)-win_size):
			score = window_check(windows[x])
			if score > num_div:	
				mask_win.append(x)
			
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")
					
		with open(outDIR+i,"w") as outfile:
			for record in infile:
				outfile.write('>' + record.id + '\n')
				seq = record.seq
				for x in xrange(len(mask_win)):
					start = mask_win[x]
					seq = seq[:start]+masked_region+seq[start+win_size:]

				outfile.write(str(seq)+'\n')
		
		if remove_by_mask(i,seq) == True:
			os.remove(outDIR+i)
		
		handle.close()




def main():		
	parser = argparse.ArgumentParser(prog="SlindingWindows", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--inDIR', help=("the input directory"), required=True)
	parser.add_argument('-o', '--outDIR', help=("the output directory"), required=True)
	parser.add_argument('-e', '--species_excluded', help=("the species not included in trimming"), default='--')
	parser.add_argument('-d', '--num_div', help=("number of divergent sites"), default=5)
	parser.add_argument('-w', '--window_size', help=("size of sliding window"), default=15)
	args = parser.parse_args()

	inDIR = args.inDIR
	outDIR = args.outDIR
	outgroup = args.species_excluded
	num_div = int(args.num_div)
	win_size = int(args.window_size)
		
	mask_process(inDIR, outDIR, outgroup, num_div, win_size)



if __name__ == "__main__":
	main()
