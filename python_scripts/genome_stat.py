#!/usr/bin/env python
"""calculate genome contents, such as total bases, GC
contents, genome coverage, N50 and N90.
generate an output fasta file with a minimum  
requirement for contig length."""


import os, sys, argparse
nucleotides = ['A', 'T', 'C', 'G']


class classFASTA:
	"""
	Parse fasta-formated sequences into lists of sequences and ids
	"""
	def __init__(self, fileFASTA):
		self.fileFASTA = fileFASTA

	def readFASTA(self):
		"""Checks for fasta by file extension"""
		file_lower = self.fileFASTA.lower()
		"""Check for three most common fasta file extensions"""
		if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
		file_lower.endswith('.fasta') or file_lower.endswith('.fna'):
			with open(self.fileFASTA, "r") as f:
				return self.ParseFASTA(f)
		else:
			print "Not in FASTA format."

	def ParseFASTA(self, fileFASTA):
		"""Gets the sequence name and sequence from a FASTA formatted file"""
		fasta_list=[]
		for line in fileFASTA:
			if line[0] == '>':
				try:
					fasta_list.append(current_dna)
				#pass if an error comes up
				except UnboundLocalError:
					#print "Inproper file format."
					pass
				current_dna = [line.lstrip('>').rstrip('\n'),'']
			else:
				current_dna[1] += "".join(line.split())
		fasta_list.append(current_dna)
		"""Returns fasa as nested list, containing line identifier \
		and sequence"""
		return fasta_list


def genome_contents(genome_seq):
	"""
	Calculate GC contents
	"""
	total = len(genome_seq)
	G = genome_seq.count('G')
	C = genome_seq.count('C')
	GC = str(round(float(G+C)/total*100,2))+'%'
	return GC


def N50_estimation(contigs, genome_size):
	"""
	Calculate N50, N90, minimum and maximum contigs, etc
	"""
	Contigs_Len = [ len(x) for x in contigs ]
	Contigs_Len.sort(reverse=True)
	Total = sum(Contigs_Len)
	Len = len(Contigs_Len)
	Max = Contigs_Len[0]
	Min = Contigs_Len[-1]
	Avg = sum(Contigs_Len)/Len
	Cov = str(round(float(Total)/genome_size*100,2))+'%'
	accumulate = 0
	
	for x in range(len(Contigs_Len)):
		accumulate += Contigs_Len[x]
		if accumulate > Total/2:
			N50 = Contigs_Len[x]
			stop = x
			L50 = x
			break
			
	for x in range(stop,len(Contigs_Len)):
		accumulate += Contigs_Len[x]
		if accumulate > Total*0.9:
			N90 = Contigs_Len[x]
			L90 = x
			break

	return Total, Cov, Len, N50, L50, N90, L90, Max, Min, Avg
	

def filter_contigs(contig_ids, contig_seqs, minLen):
	"""
	Filter out contigs smaller than specified length
	"""
	outfile1 = open("scaffold_min_"+str(minLen)+".fa", "w")
	outfile2 = open("scaffold_max_"+str(minLen)+".fa", "w")
	for x in xrange(len(contig_ids)):
		if len(contig_seqs[x]) < minLen:		
			outfile2.write('>'+contig_ids[x]+'\n')
			outfile2.write(contig_seqs[x]+'\n')
		else:		
			outfile1.write('>'+contig_ids[x]+'\n')
			outfile1.write(contig_seqs[x]+'\n')
	outfile1.close()
	outfile2.close()


def main():
	parser = argparse.ArgumentParser(prog="genome_stat.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--fasta_file', required=True)
	parser.add_argument('-gs', '--genome_size', required=True)
	parser.add_argument('-m', '--minLen', default=0)
	args = parser.parse_args()

	infile = classFASTA(args.fasta_file)
	gs = int(args.genome_size)
	minLen = int(args.minLen)
	
	######### This command returns the nested list containing sequence names
	######### and sequences to a flat list containing only sequences or ids
	contig_ids = [ x[0] for x in infile.readFASTA() ]
	contig_seqs = [ x[1] for x in infile.readFASTA() ]
	contig_seqs = [ x.upper() for x in contig_seqs ]
	
	if minLen != 0:
		filter_contigs(contig_ids, contig_seqs, minLen)

	######### concatenated the sequences
	genome = ''.join(contig_seqs[:])

	GC = genome_contents(genome)
	Total, Cov, Num, N50, L50, N90, L90, Max, Min, Avg = N50_estimation(contig_seqs,gs)
	
	print "Total bases: %d, Genome coverage: %s, GC content: %s" % (Total, Cov, GC)
	print "#Contigs: %d, N50: %d bps, L50: %d, N90: %d bps, L90: %d, \nMax: %d bps, Min: %d bps, Mean: %d bps" \
	% (Num, N50, L50, N90, L90, Max, Min, Avg)
	


if __name__ == "__main__":
	main()
	
