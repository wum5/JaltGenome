#!/usr/bin/env python
"""
This script is to concatenate sequence alignments from individual genes into
a super gene sequence alignment. The gene sequence alignments come from Prank, 
ending with ".fas".
"""

import argparse, os
from Bio import SeqIO


def rename(seqid):
	"""
	rename the sequence id in my particular project
	"""
	if seqid.startswith('CA'):
		name = 'C.annuumm'
	elif seqid.startswith('AT'):
		name = 'A.thaliatha'
	elif seqid.startswith('Jalsin'): 
		name = 'J.sinuosa'
	elif seqid.startswith('Solyc'): 
		name = 'S.lycopersicum'
	elif seqid.startswith('PGSC'): 
		name = 'S.tuberosum'
	elif seqid.startswith('NIAT'): 
		name = 'N.attenuata'
	elif seqid.startswith('Peinf'): 
		name = 'P.inflata'
	elif seqid.startswith('Peaxi'): 
		name = 'P.axillaris'
	return name


def filter(inDIR, outDIR):
	for i in os.listdir(inDIR):
		print i
		with open(outDIR+'/'+i, "w") as outfile:
			for record in SeqIO.parse(inDIR+'/'+i, "fasta"):
				geneid = rename(str(record.id))
				outfile.write('>'+geneid+'\n')
				outfile.write(str(record.seq)+'\n')

		
def main():
	parser = argparse.ArgumentParser(prog="rename_seq.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--inDIR', required=True)
	parser.add_argument('-o', '--outDIR', required=True)
	args = parser.parse_args()

	indir = args.inDIR
	outdir = args.outDIR
	filter(indir, outdir)
	

if __name__ == "__main__":
	main()
	
	