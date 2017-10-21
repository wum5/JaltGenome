#!/usr/bin/env python

"""identify the potential contaminants in genome assembly
based on the best BLAST-hitted reference, the input file
should be the BLASTN output file (-outfmt 6)"""

import os, sys, argparse

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "usage: python combine_blastout.py infile"
		sys.exit()

	infile = open(sys.argv[1], "r")
	gi_out = open("blast_gi2.txt", "w")
	query = ''

	for line in infile:
		line = line.rstrip()
		scaff = line.split('\t')[0]
		ref = line.split('\t')[1]
		score = float(line.split('\t')[11])

		# the first query sequence
		if query == '':
			query = scaff
			max_score = cur_score = score 
			best_hit = cur_hit = ref 

		# the same query
		elif scaff == query:
			# the same reference
			if ref == cur_hit:
				cur_score += score
			# different reference
			else:
				if cur_score > max_score:
					max_score = cur_score
					best_hit = cur_hit
				cur_score = score 
				cur_hit = ref 

		# the next query and return the output of last query
		elif scaff != query:
			gi = best_hit.split('|')[1]
			gi_out.write(query+'\t'+gi+'\t'+str(max_score)+'\n')
			query = scaff
			max_score = cur_score = score 
			best_hit = cur_hit = ref 


	# return the output of the last one query
	gi = best_hit.split('|')[1]
	gi_out.write(query+'\t'+gi+'\t'+str(max_score)+'\n')
	
	infile.close()
	gi_out.close()			

	
