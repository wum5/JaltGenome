#!/usr/bin/env python
"""\n This program is to perform statistics calculation of gene annotations\n"""


import argparse


def gff_stats(gff_file):
	"""
	Obtain the statistics of gene annotations from the updated gff file
	"""
	exon_num, gene_num, exon_size, cds_size, gene_size = 0, 0, 0, 0, 0

	for line in gff_file:
		if 'maker' not in line: continue
		line = line.rstrip()
		feature = line.split('\t')[2]
		start = int(line.split('\t')[3])
		end = int(line.split('\t')[4])
		size = end - start + 1	
		
		if feature == 'mRNA':
			gene_num += 1
			gene_size += size
			
		elif feature == 'exon':
			exon_num += 1
			exon_size += size
		
		elif feature == 'CDS':
			cds_size += size
		
	avg_exon_num = float(exon_num)/gene_num
	avg_exon_size = float(exon_size)/exon_num
	avg_cds_size = float(cds_size)/gene_num
	avg_transcript_size = float(exon_size)/gene_num
	intron_size = float(gene_size-exon_size)/(exon_num-gene_num)

	print "Mean number of exons per gene: %f" % (avg_exon_num)
	print "Mean exon size: %f" % (avg_exon_size)
	print "Mean transcript length: %f" % (avg_transcript_size)
	print "Mean CDS size: %f" % (avg_cds_size)
	print "Mean intron size: %f" % (intron_size)
	
	
def main():
	parser = argparse.ArgumentParser(prog="annotation_by_AED.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--gff', required=True)
	args = parser.parse_args()		
	
	infile = open(args.gff, "r")
	gff_stats(infile)
	infile.close()
	
	
	
if __name__ == "__main__":
	main()
	