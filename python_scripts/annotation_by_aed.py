#!/usr/bin/env python
"""\n This program is to extract high-confidence gene annotions from MAKER2 \
output based on the annotation edit distance (AED). AED = 1 - AC, where \
AC = (SN + SP) / 2. SN = TP / (TP + FN); SP = TP / (TP + FP); \
eAED is in terms of exons ,rather than bases. """


import argparse


def extract_value(string):
	"""\
	Extract value after '='
	"""
	value = string.split('=')[1]
	return value


def parse_maker_gff3(infile, index, AED_cutoff, eAED_cutoff, gff_out):
	"""\
	Extract high-confidence gene models from MAKER2 output using the AED/eAED cutoffs.
	"""
	if gff_out and AED_cutoff < 1:
		outfile = open(index+'.aed.'+str(AED_cutoff)+'.gff3', "w")
		outfile.write('##gff-version 3\n')
	
	elif gff_out and eAED_cutoff < 1:
		outfile = open(index+'.eaed.'+str(eAED_cutoff)+'.gff3', "w")
		outfile.write('##gff-version 3\n')
		
	# initialization
	output = []
	gene_num, mRNA_num = 0, 0
	
	for line in infile:
		if 'maker' not in line: continue
		line = line.rstrip()
		feature = line.split('\t')[2]
		annotation = line.split('\t')[8]
		
		if feature == 'gene':
			# print the last qualified gene model
			if len(output) > 1:
				gene_num += 1
				if gff_out:
					for element in output:
						outfile.write(element+'\n')	
					
			# re-initialization
			output = []
			output.append(line)
			mrna_ids = []
			
		elif feature == 'mRNA':
			AED = float(extract_value(annotation.split(';')[3]))
			eAED = float(extract_value(annotation.split(';')[4]))
			if AED <= AED_cutoff and eAED <= eAED_cutoff:
				tmp = extract_value(annotation.split(';')[0])
				mrna_ids.append(tmp)
				mRNA_num += 1
				output.append(line)
				
		elif feature in ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']:
			tmp = extract_value(annotation.split(';')[1])
			parents = tmp.split(",")
			for item in parents:
				if item in mrna_ids:
					output.append(line)
									
		else: 
			pass
		
	# print the lastest qualified gene model
	if len(output) > 1:
		gene_num += 1
		if gff_out: 
			for element in output:
				outfile.write(element+'\n')		
	
	if gff_out:
		outfile.close()		
			
	if AED_cutoff < 1:
		print "AED cutoff: %f, number of genes: %d, number of transcripts: %d\n" \
		% (AED_cutoff, gene_num, mRNA_num)

	if eAED_cutoff < 1:
		print "eAED cutoff: %f, number of genes: %d, number of transcripts: %d\n" \
		% (eAED_cutoff, gene_num, mRNA_num)


def main():
	parser = argparse.ArgumentParser(prog="annotation_by_AED.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--gff', required=True)
	parser.add_argument('--AED', default=1)
	parser.add_argument('--eAED', default=1)
	parser.add_argument('--out', action="store_true")
	args = parser.parse_args()
    
	gff_file = args.gff
	index = gff_file.split('.')[0]
	AED_cutoff = float(args.AED)
	eAED_cutoff = float(args.eAED)
	gff_out = args.out

	infile = open(gff_file, "r")
	parse_maker_gff3(infile, index, AED_cutoff, eAED_cutoff, gff_out)
	infile.close()

	

if __name__ == "__main__":
	main()
	
	