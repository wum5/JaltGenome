import os, argparse

"""
This python script is used to parse the number of ancestral segragating variants among major
specie subclades of Jaltomata. The input file should be a mvf file, in which species should
be sorted in order, such as 1-3 for group1, 4-6 for group2... 
"""		

def transform(altbases,species_num):
	changed_bases = ''
	if len(altbases) == species_num+1:
		changed_bases = altbases[1:]
	elif len(altbases) == 1:
		changed_bases = altbases * species_num
	elif '+' in altbases:
		if len(altbases.split('+')[0]) == 1:
			index = int(altbases.split('+')[1][1:])
			alt = altbases.split('+')[1][0]
			changed_bases = altbases[0]*(index-1)+alt+altbases[0]*(species_num-index)
		elif len(altbases.split('+')[0]) == 2:
			index = int(altbases.split('+')[1][1:])
			alt = altbases.split('+')[1][0]
			changed_bases = altbases[1]*(index-1)+alt+altbases[1]*(species_num-index)
	return changed_bases


def snp_call(altbases):
	
	variants = []
	
	for x in xrange(0,len(altbases)):
		if altbases[x] in ['A','T','C','G']:
			variants.append(altbases[x])
		elif altbases[x] == 'Y':
			variants.append('C')
			variants.append('T')
		elif altbases[x] == 'R':
			variants.append('A')
			variants.append('G')
		elif altbases[x] == 'W':
			variants.append('A')
			variants.append('T')
		elif altbases[x] == 'S':
			variants.append('G')
			variants.append('C')
		elif altbases[x] == 'K':
			variants.append('T')
			variants.append('G')
		elif altbases[x] == 'M':
			variants.append('C')
			variants.append('A')
		
	variants = set(variants)
	if len(variants) == 2:
		return variants


if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--mvf_file', required=True)
	parser.add_argument('-t', '--test', choices=["shared_snp", "species_hetero", "shared_hetero"], required=True)	
	args = parser.parse_args()
	infile = open(args.mvf_file, "r")

	# index boundary for lineages in each sub-clade
	gp1_index = 3   #1-3
	gp2_index = 4	#4
	gp3_index = 6	#5-6
	gp4_index = 14	#7-14
	species_num = 14

	if args.test == "shared_snp":    
	## calculate the number of shared ancestral sorting alleles (variant sites) among major species subclades
			
		outfile = open("ancestral_snps.txt", "w")
		snp_gp1, snp_gp2, snp_gp3, snp_gp4 = 0, 0, 0, 0
		snp_gp12, snp_gp13, snp_gp14, snp_gp23, snp_gp24, snp_gp34 = 0, 0, 0, 0, 0, 0
		snp_gp123, snp_gp124, snp_gp134, snp_gp234 = 0, 0, 0, 0
		snp_gp1234 = 0
		
		for line in infile:
			if ':' not in line: continue
			line = line.rstrip()
			altbases = line.split()[1]
			allbases = transform(altbases,species_num)
			if sum(map(allbases.count, ['A','T','C','G','Y','R','W','S','K','M'])) < species_num: continue
										
			# black-fruited lineages
			group1 = allbases[:gp1_index]
			gp1_var = snp_call(group1)
			if gp1_var: snp_gp1 += 1	

			# red-fruited lineages		
			group2 = allbases[gp1_index]
			gp2_var = snp_call(group2)
			if gp2_var: snp_gp2 += 1
	
			# green-fruited lineages
			group3 = allbases[gp2_index:gp3_index]
			gp3_var = snp_call(group3)
			if gp3_var: snp_gp3 += 1
		
			# orange-fruited lineages
			group4 = allbases[gp3_index:gp4_index]
			gp4_var = snp_call(group4)
			if gp4_var: snp_gp4 += 1
		
			if gp1_var and gp2_var and gp1_var == gp2_var:
				snp_gp12 += 1
			if gp1_var and gp3_var and gp1_var == gp3_var:
				snp_gp13 += 1
			if gp1_var and gp4_var and gp1_var == gp4_var:
				snp_gp14 += 1
			if gp2_var and gp3_var and gp2_var == gp3_var:
				snp_gp23 += 1
			if gp2_var and gp4_var and gp2_var == gp4_var:
				snp_gp24 += 1
			if gp3_var and gp4_var and gp3_var == gp4_var:
					snp_gp34 += 1
				
			if gp1_var and gp2_var and gp3_var and gp1_var == gp2_var == gp3_var:
				snp_gp123 += 1
			if gp1_var and gp2_var and gp4_var and gp1_var == gp2_var == gp4_var:
				snp_gp124 += 1	
			if gp1_var and gp3_var and gp4_var and gp1_var == gp3_var == gp4_var:
				snp_gp134 += 1	
			if gp2_var and gp3_var and gp4_var and gp2_var == gp3_var == gp4_var:
				snp_gp234 += 1
				
			if gp1_var and gp2_var and gp3_var and gp4_var and gp1_var == gp2_var == gp3_var == gp4_var:
				snp_gp1234 += 1	

		outfile.write("%d, %d, %d, %d\n" % (snp_gp1, snp_gp2, snp_gp3, snp_gp4))
		outfile.write("%d, %d, %d, %d, %d, %d\n" % (snp_gp12, snp_gp13, snp_gp14, snp_gp23, snp_gp24, snp_gp34))
		outfile.write("%d, %d, %d, %d\n" % (snp_gp123, snp_gp124, snp_gp134, snp_gp234))
		outfile.write("%d\n" % (snp_gp1234))
	
	
	if args.test == "species_hetero":
	## calculate heterozygous rates in each investigated species

		outfile = open("species_hetero.txt", "w")
		homozygotes = [0.0 for x in range(species_num)]
		heterozygotes = [0.0 for x in range(species_num)]
		hetero_rate = [0.0 for x in range(species_num)]	
		total_sites = [0.0 for x in range(species_num)]	

		for line in infile:
			if ':' not in line: continue
			line = line.rstrip()
			altbases = line.split()[1]
			allbases = transform(altbases,species_num)
			if sum(map(allbases.count, ['A','T','C','G','Y','R','W','S','K','M'])) < species_num: continue

			for x in xrange(0,species_num,1):			
				if allbases[x] in  ['A','T','C','G']:
					homozygotes[x] += 1.0 
				elif allbases[x] in ['Y','R','W','S','K','M']:
					heterozygotes[x] += 1.0 
					
		total_sites = [x+y for x, y in zip(heterozygotes, homozygotes)]
		hetero_rate = [x/y for x, y in zip(heterozygotes, total_sites)]
		
		outfile.write("number of heterozygous sites\t"+"\t".join(str(x) for x in heterozygotes)+"\n")
		outfile.write("number of total sites\t"+"\t".join(str(x) for x in total_sites)+"\n")
		outfile.write("rate of heterozygous sites\t"+"\t".join(str(x) for x in hetero_rate)+"\n")
		

	if args.test == "shared_hetero":
	## calculate the number of heterozygous sites that are also sorting in at least other two more major species subclades

		outfile = open("shared_hetero.txt", "w")
		total_hetero = [0.0 for x in range(species_num)]	
		shared_hetero = [0.0 for x in range(species_num)]	

		for line in infile:
			if ':' not in line: continue
			line = line.rstrip()
			altbases = line.split()[1]
			allbases = transform(altbases,species_num)
			if sum(map(allbases.count, ['A','T','C','G','Y','R','W','S','K','M'])) < species_num: continue
		
			gp1_bases = allbases[:gp1_index]
			gp2_bases = allbases[gp1_index:gp2_index]
			gp3_bases = allbases[gp2_index:gp3_index]
			gp4_bases = allbases[gp3_index:]
			
			for x in xrange(species_num):
				if x < gp1_index and allbases[x] in ['Y','R','W','S','K','M']:
					total_hetero[x] += 1
					variants = snp_call(allbases[x]) 
					if [snp_call(gp2_bases)==variants,snp_call(gp3_bases)==variants,snp_call(gp4_bases)==variants].count(True)>=1:  
						shared_hetero[x] += 1	

				elif gp1_index <= x < gp2_index and allbases[x] in ['Y','R','W','S','K','M']:
					total_hetero[x] += 1
					variants = snp_call(allbases[x]) 
					if [snp_call(gp1_bases)==variants,snp_call(gp3_bases)==variants,snp_call(gp4_bases)==variants].count(True)>=1:
						shared_hetero[x] += 1				

				elif gp2_index <= x < gp3_index and allbases[x] in ['Y','R','W','S','K','M']:
					total_hetero[x] += 1
					variants = snp_call(allbases[x]) 
					if [snp_call(gp1_bases)==variants,snp_call(gp2_bases)==variants,snp_call(gp4_bases)==variants].count(True)>=1:
						shared_hetero[x] += 1			

				elif gp3_index <= x and allbases[x] in ['Y','R','W','S','K','M']:
					total_hetero[x] += 1
					variants = snp_call(allbases[x]) 
					if [snp_call(gp1_bases)==variants,snp_call(gp2_bases)==variants,snp_call(gp3_bases)==variants].count(True)>=1:
						shared_hetero[x] += 1	

		shared_rate = [x/y for x, y in zip(shared_hetero, total_hetero)]

		outfile.write("number of shared heterozygous sites\t"+"\t".join(str(x) for x in shared_hetero)+"\n")
		outfile.write("number of total heterozygous sites\t"+"\t".join(str(x) for x in total_hetero)+"\n")
		outfile.write("rate of shared heterozygous sites\t"+"\t".join(str(x) for x in shared_rate)+"\n")
		

	infile.close()
	outfile.close()	