#!/usr/bin/env python
"""\n This program is prepare the gaf formatted file used for GO annotation,
the input is the mapping file of gene id and GO id. Also, the "go.obo" is required. """


import argparse


def alt_go(goids, gofile):
	'''check alternative go annotation'''
	alt_goids = []
	for x in xrange(len(goids)):
		print goids[x]
		with open(gofile, "r") as f:
			for line in f:		
				line = line.rstrip()
				if "alt_id: "+goids[x] == line:
					#print line
					alt_goids.append(goids[x])
					#print goids[x]
					break
					
	return alt_goids				


def go_search(geneid, goid, database, taxonid, date, gofile):
	dbid = 'db.'+geneid
	with open(gofile, "r") as f:
		for line in f:
			aspect = ''
			line = line.rstrip()
			
			if "id: "+goid == line:
				line = f.next()
				line = line.rstrip()
				if line.startswith('name: '):
					genename = line.split('name: ')[1]
				line = f.next()
				line = line.rstrip()
				if line.startswith('namespace:'):
					if 'biological_process' in line:
						aspect = 'P'
					elif 'molecular_function' in line:
						aspect = 'M'
					elif 'cellular_component' in line:
						aspect = 'C'
					
				output = 'UniProtKB\t'+geneid+'\t0\t'+goid+'\tPMID:0000000\tISO\t0\t'+\
				aspect+'\t'+genename+'\t0\tprotein\t'+'taxon:'+taxonid+'\t'+date+'\t'+database+'\n'
				print output
				return output


	
	
def main():
	parser = argparse.ArgumentParser(prog="prepare_go_gaf.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input_file', help="a file mapping geneid and goid", required=True)
	parser.add_argument('-g', '--go_file', help="the go.obo file", required=True)
	parser.add_argument('-d', '--database', help="where the data come from", required=True)
	parser.add_argument('-x', '--taxon_id', help="the taxon id number", required=True)
	parser.add_argument('-t', '--date', help="the day the database created (e.g. 20170431)", required=True)
	args = parser.parse_args()
    
	gofile = args.go_file
	database = args.database
	taxonid = args.taxon_id	
	date = args.date
	outfile = open('goa_'+database+'.gaf', "w")
	
	# check alternative go annotation
	'''
	infile = open(args.input_file, "r")
	goids = []
	for line in infile:
		line = line.rstrip()
		if 'GO:' in  line:
			goids.append(line.split('\t')[1])
			
	print len(goids)
	goids = list(set(goids))
	print len(goids)
	infile.close()
	
	alt_goids = alt_go(goids, gofile)
	print alt_goids
	'''
	
	infile = open(args.input_file, "r")
	for line in infile:
		line = line.rstrip()
		
		if 'GO:' in  line:
			geneid = line.split('\t')[0]
			goid = line.split('\t')[1]
			output = go_search(geneid, goid, database, taxonid, date, gofile)
			outfile.write(output)
			
	outfile.close()	
	infile.close()	
	
	


if __name__ == "__main__":
	main()	
	
	
	