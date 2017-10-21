import sys, os

"""
Sub-package tons of input files in a main directory into multiple subdirectories
based on the maximum number of input files per subdirectories
"""


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "usage: python directory_subpackage.py inDIR max_num_files FileEnding"
		sys.exit()
		
	inDIR = sys.argv[1]+"/"
	max_num_file = int(sys.argv[2])
	FileEnding = sys.argv[3]
	alignment_num = 0
	subdir_tag = 1
	os.system("mkdir "+inDIR+"sub_folder_"+str(subdir_tag))
	subDIR = "sub_folder_"+str(subdir_tag)
	
	for i in os.listdir(inDIR):
		if alignment_num < max_num_file and FileEnding in i:  #change the file_ending name if needed
			os.system("mv "+inDIR+i+" "+inDIR+subDIR)			
			alignment_num += 1
		elif FileEnding in i:
			alignment_num = 0
			subdir_tag += 1
			os.system("mkdir "+inDIR+"sub_folder_"+str(subdir_tag))
			subDIR = "sub_folder_"+str(subdir_tag)
			os.system("mv "+inDIR+i+" "+inDIR+subDIR)			
			alignment_num += 1
			
