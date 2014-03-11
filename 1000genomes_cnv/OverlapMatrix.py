import sys

'''
Class that compares CNVs from two subpopulations, 
combines the same size/location cnvs into one, copies all of theem to single file.
Once all cnvs are combined use overlap_final.py to do a 80 percent reciprocal overlap 

author: gerikson
date: March 11 2014
'''

def overlap(line1, line2):
	start1 = int(line1[1])
	start2 = int(line2[1])
	end1 = int(line1[2])
	end2 = int(line2[2])
	#Extracting only the common varriants
	if start1 == start2 and end1 == end2:
		#print "Category 0"
		outfile.write("\t".join(line1) + "\t"+ line2[5] + "\t" + line2[6] + "\t" + line2[7] + "\n")
		#outfile.write("\t".join(line1[:8]) + "\t"+ line2[5] + "\t" + line2[6] + "\t" + line2[7] + "\n")
		return True
	else:
		return False


def indexFile(file_name):
	dict_name = {}

	for line in file_name:
		line = line.strip()
		if not line.startswith("Chr"):
			'''
			temp_line = line.split("\t")
			dict_key = temp_line[0]+ "_" + temp_line[1]+"_"+temp_line[2]
			dict_name[dict_key] = line
			'''
			dict_name[line] = line
	return dict_name


inFilename = str(sys.argv[1])
filenameTwo = str(sys.argv[2])
outfilename = str(sys.argv[3])

inFile = open(inFilename)
fileTwo = open(filenameTwo)
outfile = open(outfilename, 'w')

if __name__ == "__main__":

	file1_index = {}
	file2_index = {}
	#each time we find an overlap we have to remove those keys so we don't add them to the file twice
	file1_values_to_remove = {}
	file2_values_to_remove = {}

	file1_index = indexFile(inFile)
	file2_index = indexFile(fileTwo)


	#outfile.write("Chromosome\tBegin\tEnd\tVarType\tFilter\tCEU_ID\tCEU_Sample_Count\tCEU_Total_Sample\tFIN_ID\tFIN_Sample_Count\tFIN_Total_Sample\tGBR_ID\tGBR_Sample_Count\tGBR_Total_Sample\tIBS_ID\tIBS_Sample_Count\tTSI_Total_Sample\tTSI_ID\tTSI_Sample_Count\tTSI_Total_Sample\n")
	#outfile.write("Chromosome\tBegin\tEnd\tVarType\tFilter\tASW_ID\tASW_Sample_Count\tASW_Total_Sample\tLWK_ID\tLWK_Sample_Count\tLWK_Total_Sample\tYRI_ID\tYRI_Sample_Count\tYRI_Total_Sample\n")
	counter1 = 0 
	for line_index in file1_index.keys():
		counter2 = 0
		counter1 = counter1 + 1
		#print "counter 1: " + str(counter1)
		l1 = file1_index[line_index]
		line1 = l1.split("\t")
		for line2_index in file2_index.keys():
			#counter2 = counter2 + 1
			#print "counter 2: " + str(counter2)
			l2 = file2_index[line2_index]
			line2 = l2.split("\t")
			# Make sure there is a chance of overlap, begin of line 1 needs to be smaller then line2
			if (line1[0] == line2[0]):
				found_overlap = overlap(line1, line2)
				if (found_overlap):
					if not file1_values_to_remove.has_key(line_index):
						file1_values_to_remove[line_index] = "w"
					if not file2_values_to_remove.has_key(line2_index):
						file2_values_to_remove[line2_index] = "w"	
					#print "found!"
					#print l1
					#print l2


	#add the rest of the lines that didn't overlap
	count1 = 0
	count2 = 0
	for line_index in file1_index.keys():
		if file1_values_to_remove.has_key(line_index):
			count1 = count1 + 1	
		else:
			line1 = file1_index[line_index].strip().split("\t")
			outfile.write("\t".join(line1) + "\t-\t-\t-\n")
			#outfile.write("\t".join(line1[:8]) + "\t-\t-\t-\n")

	for line2_index in file2_index.keys():
		if file2_values_to_remove.has_key(line2_index):
			count2 = count2 + 1
			#print "overlaped"
		else:
			line2 = file2_index[line2_index].strip().split("\t")
			#outfile.write("\t".join(line2[:5]) + "\t-\t-\t-\t" + "\t".join(line2[5:8]) + "\n")
			#outfile.write("\t".join(line2[:5]) + "\t-\t-\t-\t-\t-\t-\t" + "\t".join(line2[5:8]) + "\n")
			outfile.write("\t".join(line2[:5]) + "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t" + "\t".join(line2[5:8]) + "\n")
			#outfile.write("\t".join(line2[:5]) + "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t" + "\t".join(line2[5:]) + "\n")	

	print "Number of lines from file 1 overlaped = " + str(count1)
	print "Number of lines from file 2 overlaped = " + str(count2) 		
	inFile.close()
	fileTwo.close()
	outfile.close()
