'''
Class that compares the final combined 1000G CNV file for overlaps within the file

author: gerikson
Date: March 10, 2014

'''

import sys

def overlap(line1, line2):
	start1 = int(line1[1])
	start2 = int(line2[1])
	end1 = int(line1[2])
	end2 = int(line2[2])
	len1 = end1 - start1
	len2 = end2 - start2
	'''
	cnv1 -------
	cnv2 -------
	'''
	if start1 == start2 and end1 == end2:
		#outfile.write("\t".join(line1) + "\t"+ line2[5] + "\t" + line2[6] + "\t" + line2[7] + "\n")
		return True
		'''
		cnv1  ---------------------
		cnv2  		 ----------------------
		'''	
	elif start1 <= start2 and end1 <= end2 and start2 < end1:
		size_overlap = end1 - start2
		if float(size_overlap)/float(len1) >= 0.8 and float(size_overlap)/float(len2) >= 0.8:
			return True
		'''
		cnv1           ---------------------
		cnv2   ----------------------
		'''
	elif start2 <= start1 and end2 <= end1 and start1 < end2:
		#print "Category 2"
		size_overlap = end2 - start1
		if float(size_overlap)/float(len1) >= 0.8 and float(size_overlap)/float(len2) >= 0.8:
			return True
		'''
		cnv1        ------------------
		cnv2    ----------------------------
		'''    
	elif start1 >= start2 and end1 <= end2:
		if float(len1)/float(len2) >= 0.8:
			return True

		'''
		cnv1 --------------------------------
		cnv2          -----------------
		'''    
	elif start1 <= start2 and end1 >= end2:
		if float(len2)/float(len1) >= 0.8:
			return True
	else:
		return False

def indexFile(file_name):
	dict_name = {}

	for line in file_name:
		line = line.strip()
		if not line.startswith("Chr"):
			#temp_line = line.split("\t")
			#dict_key = temp_line[0]+ "_" + temp_line[1]+"_"+temp_line[2]
			#dict_name[dict_key] = line
			dict_name[line] = line
	return dict_name

def extract_max_overlap(dict_of_overlaps, non_overlap_counter):
	temp_counter = int(non_overlap_counter)
	count = 0
	begin_array = []
	end_array = []
	new_line = []
	chrom = ""
	vartype = ""
	for line in dict_of_overlaps.keys():
		l1 = line.split("\t")
		chrom = l1[0]
		vartype = l1[3]
		begin_array.append(l1[1])
		end_array.append(l1[2])
	begin = min(begin_array)
	#print str(begin)
	end = max(end_array)
	#print str(end)

	new_line.append(chrom)
	new_line.append(str(begin))
	new_line.append(str(end))

	#population = "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t"
	population = "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t"
	#population = "-\t-\t-\t-\t-\t-\t-\t-\t-"
	pop = population.split("\t")
	for line in dict_of_overlaps.keys():
		l1 = line.split("\t")
		for index, var in enumerate(l1[5:]):
			if var != "-":
				if pop[index] == "-":
					pop[index] = var
				else:
					temp = pop[index]
					pop[index] = temp + '///' + var

	#overlaps.write("\t".join(new_line) + "\t" + vartype + "\tPASS\t" + "\t".join(pop) + "\n")
	outfile.write("\t".join(new_line) + "\t" + vartype + "\tPASS\t" + "\t".join(pop) + "\n")
	return temp_counter

def min(arr):
	minim = sys.maxint
	for num in arr:
		if int(num) < minim:
			minim = num
	#print str(minim)
	return minim

def max(arr):
	maxim = 0
	for num in arr:
		if int(num) > maxim:
			maxim = num
	#print str(maxim)
	return maxim



inFilename = str(sys.argv[1])
outfilename = str(sys.argv[2])
#ov = str(sys.argv[3])

inFile = open(inFilename)
outfile = open(outfilename, 'w')
#overlaps = open(ov, 'w')

if __name__ == "__main__":

	non_overlap_counter = 0
	file1_index = {}

	#each time we find an overlap we have to remove those keys so we don't add them to the file twice
	file1_values_to_remove = {}

	file1_index = indexFile(inFile)

	#outfile.write("Chromosome\tBegin\tEnd\tVarType\tFilter\tASW_ID\tASW_Sample_Count\tASW_Total_Sample\tLWK_ID\tLWK_Sample_Count\tLWK_Total_Sample\n")
	#outfile.write("Chromosome\tBegin\tEnd\tVarType\tFilter\tASW_ID\tASW_Sample_Count\tASW_Total_Sample\tLWK_ID\tLWK_Sample_Count\tLWK_Total_Sample\tYRI_ID\tYRI_Sample_Count\tYRI_Total_Sample\n")
	counter1 = 0 
	counter2 = 0
	for line_index in file1_index.keys():
		dict_of_overlaps = {}
		overlap_array = []
		counter1 = counter1 + 1
		l1 = file1_index[line_index]
		line1 = l1.split("\t")
		for line2_index in file1_index.keys():
			if line_index != line2_index:
				l2 = file1_index[line2_index]
				line2 = l2.split("\t")
				# Make sure this is the same chromosome
				if (line1[0] == line2[0]):
					found_overlap = overlap(line1, line2)
					if (found_overlap):
						if not file1_values_to_remove.has_key(line_index):
							file1_values_to_remove[line_index] = "w"
							dict_of_overlaps[line_index] = line_index
						if not file1_values_to_remove.has_key(line2_index):
							file1_values_to_remove[line2_index] = "w"
							dict_of_overlaps[line2_index] = line2_index

		#Go through each element of the dictionary of overlaps and overlap it with each entry again, this way we will 
		#get the overlaps from the perspective of each entry					
		if len(dict_of_overlaps) > 1:
			
			temp_dict = {}
			for line in dict_of_overlaps.keys():
				#print line
				temp_dict[line] = "w"
				for line_index in file1_index.keys():
					#if it's the current line, or we already overlaped this line, skip it
					if not temp_dict.has_key(line_index):
						l1 = line.strip().split("\t")
						l2 = line_index.strip().split("\t")
						found_overlap = False
						#check to see if it's on the same chromosome
						if l1[0] == l2[0]:
							found_overlap = overlap(l1, l2)
						if (found_overlap):
							if not dict_of_overlaps.has_key(line_index):
								#print line_index
								counter2 = counter2 + 1
								print "new line: " + line_index
								if not file1_values_to_remove.has_key(line_index):
									file1_values_to_remove[line_index] = 'w'
								dict_of_overlaps[line_index] = line_index	
			temp_dict.clear()					
			
			non_overlap_counter = extract_max_overlap(dict_of_overlaps, non_overlap_counter)

		dict_of_overlaps.clear()
	
	count1 = 0
	for line_index in file1_index.keys():
		if file1_values_to_remove.has_key(line_index):
			count1 = count1 + 1	
		else:
			line1 = file1_index[line_index].strip().split("\t")
			outfile.write("\t".join(line1) + "\n")
	

	print "Number of lines from file 1 overlaped = " + str(count1)	
	print "newly overlaped found = " + str(counter2)
	print "Lines that fell off the final overlap = " + str(non_overlap_counter)
	inFile.close()
	outfile.close()
