'''
Script that fixes tabix file by replacing missing alleles with the most frequent allele
author: Galina Erikson
'''


import sys

def frequent_char(text, frequent_allele):
	#print "text is: " + text
	charset = set(sorted(text))
	textTemp = text.replace(" ", "")  # add this to remove spaces
	#print textTemp
	#textTemp = ''.join(text)
	#temp_maxcount = 0
	maxcount = 0
	maxchar = ""
	for item in charset:
		charcount = textTemp.count(item)
		#print "Item is: " + item
		if frequent_allele == 'N/A' or item != frequent_allele:
			if charcount > maxcount:
				maxcount = charcount
				maxchar = item


	return maxchar

inFilename = str(sys.argv[1])
outfilename = str(sys.argv[2])

inFile = open(inFilename) 
outFile = open(outfilename, 'w')
inLine = inFile.readline()
counter = 0
count_al = 0
while inLine != "":

	counter = counter + 1
	tempLine = inLine.split("\t")
	text = tempLine[4:]
	tempText = "".join(text)
	tempText = tempText.replace(" ", "")
	#Find the two most comon alleles
	freq_allele1 = frequent_char(tempText, 'N/A')
	#print "Most frequent " + freq_allele1
	freq_allele2 = frequent_char(tempText, freq_allele1)
	#print "Second most frequent " + freq_allele2
	if counter%1000 == 0:
		print str(counter)
	for index, item in enumerate(text):
		item_list = item.split(' ')
		for allele in item_list:
			#print "allele is: " + allele
			allele = allele.strip()
			if allele != freq_allele1 and allele != freq_allele2:
				count_al = count_al + 1
				#print "line number " + str(counter) + " " + item
				text[index] = "0 0"
	outFile.write("\t".join(tempLine[:4]) + "\t" + "\t".join(text))
	outFile.flush()
	inLine = inFile.readline()
	
print "Missing alleles " + str(count_al)
inFile.close()
outFile.close()

	




