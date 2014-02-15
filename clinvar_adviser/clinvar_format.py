'''
Script that formats clinvar inFileauthor: gerikson
date: 14 February 2014
'''


import sys

inFilename = str(sys.argv[1])
outfilename = str(sys.argv[2])

inFile = open(inFilename) 
outFile = open(outfilename, 'w')
inLine = inFile.readline()
counter = 0
while inLine != "":
	if counter == 0:
		counter = counter + 1
		inLine = inFile.readline()
		continue
	try:
		counter = counter +1
		if (counter%10000):
			print str(counter)	
		inLine = inLine.strip().split("\t")
		tempL = inLine[0].split("_")
		value = inLine[1].split(")~")
		v1 = value[0].split("(")
		gene = v1[1]
		outFile.write(tempL[0] + "\t" + tempL[1] + "\t" + tempL[2] + "\t" + tempL[3] + "\t" + tempL[4] + "\t" + tempL[5] + "\t" + inLine[1] + "\t" + gene + "\n")
		outFile.flush()
		inLine = inFile.readline()
	except:
		print inLine
inFile.close()
outFile.close()
