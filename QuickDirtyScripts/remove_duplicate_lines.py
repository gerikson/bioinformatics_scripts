'''
Script that removes duplicate lines from fiel
Takes 3 arguments:
1. Input inFilename
2. output folder path

'''

import sys

inFilename = str(sys.argv[1])
outpath = str(sys.argv[2])


inFile = open(inFilename) 

outfilename = outpath + "/filtered_" + inFilename
outFile = open(outfilename, 'w')
previous_line = ""
inLine = inFile.readline()
while inLine != "":
	inLine = inLine.strip()
	if inLine != previous_line:
		outFile.write(inLine + "\n")
		previous_line = inLine
		inLine = inFile.readline()
	else:	
		inLine = inFile.readline()

inFile.close()
outFile.close()
