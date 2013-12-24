# Script tha splits a file by chromosome
# Galina Erikson
# 23 December 2013

import sys

def splitFile(line, chromColumn, chrom, outFile):
	ch = ""
	if line[chromColumn].startswith("chr"):
		ch = line[chromColumn][3:]
	else:
		ch = line[chromColumn]
	
	if ch != chrom:
		chrom = ch
		outFile.close()
		outfilename = outpath + "/chr" + chrom + "_" + inFilename
		outFile = open(outfilename, 'w')
	tempL = "\t".join(line)
	outFile.write(tempL + "\n")
	outFile.flush()	
		
'''
Takes 3 arguments:
1. Input inFilename
2. output folder path
3. column number where the chromosome is located
'''

inFilename = str(sys.argv[1])
outpath = str(sys.argv[2])
chromColumn = int(sys.argv[3])

inFile = open(inFilename) 
inLine = inFile.readline()
if inLine.startswith("variantaccession"):
	inLine = inFile.readline()

chrom = "1"
outfilename = outpath + "/chr" + chrom + "_" + inFilename
outFile = open(outfilename, 'w')

while inLine != "":
	Line = inLine.split('\t')
	ch = ""
	if Line[chromColumn].startswith("chr"):
		ch = Line[chromColumn][3:]
	else:
		ch = Line[chromColumn]
	
	if ch != chrom:

		chrom = ch
		outFile.close()
		outfilename = outpath + "/chr" + ch + "_" + inFilename
		outFile = open(outfilename, 'w')

	tempL = "\t".join(Line)
	outFile.write(tempL)
	#splitFile(Line, chromColumn, chrom, outFile)
	inLine = inFile.readline()

inFile.close()
outFile.close()




