'''
Script that converts a file generated by CNVnator to the pipeline input file 
Author: Galina Erikson
September 10th 2013
'''


#!/usr/bin/python
import sys, json

def transform_file(inLine):
    '''
    put the tranformation in this function eventually
    '''
    
inFilename = str(sys.argv[1])
outFilename = str(sys.argv[2])

inFile = open(inFilename)
outFile = open(outFilename, 'w')

header = "Haplotype       Chromosome      Begin   End     VarType Reference       Allele  Notes"
outFile.write(header + "\n")
inLine = inFile.readline()

while inLine.strip() != '':
    inLin = inLine.strip().split('\t')
    
    # make sure there are enough columns in file
    if len(inLin) < 8:
        print "Error! Not enough values in line"
        inLine = inFile.readline()
        continue
        
    #extract the CNV type    
    if (inLin[0] == "deletion"):
        CNV_type = "loss"
    elif (inLin[0] == "duplication"):
        CNV_type = "gain"
    else:
        print "Error! Neither deletion or duplication " + str(inLine[0])
        inLine = inFile.readline()
        continue
        
        
    # extract the chromosome
    coordinates = inLin[1].split(":")
    chr = "chr" + str(coordinates[0])
    coord = coordinates[1].split("-")
    begin = coord[0]
    end = coord[1]
    
    #adding q0, pvalues, lenght, normalized_rd to the Notes column in json format
    js = {}
    js['q0'] = inLin[8]
    js['p-val1'] = inLin[4]
    js['p-val2'] = inLin[5]
    js['p-val3'] = inLin[6]
    js['p-val4'] = inLin[7]
    js['length'] = inLin[2]
    js['normalized_rd'] = inLin[3]
    
    outFile.write("-" + "\t" + chr + "\t" + begin + "\t" + end + "\t" + CNV_type + "\t" + "-" + "\t" + "-" + "\t"  + json.dumps(js) + "\n")
    inLine = inFile.readline()
    continue

inFile.close()
outFile.close()