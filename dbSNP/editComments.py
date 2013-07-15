#!/usr/bin/python
import  os, sys 

'''
Script that takes /gpfs/home/gerikson/dbSNP/snp137_annotation/annotClean and edits the comments (everything past the question mark)
'''

dbSNP_file = str(sys.argv[1])
parsed_file = str(sys.argv[2])


infile = open(dbSNP_file, 'r')
outfile = open(parsed_file, 'w')


lin = infile.readline()
#counter_allele = 0
while lin:
    #counter_allele = counter_allele + 1
    line = lin.strip().split("\t")
    com = line[7].strip().split("?")
    outfile.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t"  + line[5] + "\t" + line[6] + "\t" + com[0] + "\n")
    lin = infile.readline()
    continue

infile.close()
outfile.close()
   
