#!/usr/bin/python
import  os, sys 

'''
Script that reads in a file and removes all of the lines that just contain N/A 
'''

dbSNP_file = str(sys.argv[1])
parsed_file = str(sys.argv[2])


infile = open(dbSNP_file, 'r')
outfile = open(parsed_file, 'w')


lin = infile.readline()
#counter_allele = 0
while lin:
    #counter_allele = counter_allele + 1
    line = lin.strip()
    if line != "N/A":   
        outfile.write(line + "\n")
    lin = infile.readline()
    continue

infile.close()
outfile.close()
   