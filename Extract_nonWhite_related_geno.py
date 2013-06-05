#!/usr/bin/python
import  os, sys 

'''
Script that takes two files: 1. reduced genotype
2.all genotypes
and extract the genotypes from file2 that were not found in file1
python compare_allele_to_existing_tabix.py /path/to/reduced_file /path/to/all_genotypes /path/to/output/file
'''

reduced_file = str(sys.argv[1])
allgeno_file = str(sys.argv[2])
resultfile = str(sys.argv[3])

infile = open(reduced_file, 'r')
infiletwo = open(allgeno_file, 'r')
outfile = open(resultfile, 'w')

'''
read file 1, store into an array
'''

lin = infile.readline()
line = lin.strip().split("\t")

allgeno = infiletwo.readline()
allgen = allgeno.strip().split("\t")

'''
for each element in the second file
'''
counter = 0
for j in range(0, len(allgen)):
    '''
    compare it to each element of the first file
    '''
    '''
    counter that will determine if a genotype was found or not
    '''
    found = 0
    for i in range(0, len(line)):
        
        if line[i]==allgen[j]:
            found = 1
            
    '''
    check if the genotype was found if the reduced dataset, if not found, copy it to the result file
    '''
     
    if found == 0:
        counter = counter +1
        outfile.write(allgen[j] + "\t")

print "Number of genotypes excluded is: " + str(counter)

infile.close()
infiletwo.close()
outfile.close()                
