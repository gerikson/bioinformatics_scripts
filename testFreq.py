#!/usr/bin/python
import  os, sys, time 

'''
Script that takes the new tabix file and counts how many variants previously had an allele frequency but there is no allele freq in the AF325 dataset

'''

allele_file = str(sys.argv[1])
resfilename = str(sys.argv[2])

infile = open(allele_file, 'r')
outres = open(resfilename, 'w')

lin = infile.readline()

total_var = 0
missing_AF325 = 0
missing_CG50 = 0
different_count = 0

while lin:
    total_var = total_var+1
    if total_var % 10000 == 0:
	print str(total_var) + time.asctime( time.localtime(time.time()) ) 	
    line = lin.strip().split("\t")
    if line[20] == 'N/A' and line[19] !='N/A':
        missing_AF325 = missing_AF325 + 1
	outres.write(lin)
    elif line[19] =='N/A' and line[20] != 'N/A':
        missing_CG50 = missing_CG50 + 1
    elif line[19] != line[20]:
        different_count  = different_count + 1
    lin = infile.readline()
    continue
		
infile.close()
print "Total number of variants " + str(total_var)
a = missing_AF325/total_var
print "Variants previously found in wellderly missing in the new dataset " + str(missing_AF325) + "  " + str(a)
b = missing_CG50/total_var
print "New variants found in the new dataset " + str(missing_CG50) + "  " + str(b)
c = different_count/total_var
print "Variants that have a different AF " + str(different_count) + "  " + str(c)

 
        
            
