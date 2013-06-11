#!/usr/bin/python
import  os, sys, time 

'''
Script that takes a variant file a verifies it with the hg19 fasta chromosome to see if if the variant file is correct

'''

allele_file = str(sys.argv[1])
chrom_fast = str(sys.argv[2])
bad_variants = str(sys.argv[3])
chrom = str(sys.argv[4])

infile = open(allele_file, 'r')
chromf = open(chrom_fast, 'r')
outfile = open(bad_variants, 'w')

lin = infile.readline()
#chromFile = chromf.readLine()
'''
reading the chromosome file
'''
seq = chromf.read().replace(">","").replace(chrom,"").replace("\n","").upper()
good_var = 0
bad_var=0
counter = 0

while lin:
    counter = counter + 1
    if counter == 1:
	lin = infile.readline()
	continue	
    if counter % 1000 == 0:
        print str(counter) + time.asctime( time.localtime(time.time()) )
    line = lin.strip().split("\t")
    begin = int(line[1])
    ref = line[4]
    if seq[begin] == ref[0]:
        good_var = good_var  + 1
    else:
        if ref[0] == "-":
            good_var = good_var + 1
        else:
            bad_var = bad_var + 1
            outfile.write(lin + "\n")
    lin = infile.readline()        

print "Good variants are " + str(good_var)
print "Bad variants are " + str(bad_var)

infile.close()
outfile.close()
chromf.close()
