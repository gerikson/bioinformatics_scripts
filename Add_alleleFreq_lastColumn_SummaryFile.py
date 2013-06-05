#!/usr/bin/python
import  os, sys 

'''
Script that takes the NEW wellderly allele frequency, compare and adds it to a preexisting SG Adviser 
SUMMARY file as the last column
'''

allele_file = str(sys.argv[1])
tabix_file = str(sys.argv[2])
resultfile = str(sys.argv[3])

infile = open(allele_file, 'r')
tabixfile = open(tabix_file, 'r')
outfile = open(resultfile, 'w')

'''
Create a dictionary with the allele frequency file
'''
lin = infile.readline()
allele_dict = dict()
counter_allele = 0

while lin:
    counter_allele = counter_allele + 1
    line = lin.strip().split("\t")
    templ = line[0]+"_"+line[1]+"_"+line[2]+"_"+line[3]+"_"+line[4]+"_"+line[5]
    allele_dict[templ] = line[6]
    '''
    test, is anything happening
    '''
    if counter_allele % 10000 == 0:
        print "Creating a dictionary of  allele frequency" + str(counter_allele)
    lin = infile.readline()

 
'''
Read in the old tabix file and check each line if it's present in the allele frequency, if true, add new frequency
else add N/A
'''

lin = tabixfile.readline()
counter_tabix = 0

while lin:
    lin = lin.strip()	
    counter_tabix = counter_tabix +1
    if counter_tabix % 10000 == 0:
        print "Checking tabix file " + str(counter_tabix)	
    line = lin.split("\t")
    '''
    If this is the first line, skip
    '''
    if line[0] == 'Haplotype':
        outfile.write(lin + '\t' + "Wellderly_AF325" + "\n")
        lin = tabixfile.readline()
        continue
    else :
        templ = line[1]+"_"+line[2]+"_"+line[3]+"_"+line[4]+"_"+line[5] + "_" + line[6]
        '''
        check if that value exist in the allele dictionary
        '''
        if templ in allele_dict:
             outfile.write(lin + '\t' + allele_dict[templ] + "\n")
             '''
             We might have to remove this value from dictionary
             '''
             # del allele_dict[key]
             lin = tabixfile.readline()
             continue
        else:
              outfile.write(lin + '\t' + "N/A" + "\n")
              lin = tabixfile.readline()
              continue

'''
need to iterate through the dictionary and print the lines that 
'''

infile.close()
outfile.close()
tabixfile.close()
    
