#!/usr/bin/python
import  os, sys 

'''
Script that takes the correctly parsed rsID file, saves the rsID in a dictionary and then loads in the file of unparsed rsID and checks line by line if the rsID was previously found, if not, it saves it to a file.
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
    templ = line[7]
    allele_dict[templ] = "s"
    '''
    test, is anything happening
    '''
    if counter_allele % 10000 == 0:
        print "Creating a dictionary of  allele frequency" + str(counter_allele)
    lin = infile.readline()

 
'''
Read in file with the rsID that were not parsed correctly, and see if they are found in the allele dictionary
'''

lin = tabixfile.readline()
counter_tabix = 0
found = 0
not_found = 0

while lin:
    if found % 10000 == 0 or   not_found % 10000 ==0
         print "Found " + str(found)
         print "Not found " + str(not_found)
         
    lin = lin.strip()	
    line = lin.strip().split("\t")
    templ = line[7]
        if templ in allele_dict:
             found++
             lin = tabixfile.readline()
             continue
        else:
              not_found++  
              outfile.write(lin + '\t' + "N/A" + "\n")
              lin = tabixfile.readline()
              continue
             
infile.close()
outfile.close()
tabixfile.close()
    
