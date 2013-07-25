#!/usr/bin/python
import  os, sys
from Bio.Seq import Seq 

'''
Script that takes the dbSNP file from /gpfs/home/gerikson/dbSNP/snp137.txt and parses it for SG Adviser annotation
'''

dbSNP_file = str(sys.argv[1])
parsed_file = str(sys.argv[2])
weird_variants = str(sys.argv[3])

infile = open(dbSNP_file, 'r')
outfile = open(parsed_file, 'w')
out_bad = open(weird_variants, 'w')

lin = infile.readline()
counter_allele = 0
bad_counter = 0

#begin = ""
#endpos = ""
#start = 0
#end = 0
#re = ""
#var = ""    

def transform(line):
        found = 0
	begin = 0

        '''
        Need to do some transformations     
        '''
                        #if len(var) == len(ref) and len(var) > 1:
                    # handles messy cases like:
                    # ref =   'TGCCGTGGCGGCA'
                    # vars = ['TGCCGTGGCGGAA','TGCCGTGGTGGCA'] 
                    #  TGCCGTGGCGGCA
                    #             C
                    #  TGCCGTGGCGGAA
                    #             A
                    #  ************
                    #  TGCCGTGGCGGCA
                    #          C
                    #  TGCCGTGGTGGCA
                    #          T
                    #  ************        
        try:	
	    begin = int(line[2])
            endpos = int(line[3])
	except:
	    #bad_counter = bad_counter + 1
            #print "Counter = " + str(counter_allele)
            #print "Bad counter = " + str(bad_counter)
            #out_bad.write(lin + "\n")
            #break
 	    endpos = begin + len(line[8])	
	start = 0
        end = 0
        re = line[8]
        ref = line[9].split("/")
         
        for var in ref:    
            if var != "-":                 
                for i in range(0,min( len(var),len(re)) ):
                            if var[i] == re[i]:
                                start += 1
                            else:
                                break    
                re = re[start:]
                var = var[start:]
                begin += start
        
                    # check the tail    
                for i in range(1,min( len(var),len(re)) ):
                    if var[-i] == re[-i]:
                        end -= 1
                    else:
                        break
                           
                if end == 0:
                    # if end is > -1 then we have no matching tail
                    pass
                else:
                    var = var[:end]
                    re = re[:end]
                    
                endpos +=end        
                
                if len(re) == 0:           
                     found = 1
                     outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + str(begin) + "\t" + str(endpos) + "\t" + "ins" + "\t"  + "-" + "\t" + var + "\t" + line[4] + "\n")
                elif len(var) == 0:
                     found =1
                     outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + str(begin) + "\t" + str(endpos) + "\t" + "del" + "\t"  + re + "\t" + "-" + "\t" + line[4] + "\n")
                elif len(var) == 1 and len(re) == 1 and var != re:
                     found = 1
                     outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + str(begin) + "\t" + str(endpos) + "\t" + "snp" + "\t"  + re + "\t" + var + "\t" + line[4] + "\n")
                elif var != re:
                     found = 1
                     outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + str(begin) + "\t" + str(endpos) + "\t" + "delins" + "\t"  + re + "\t" + var + "\t" + line[4] + "\n")

        if found == 0:
            bad_counter = bad_counter + 1
            print "Counter = " + str(counter_allele)
            print "Bad counter = " + str(bad_counter)
            out_bad.write(lin + "\n")                       

    
while lin:
    counter_allele = counter_allele + 1
    line = lin.strip().split("\t")
    if line[11] == "single":
        found = 0
        #this might be something other then a snp
        if len(line[8]) > 1:
            transform(lin)
            lin = infile.readline()
            continue 
            
        ref = line[9].split("/")
	for var in ref:
	    '''
            Check to see if this is on the negative strand
            '''
            if line[6] == "-":
		#Transform var to positive strand
		var = str(Seq(var).reverse_complement())		
            if var != line[8]:
               found = 1
	       outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + "snp" + "\t"  + line[8] + "\t" + var + "\t" + line[4] + "\n") 		
         
        if found == 0:
            bad_counter = bad_counter + 1
            print "Counter = " + str(counter_allele)
            print "Bad counter = " + str(bad_counter)
            out_bad.write(lin + "\n")              
               		
        lin = infile.readline()
        continue   

    elif line[11] == "deletion":
            found = 0
            # if there is no "-" in 10th column this might be something else other then a deletion
            if not "-" in line[9]:
                transform(lin)
                lin = infile.readline()
                continue 
            
	    ref = line[9].split("/")	
            for var in ref:
             
                '''
                Check to see if this is on the negative strand
                '''
                if line[6] == "-" and var != "-":
                   #Transform var to positive strand
                   var = str(Seq(var).reverse_complement())        
                if var != "-":
                   found = 1
                   outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + "del" + "\t"  + "-" + "\t" + var + "\t" + line[4] + "\n")
                   
            if found == 0:
               bad_counter = bad_counter + 1
               print "Counter = " + str(counter_allele)
               print "Bad counter = " + str(bad_counter)
               out_bad.write(lin + "\n")        
               
            lin = infile.readline()
            continue
               
    elif line[11] == "insertion":
        found = 0
        # if there is no "-" in 10th column this might be something else other then a deletion
        if not "-" in line[9]:
            transform(lin)
            lin = infile.readline()
            continue
 
        ref = line[9].split("/")
        for var in ref:    
            '''
            Check to see if this is on the negative strand
            '''
            if line[6] == "-" and var != "-":  
                #Transform var to positive strand
                var = str(Seq(var).reverse_complement())        
            if var != "-":
                found = 1
                outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + "ins" + "\t"  + "-" + "\t" + var + "\t" + line[4] + "\n")
	if found == 0:
            bad_counter = bad_counter + 1
            print "Counter = " + str(counter_allele)
            print "Bad counter = " + str(bad_counter)
            out_bad.write(lin + "\n")
            
        lin = infile.readline()
	continue
    
    elif line[11] == "in-del" or line[11] == "mnp" or line[11] == "mixed":
	    transform(line) 
            lin = infile.readline()
            continue    
	
    else:
	bad_counter = bad_counter + 1
        print "Counter = " + str(counter_allele)
        print "Bad counter = " + str(bad_counter)
        out_bad.write(lin + "\n")
        lin = infile.readline()    
	continue

infile.close()
outfile.close()
out_bad.close()   
