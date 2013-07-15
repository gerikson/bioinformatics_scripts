#!/usr/bin/python
import  os, sys 

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
while lin:
    counter_allele = counter_allele + 1
    line = lin.strip().split("\t")
    if line[11] == "single":
        ref = line[9].split("/")
	if len(ref)>2:
	    for var in ref:
		if var != line[7]:
		    outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + "snp" + "\t"  + line[7] + "\t" + var + "\t" + line[4] + "\n")
		    lin = infile.readline()
		    continue 				
	    '''
	    bad_counter = bad_counter + 1
	    print "Counter = " + str(counter_allele)
	    print "Bad counter = " + str(bad_counter)		 	
	    out_bad.write(lin + "\n")
	    lin = infile.readline()
            continue
	    '''
	else:
            outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + "snp" + "\t"  + ref[0] + "\t" + ref[1] + "\t" + line[4] + "\n")
	    lin = infile.readline()
	    continue
    elif line[11] == "deletion":
	ref = line[9].split("/")	
	if len(ref)>2:
            bad_counter = bad_counter + 1
            print "Counter = " + str(counter_allele)
            print "Bad counter = " + str(bad_counter)
            out_bad.write(lin + "\n")
            lin = infile.readline()
            continue
	else:
            outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + "del" + "\t"  + ref[1] + "\t" + ref[0] + "\t" + line[4] + "\n")
	    lin = infile.readline()
	    continue
    elif line[11] == "insertion":
        ref = line[9].split("/")
	if len(ref)>2:
            bad_counter = bad_counter + 1
            print "Counter = " + str(counter_allele)
            print "Bad counter = " + str(bad_counter)
            out_bad.write(lin + "\n")
            lin = infile.readline()
            continue
        else:
            outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + "ins" + "\t"  + ref[0] + "\t" + ref[1] + "\t" + line[4] + "\n")
	    lin = infile.readline()
	    continue
    elif line[11] == "in-del":
        ref = line[9].split("/")
        if len(ref)>2:
            bad_counter = bad_counter + 1
            print "Counter = " + str(counter_allele)
            print "Bad counter = " + str(bad_counter)
            out_bad.write(lin + "\n")
            lin = infile.readline()
            continue
	elif len(ref)<2:
            bad_counter = bad_counter + 1
            print "Counter = " + str(counter_allele)
            print "Bad counter = " + str(bad_counter)
            out_bad.write(lin + "\n")
	    print lin	
            lin = infile.readline()
            continue
        else:
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
	    begin = int(line[2])
	    endpos = int(line[3])
	    start = 0
	    end = 0
	    re = line[8]
	    var = ref[1]	
	     		 				
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
            #start = str(start)
            #end = str(end)
	    
	    if len(re) == 0:		 	
            	outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + str(begin) + "\t" + str(endpos) + "\t" + "ins" + "\t"  + "-" + "\t" + var + "\t" + line[4] + "\n")
                lin = infile.readline()
                continue
            elif len(var) == 0:
                outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + str(begin) + "\t" + str(endpos) + "\t" + "del" + "\t"  + re + "\t" + "-" + "\t" + line[4] + "\n")
                lin = infile.readline()
                continue
	    elif len(var) == 1 and len(re) == 1:
                outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + str(begin) + "\t" + str(endpos) + "\t" + "snp" + "\t"  + re + "\t" + var + "\t" + line[4] + "\n")
                lin = infile.readline()
                continue	
            else:
                outfile.write(str(counter_allele) + "\t" + line[1] + "\t" + str(begin) + "\t" + str(endpos) + "\t" + "delins" + "\t"  + re + "\t" + var + "\t" + line[4] + "\n")
                lin = infile.readline()
                continue    	
	
    #elif line[10] == "genomic":
#	print lin + "\n"
#	lin = infile.readline()
#	continue
        #ref = line[9].split("/")
        #outfile.write(counter_allele + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + "delins" + "\t"  + ref[1] + "\t" + ref[0] + "\n")
    else:
	bad_counter = bad_counter + 1
        print "Counter = " + str(counter_allele)
        print "Bad counter = " + str(bad_counter)
        out_bad.write(lin + "\n")
        lin = infile.readline()    
	continue

infile.close()
outfile.close()
   
