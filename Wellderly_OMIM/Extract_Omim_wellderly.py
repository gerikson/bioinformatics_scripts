import os, sys,math, time, pymongo,json,bson

HOME = "/gpfs/group/stsi/methods/annotation/sg-adviser/dev/"
SCHORK_HOME = "/gpfs/home/nschork/"

sys.path.append(SCHORK_HOME + "biopython-1.54/")
from collections import defaultdict

print "START: "
print time.asctime( time.localtime(time.time()) )
varfilename = str(sys.argv[1])
resfilename =str(sys.argv[2])


print "Extracting OMIM variants: ", varfilename

outres = open(resfilename, 'w')

invar = open(varfilename)
line = invar.readline()
counter = 0

while line != "":
    counter = counter + 1
    if counter % 100000 == 0:
        print str(counter)
        if line.startswith("chrom"):
            line = invar.readline()
            continue
    lin = line.strip().split("\t")
    omim = lin[60]
   	
    #if omim.isdigit() or omim.isupper or omim.islower:
    for c in omim:
	 if c.isalnum() or c.isalpha():
         	#print omim
	 	outres.write(line)
	 	#outres.write("\n")
         	#line = invar.readline()
         	#continue
     		break

    line = invar.readline()
    continue

outres.close()

print "END: ", time.time()
print time.asctime( time.localtime(time.time()) )

sys.stdout.flush()
