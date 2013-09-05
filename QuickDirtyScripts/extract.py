import os, sys,math, time

print "START: "
print time.asctime( time.localtime(time.time()) )

"""
Takes 3 arguments: 
"""

varfilename = str(sys.argv[1])
outfile = str(sys.argv[2])
CUR_CHROM = 'chr1'


print "Chromosome is: "+CUR_CHROM

print "CHECKING FOR PRE-EXISTING VARIANTS: ", varfilename

inserted_variants = 0
skiped_var = 0
non_ascii_var = 0
invar = open(varfilename)
outvar = open(outfile, "w")



line = invar.readline()

while line:
	"""
	If this is the first line skip
	"""	
	lin = line.strip().split("\t")
	
	if (lin[1]=="Chromosome"):		               
	        outvar.write(line)
                outvar.write("\n")
		line = invar.readline()
		continue	
	
	beg = long(lin[2])
	
	if (lin[1] == "chr2"):
			break
			
	num = 157577520	
	
	if lin[1] == "chr1" and beg > num:
		outvar.write(line)
		outvar.write("\n")
		line = invar.readline()
		continue
	else:
		line = invar.readline()
                continue
