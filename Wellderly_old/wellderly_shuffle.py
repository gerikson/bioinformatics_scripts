import sys, gzip, random, datetime
import multiprocessing
import subprocess

def shuffle(line):
	#print line
	#line.replace("|", "/")
	#line = line.split("\t")
	line = line.split()
	if len(line) < 543:
		print "line size is: " + str(len(line))
		print line
		#print line
	else:
		line[8] = line[8] + ":CG"
		genotype = line[9:]
		genotype = find_half(genotype)
		random.shuffle(genotype)
		line = "\t".join(line[:9]) + "\t" + "\t".join(genotype)
		outfile.write(line + "\n")
		outfile.flush()

def find_half(genotype):
	for index, gen in enumerate(genotype):
		gen_split = gen.split(":")
		haplotype = gen_split[0]
		# If genotype contains only one dot, replace it with ./. and add the genotype at the end
		if haplotype.count(".") == 1:
			#print "found" + haplotype
			gen_split[0] = "./."
			#print gen_split[0]
		final_genotype = ":".join(gen_split) + ":" + haplotype
		genotype[index] = final_genotype
	return genotype

inFilename = str(sys.argv[1]) 
outfilename = str(sys.argv[2])

infile = gzip.open(inFilename, 'rb')
outfile = gzip.open(outfilename, 'wb')

print 'SHUFFLE START'
print datetime.datetime.now().time()

#inLine = infile.readline()
#counter = 0
#count_al = 0

#while inLine != "":
for inLine in infile:
	'''
	count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=(count - 3))
    results = []
    r = pool.map_async(calculate, tasks, callback=results.append)
    r.wait() # Wait on the results
    print results
	'''
	if inLine.startswith("#"):
		outfile.write(inLine + "\n")
		inLine = infile.readline()
		continue
	inLine = inLine.strip()
	shuffle(inLine)
	'''
	try:
		shuffle(inLine)
	except:
		print "Error"
		#print inLine
	'''
	#pool = 
	#inLine = infile.readline()

print 'SHUFFLE End'
print datetime.datetime.now().time()

infile.close()
outfile.close()