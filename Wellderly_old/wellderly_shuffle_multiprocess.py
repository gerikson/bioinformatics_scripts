'''
Script that shuffles the alleles in a vcf file using multiprocessing
author: gerikson
date: 14 February 2014
'''

import sys, gzip, random, datetime, traceback, os, math, pysam, pickle, re
import multiprocessing as mp
import subprocess

def shuffle(line_sampleIDdict):
	line, w_indeces = line_sampleIDdict
	line = line.split()
	'''
	if len(line) < 10:
		print "line size is: " + str(len(line))
		line = ""
	'''
	#else:
	line[7] = "."
	line[8] = line[8] + ":CG"
	genotype = line[9:]
	temp_gen = 0
	white_gen = []
	for i, gen in enumerate(genotype):
		#check if the white indeces dictionary if we have this index, if so, store the genotype
		if i in w_indeces.keys():
			white_gen.append(gen)
			temp_gen = temp_gen + 1
	#print str(temp_gen)
	white_gen = find_half(white_gen)
	#print len(white_gen)
	random.shuffle(white_gen)
	line = "\t".join(line[:9]) + "\t" + "\t".join(white_gen)
	'''
	if temp_gen != count_gen:
		print "Error, not all genotypes were found!!!!!"
	'''
	return line

def find_half(genotype):
	for index, gen in enumerate(genotype):
		gen_split = gen.split(":")
		haplotype = gen_split[0]
		# If genotype contains only one dot, replace it with ./. and add the genotype at the end
		if haplotype.count(".") == 1:
			gen_split[0] = "./."
		final_genotype = ":".join(gen_split) + ":" + haplotype
		genotype[index] = final_genotype
	return genotype

def extract_white_index():
	for line in white:
		line = line.strip()
		line = line + "-ASM"
		#put the genotypes into a dictionary witn the whites being the index
		white_index[line] = "w"

inFilename = str(sys.argv[1]) 
outfilename = str(sys.argv[2])

infile = gzip.open(inFilename, 'rb')
outfile = gzip.open(outfilename, 'wb')
white = open("/gpfs/home/gerikson/CNV_pipeline/Wellderly_scripts/wellderly_ids_nonRelated_european_453_2014_ref_snps.txt", 'r')

white_index = {}
extract_white_index()
w_indeces = {}

print 'SHUFFLE START'
print datetime.datetime.now().time()


def main():

	count_gen = 0
	count = mp.cpu_count()
	print str(count)
	pool = mp.Pool(processes=(count - 3))
	chunk = []
	for inLine in infile:

		#if this is the line with the genotypes name
		if inLine.startswith("#CHROM"):
			inLine = inLine.rstrip("\n")
			#print "Chrom line!"
			inLine = inLine.split("\t")
			string_gen = []
			#extract only the white genotypes
			genotype = inLine[9:]
			for i, gen in enumerate(genotype):
				#check if the white dictionary if we have this key, store the index into a separate dictionary
				if white_index.has_key(gen):
					string_gen.append(gen)
					w_indeces[i] = "sf"
					count_gen = count_gen + 1

			outfile.write("\t".join(inLine[:9]) + "\t" + "\t".join(string_gen) + "\n")
			outfile.flush()
			continue

		if inLine.startswith("##"):

			outfile.write(inLine)
			outfile.flush()
			continue

		inLine = inLine.strip()
		chunk.append((inLine, w_indeces))

		if not len(chunk)%10000:
			#print "chunk 1000"
			results = []
			r = pool.map_async(shuffle, chunk, callback=results.extend)
			#results.extend(r)
			r.wait()
			chunk = []
			
			for res in results:
				outfile.write(res + "\n")
				outfile.flush()
			results = []	
			
			#print results
			#results.extend(r)
	results = []		
	r = pool.map_async(shuffle, chunk, callback=results.extend)
	r.wait()
	#results.extend(r)
	#r.wait()
	chunk = []
	
	for res in results:
		outfile.write(res + "\n")
		outfile.flush()
	
	#print results

	pool.close()
	pool.join()
	#print(results)
	print "SHUFFLE End"
	print datetime.datetime.now().time()
	print "Total number of genotypes found is: " + str(count_gen)
	infile.close()
	outfile.close()

	'''
    largefile = 'test.dat'
    num_chunks = 10
    results = []
    with open(largefile) as f:
        reader = csv.reader(f)
        chunks = itertools.groupby(reader, keyfunc)
        while True:
            # make a list of num_chunks chunks
            groups = [list(chunk) for key, chunk in
                      itertools.islice(chunks, num_chunks)]
            if groups:
                result = pool.map(worker, groups)
                results.extend(result)
            else:
                break
    pool.close()
    pool.join()
    print(results)
    '''


if __name__ == '__main__': 
	main()

