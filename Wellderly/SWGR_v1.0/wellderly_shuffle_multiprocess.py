'''
Script that shuffles the alleles in a vcf file using multiprocessing
and converts phased genotype calls to unphased

Ordered command-line arguments:
input file path - vcf format +/- gzip, e.g. (/path/to/in.vcf.gz)
output file path - vcf format +/- gzip, e.g. (/path/to/target_dir/out.vcf)
subset file path - newline separated file with subset sample_ids


author: gerikson
date: 14 February 2014
'''

import sys, gzip, random, datetime, traceback, os, math, pysam, pickle, re
import multiprocessing as mp
import subprocess

def shuffle(line_sampleIDdict):
	''' This function permutes the index assignments
	    for each sample's genotype/quality call '''
	line, ea_indices = line_sampleIDdict
	line = line.split()
	line[7] = "."  #removing INFO
	line[8] = line[8]  #FORMAT
	genotype = line[9:]  #Sample genotype/quality calls
	temp_gen = 0
	ea_gen = []
	for i, gen in enumerate(genotype):
		#check if sample in ea indices dictionary if we have this index, if so, store the genotype
		if i in ea_indices.keys():
			ea_gen.append(gen)
			temp_gen = temp_gen + 1
	ea_gen = find_half(ea_gen)
	random.shuffle(ea_gen)
	line = "\t".join(line[:9]) + "\t"  "\t".join(ea_gen)
	return line
	
def find_half(genotype):
	''' This function replaces half-calls with missing (".") 
	    calls '''
	for index, gen in enumerate(genotype):
		gen_split = gen.split(":")
		haplotype = gen_split[0]
		# If genotype contains only one dot, replace it with ./.
		if haplotype.count(".") == 1:
			gen_split = [".,." if "," in i else "." for i in gen_split[1:] ]
			gen_split.insert(0,"./.")
		final_genotype = ":".join(gen_split) 
		genotype[index] = final_genotype
	return genotype

def extract_ea_index():
	for line in ea:
		line = line.strip()
		line = line + "-ASM"
		#put the genotypes into a dictionary witn the eas being the index
		ea_index[line] = "w"


inFilename = str(sys.argv[1]) 
outfilename = str(sys.argv[2])
eafilename = str(sys.argv[3])

def gzip_call(filename, read_write):
	if ".gz" in filename:
		return gzip.open(filename,read_write)
	return open(filename,read_write)



infile = gzip_call(inFilename, 'rb')
outfile = gzip_call(outfilename, 'wb')
ea = gzip_call(eafilename, 'rb')

ea_index = {}
extract_ea_index()
ea_indices = {}  #european ancestry index position

print 'SHUFFLE START'
print datetime.datetime.now().time()


def main():

	count_gen = 0
	count = mp.cpu_count()  #counting total processors
	print str(count)
	pool = mp.Pool(processes=(count - 3))  #throttling processor use to total -3
	chunk = []  #iterator passed to multiprocessing pool
	for inLine in infile:  #iterate through file lines

		
		if inLine.startswith("#CHROM"):  #header line
			inLine = inLine.rstrip("\n")
			#print "Chrom line!"
			inLine = inLine.split("\t")
			string_gen = []
			#extract only the ea genotypes
			genotype = inLine[9:]
			for i, gen in enumerate(genotype):  #iterate through sample genotypes
				#check if the ea dictionary has key, store the index into a separate dictionary
				if ea_index.has_key(gen):  #check if sample in european ancestry subset
					string_gen.append(gen)  #save genotype string
					ea_indices[i] = "sf"  #saving ea index positions as dictionary keys
					count_gen = count_gen + 1

			outfile.write("\t".join(inLine[:9]) + "\t" + "\t".join(string_gen) + "\n")
			outfile.flush()
			continue

		if inLine.startswith("##"):

			outfile.write(inLine)
			outfile.flush()
			continue

		inLine = inLine.replace("|","/")  #replacing phasing information
		inLine = inLine.strip()
		chunk.append((inLine, ea_indices))

		if not len(chunk)%10000:  #multi-processing every 10000 lines
			results = []
			r = pool.map_async(shuffle, chunk, callback=results.extend)
			r.wait()
			chunk = []
			
			for res in results:  #writing results to file
				outfile.write(res + "\n")
				outfile.flush()
			results = []	
		
	results = []  #multi-processing last set of lines
	r = pool.map_async(shuffle, chunk, callback=results.extend)
	r.wait()
	chunk = []
	
	for res in results:  #writing results to file
		outfile.write(res + "\n")
		outfile.flush()

	pool.close()
	pool.join()
	print "SHUFFLE End"
	print datetime.datetime.now().time()
	#print "Total number of genotypes found is: " + str(count_gen)
	infile.close()
	outfile.close()





if __name__ == '__main__': 
	main()

