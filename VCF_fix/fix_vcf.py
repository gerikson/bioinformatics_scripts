'''
Script that parses the VCF file and separates the complex allele in 2 separate lines, edits the genotype
(New vcf files will be easier to combine across multiple samples)
author: gerikson
date: 12 November 2014
'''

import sys, gzip, random, datetime, traceback, os, math, pysam, pickle, re
import multiprocessing as mp
import subprocess

class vcf():
	chrom = ""
	begin = ""
	ref = ""
	alt = ""
	varType = ""

	def _init_(self):
		self.chrom = ""
		self.begin = ""
		self.ref = ""
		self.alt = ""
		self.varType = ""
	
	
	def _init_(self, begin, ref, alt, var):
		self.begin = begin
		self.ref = ref
		self.alt = alt
		self.varType = var
		return self

	def getVartype(self):
		return self.varType

	def getStart(self):
		return self.begin	

	def getReference(self):
		return self.ref

	def getAllele(self):
		return self.alt

def fix_complex(line):
	#print "lah"
	alleles = line[4].split(",")

	vcfObject1 = vcf()
	vcfObject2 = vcf()
	#print "line[3] = " + line[3] 
	#print "1st allele = " + alleles[0]
	#print "2d allele = " + alleles[1]
	if len(line[3]) > 1:
		if len(alleles[0]) > 1:
			#print "trimming"
			vcfObject1 = vcfTrim(line[1], line[3], alleles[0])
		else:
			vcfObject1 = vcfObject1._init_(line[1], line[3], alleles[0], "del")

		if len(alleles[1]) > 1:
			#print "trimming"
			vcfObject2 = vcfTrim(line[1], line[3], alleles[1])
		else:
			vcfObject2 = vcfObject2._init_(line[1], line[3], alleles[1], "del")
	elif len(line[3]) == 1:	
		#print "Ref has a lenght of 1"
		if len(alleles[0]) == 1:
			vcfObject1 = vcfObject1._init_(line[1], line[3], alleles[0], "snp")
		elif len(alleles[0]) > 1:
			#print "1st alle is > 1"
			vcfObject1 = vcfObject1._init_(line[1], line[3], alleles[0], "ins")
			
			#print "this vartype is " + str(vcfObject1.getVartype())

		if len(alleles[1]) == 1:
			vcfObject2 = vcfObject2._init_(line[1], line[3], alleles[1], "snp")
		elif len(alleles[1]) > 1:
			vcfObject2 = vcfObject2._init_(line[1], line[3], alleles[1], "ins")	

	#print "herro"	
	

	if vcfObject1.getVartype() == vcfObject2.getVartype() and vcfObject1.getStart() == vcfObject2.getStart():
		if vcfObject1.getReference() == vcfObject2.getReference():
			#print "Same Vartype and pos!"
			outfile.write(line[0] + "\t" + vcfObject1.getStart() + "\t" + line[2] + "\t" + vcfObject1.getReference() + "\t" + vcfObject1.getAllele() + "," + vcfObject2.getAllele() + "\t" + line[5] + "\t" + line[6] + "\t" + line[7]+ "\t" + line[8]+ "\t" + line[9] + "\t" + vcfObject1.getVartype() + "\n")
		else:
			if vcfObject1.getAllele() == vcfObject2.getAllele() and vcfObject1.getVartype() == 'del':
				#These are deletions of different length that start at the same position, comma separate them in reference
				outfile.write(line[0] + "\t" + vcfObject1.getStart() + "\t" + line[2] + "\t" + vcfObject1.getReference() + "," + vcfObject2.getReference() + "\t" + vcfObject1.getAllele() + "\t" + line[5] + "\t" + line[6] + "\t" + line[7]+ "\t" + line[8]+ "\t" + line[9] + "\t" + vcfObject1.getVartype() + "\n")
			else:
				print "something terible has happened: " + vcfObject1.getStart() + "\t" + vcfObject1.getReference() + "," + vcfObject2.getReference() + "\t" + vcfObject1.getAllele() + "\t" + vcfObject1.getAllele()

	else:
		g = line[9].split(":")
		genotype = g[0]
		if genotype == "2|1":
			g1 = "X|1"
			g2 = "1|X"
		elif genotype == "1|2":
			g1 = "1|X"
			g2 = "X|1"
		else:
			g1 = "1/X"
			g2 = "X/1"

		outfile.write(line[0] + "\t" + vcfObject1.getStart() + "\t" + line[2] + "\t" + vcfObject1.getReference() + "\t" + vcfObject1.getAllele() + "\t" + line[5] + "\t" + line[6] + "\t" + line[7]+ "\t" + line[8]+ "\t" + g1+ ":" + ":".join(g[1:]) + "\t" + vcfObject1.getVartype() + "\n")
		outfile.write(line[0] + "\t" + vcfObject2.getStart() + "\t" + line[2] + "\t" + vcfObject2.getReference() + "\t" + vcfObject2.getAllele() + "\t" + line[5] + "\t" + line[6] + "\t" + line[7]+ "\t" + line[8]+ "\t" + g2+ ":" + ":".join(g[1:]) + "\t" + vcfObject2.getVartype() + "\n")



def vcfTrim(position, ref, alt):
	# check the tail   
	end = 0 
	#print "orig_poz: " + position + " " + ref + " " + alt
	for i in range(1,min( len(ref),len(alt)) ):
		if ref[-i] == alt[-i]:
			end = end + 1
		else:
			break

	if end == 0:
		# if end is > -1 then we have no matching tail
		pass
	else:
		# we have to leave one common element (no dashes)
		if end == len(ref) or end == len(alt):
			if end == 1:
				pass
			else:
				end = end - 1
				alt_end = len(alt) - end
				alt = alt[:alt_end]
				ref_end = len(ref) - end
				ref = ref[:ref_end]
		else:
			alt_end = len(alt) - end
			alt = alt[:alt_end]
			ref_end = len(ref) - end
			ref = ref[:ref_end]

	start = 0
	for i in range(0,min( len(alt),len(ref)) ):
		if alt[i] == ref[i]:
			start += 1
		else:
			break
	#To be consistent with vcf files we will have to leave at least 1 element, no dashes		
	if start == len(alt) or start == len(ref):
		start = start - 1

	# trim the head 
	ref = ref[start:]
	
	alt = alt[start:]
	position = int(position) + start

	#get the varType
	if len(ref) == 1 and len(alt) == 1:
		varType = "snp"
	elif len(ref) == 1 and len(alt) > 1:
		varType = "ins"
	elif len(ref) > 1 and len(alt) == 1:
		varType = "del"
	else:
		varType = "sub"
		print "this is a sub" + str(position) + "\t" + ref + "\t" + alt + "\t"

	#print "end poz: " + str(position) + " " + ref + " " + alt

	vcfObject = vcf()	
	vcfObject = vcfObject._init_(str(position), ref, alt, varType)	
	
	return vcfObject



inFilename = str(sys.argv[1]) 
outfilename = str(sys.argv[2])

infile = gzip.open(inFilename, 'rb')
outfile = gzip.open(outfilename, 'wb')


print 'VCF fix start'
print datetime.datetime.now().time()


def main():
	count_temp = 0
	complex_count = 0
	complex_snp = 0
	simple_snp = 0
	simple_ins = 0
	simple_del = 0
	for inLine in infile:
		count_temp = count_temp + 1
		#print str(count_temp)
		if inLine.startswith("#"):
			print "header"
			outfile.write(inLine)
		else:
			count_temp = count_temp + 1
			line = inLine.strip().split("\t")
			if ',' in line[4] and len(line[4]) > 3:
				#print inLine
				fix_complex(line)
				complex_count = complex_count + 1
			elif ',' in line[4] and len(line[4]) == 3:
				#print inLine
				inLine = inLine.strip()
				outfile.write(inLine + "\t" + "snp" + "\n")
				complex_snp = complex_snp + 1
			else:
				inLine = inLine.strip()
				line = inLine.strip().split("\t")
				#get the varType
				if len(line[3]) == 1 and len(line[4]) == 1:
					varType = "snp"
					simple_snp = simple_snp + 1
				elif len(line[3]) == 1 and len(line[4]) > 1:
					varType = "ins"
					simple_ins = simple_ins + 1
				elif len(line[3]) > 1 and len(line[4]) == 1:
					varType = "del"
					simple_del = simple_del + 1
				else:
					print "what is this? " + inLine

				outfile.write(inLine + "\t" + varType + "\n")	



	print "Total variant: " + str(count_temp)
	print "complex variants: " + str(complex_count)
	print "complex snp " + str(complex_snp)
	print "simple snp " + str(simple_snp)
	print "simple ins " + str(simple_ins)
	print "simple del " + str(simple_del)

	infile.close()
	outfile.close()

if __name__ == '__main__': 
	main()