"""
Script that extracts the AF from a 1000Genomes CNV final aggregate file and gives each CNV a new ID

Takes 2 arguments: the input filename and the output file
author: gerikson
date: March 11, 2014
"""

import sys

inFilename = str(sys.argv[1])
outfilename = str(sys.argv[2])

inFile = open(inFilename)
outfile = open(outfilename, 'w')

def max(arr):
	maxim = 0
	for num in arr:
		if int(num) > maxim:
			maxim = num
	#print str(maxim)
	return maxim

def EUR (counter, line):
	if line[6] == "-":
			CEU = 0
	else:
		CEU = int(line[6])

	if line[9] == "-":
		FIN = 0
	else:
		FIN = int(line[9])
	if line[12] == "-":
		GBR = 0
	else:
		if "///" in line[12]:
			temp_l = line[12].split("///")
			GBR = int(temp_l[0])
		else:
			GBR = int(line[12])
	if line[15] == "-":
		IBS = 0
	else:
		IBS = int(line[15])
	if line[18] == "-":
		TSI = 0
	else:
		TSI = int(line[18])
	total_allele = CEU + FIN + GBR + IBS + TSI
	total_population = 341
	outfile.write("\t".join(line[:5]) + "\tEUR_CNVnator_" + str(counter) + "\t" + str(total_allele) + "\t" + "341" + "\t" + "\t".join(line[5:]) + "\n")	

def ASN (counter, line):
	if line[6] == "-":
		CHB = 0
	else:
		CHB = int(line[6])

	if line[9] == "-":
		CHS = 0
	else:
		if "///" in line[9]:
			temp_l = line[9].split("///")
			CHS = int(temp_l[0])
		else:
			CHS = int(line[9])		

	if line[12] == "-":
		JPT = 0
	else:
		JPT = int(line[12])
	total_allele = CHB + CHS + JPT
	total_population = 341
	outfile.write("\t".join(line[:5]) + "\tASN_CNVnator_" + str(counter) + "\t" + str(total_allele) + "\t" + "253" + "\t" + "\t".join(line[5:14]) + "\n")	

def AMR(counter, line):
	if line[6] == "-":
		CLM = 0
	else:
		if "///" in line[6]:
			temp_l = line[6].split("///")
			CLM = int(temp_l[0])
		else:
			CLM = int(line[6])

	if line[9] == "-":
		MXL = 0
	else:
		if "///" in line[9]:
			temp_l = line[9].split("///")
			MXL = int(temp_l[0])
		else:
			MXL = int(line[9])		

	if line[12] == "-":
		PUR = 0
	else:
		PUR = int(line[12])
	total_allele = CLM + MXL + PUR
	total_population = 341
	#outfile.write("\t".join(line[:5]) + "\tAMR_CNVnator_" + str(counter) + "\t" + str(total_allele) + "\t" + "158" + "\t" + "\t".join(line[5:14]) + "\n")	
	outfile.write("\t".join(line[:5]) + "\tAFR_CNVnator_" + str(counter) + "\t" + str(total_allele) + "\t" + "210" + "\t" + "\t".join(line[5:14]) + "\n")			

def aggregate(counter, line):
	if line[6] == "-":
		AFR = 0
	else:
		if "///" in line[6]:
			temp_l = line[6].split("///")
			AFR = max(temp_l)
		else:
			AFR = int(line[6])

	if line[9] == "-":
		AMR = 0
	else:
		if "///" in line[9]:
			temp_l = line[9].split("///")
			AMR = max(temp_l)
		else:
			AMR = int(line[9])		

	if line[12] == "-":
		ASN = 0
	else:
		if "///" in line[12]:
			temp_l = line[12].split("///")
			ASN = max(temp_l)
		else:
			ASN = int(line[12])

	if line[15] == "-":
		EUR = 0
	else:
		if "///" in line[15]:
			temp_l = line[15].split("///")
			EUR = max(temp_l)
		else:
			EUR = int(line[15])
	
	total_allele_freq = (float(AFR) + float(AMR) + float(ASN) + float(EUR))/962.0
	AFR_AF = float(AFR)/210.0
	AMR_AF = float(AMR)/158.0
	ASN_AF = float(ASN)/253.0
	EUR_AF = float(EUR)/341.0

	outfile.write("\t".join(line[:5]) + "\t1000G_CNVnator962_" + str(counter) + "\t" + str(total_allele_freq) + "\t" + line[5] + "\t" + str(AFR_AF) + "\t" + line[8]+ "\t" + str(AMR_AF) + "\t" + line[11] + "\t" + str(ASN_AF) + "\t" + line[14] + "\t" + str(EUR_AF)+"\n")			


if __name__ == "__main__":
	outfile.write("Chromosome\tBegin\tEnd\tVarType\tFilter\tAggregate_ID\tAggregate_AF\tAFR_ID\tAFR_AF\tAMR_ID\tAMR_AF\tASN_ID\tASN_AF\tEUR_ID\tEUR_AF\n")
	#outfile.write("Chromosome\tBegin\tEnd\tVarType\tFilter\tEUR_ID\tEUR_Sample_count\tEUR_Total_Sample\tCEU_ID\tCEU_Sample_Count\tCEU_Total_Sample\tFIN_ID\tFIN_Sample_Count\tFIN_Total_Sample\tGBR_ID\tGBR_Sample_Count\tGBR_Total_Sample\tIBS_ID\tIBS_Sample_Count\tTSI_Total_Sample\tTSI_ID\tTSI_Sample_Count\tTSI_Total_Sample\n")
	#outfile.write("Chromosome\tBegin\tEnd\tVarType\tFilter\tASN_ID\tASN_Sample_count\tASN_Total_Sample\tCHB_ID\tCHB_Sample_Count\tCHB_Total_Sample\tCHS_ID\tCHS_Sample_Count\tCHS_Total_Sample\tJPT_ID\tJPT_Sample_Count\tJPT_Total_Sample\n")
	#outfile.write("Chromosome\tBegin\tEnd\tVarType\tFilter\tAFR_ID\tAFR_Sample_count\tAFR_Total_Sample\tASW_ID\tASW_Sample_Count\tASW_Total_Sample\tLWK_ID\tLWK_Sample_Count\tLWK_Total_Sample\tYRI_ID\tYRI_Sample_Count\tYRI_Total_Sample\n")
	counter = 0
	counter = 0
	for line in inFile:
		counter = counter + 1
		line = line.strip().split("\t")
		aggregate(counter, line)
		#AMR(counter, line)
		#EUR(counter, line)