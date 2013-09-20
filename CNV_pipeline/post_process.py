#Post annotation process
#combines the files from each step into one single summary file
#Author: Galina Erikson
#September 18 2013

import os, sys, datetime

directory = str(sys.argv[1])
original_file = str(sys.argv[2])

print 'POST PROCESS START'
print datetime.datetime.now().time()

OUTP_FILENAME = directory + "/SUMMARY_"+original_file
CONSERVATIONDIR = directory + "/temp/conservation/"
OVERLAP = directory + "/temp/overlap/"
Weld_OVERLAP = directory + "/temp/overlap_w/"
sum_conservation = CONSERVATIONDIR + "SUMMARY_cons.txt"
sum_overlap = OVERLAP + "SUMMARY_overlap.txt"
weld_sum_overlap = Weld_OVERLAP + "SUMMARY_overlap_w.txt"

outp=open(OUTP_FILENAME, "w")
outp.write("Haplotype" + '\t' + "Chromosome" + '\t' + "Begin" + '\t' + "End" + '\t' + "VarType" + '\t' + "Reference" + '\t' + "Allele" + '\t' + "Notes" + '\t' + "Gene ID" + '\t' + "Gene_Type" + '\t' + "Location" + '\t' + "Distance" + '\t' + "Coding_Impact" + '\t' + "Protein_Pos" + '\t' + "Original_AA" + '\t' + "Allele_AA\tProp_Cons_Affected_Inside\tProp_Cons_Affected_Outside\tExonic_Bases_Inside_CNV\tExonic_Bases_Outside_CNV\tKnown_CNV(80%_Overlap_ID)\tKnown_CNV(80%_Overlap_AF_Gain)\tKnown_CNV(80%_Overlap_AF_Loss)\tKnown_CNV(100%_Overlap_ID)\tKnown_CNV(100%_Overlap_AF_Gain)\tKnown_CNV(100%_Overlap_AF_Loss)\tKnown_Wellderly(80%_Overlap_ID)\tKnown_Wellderly(80%_Overlap_AF_Gain)\tKnown_Wellderly(80%_Overlap_AF_Loss)\tKnown_Wellderly(100%_Overlap_ID)\tKnown_Wellderly(100%_Overlap_AF_Gain)\tKnown_Wellderly_CNV(100%_Overlap_AF_Loss)\n")

files_conservation = []

'''
put filenames in an array, sort the names and concatenate each separate step into one file
'''

command = "cat "
for file in os.listdir(CONSERVATIONDIR):
	files_conservation.append(file)
files_conservation.sort()
for file in files_conservation:
	command = command + " " + CONSERVATIONDIR  +  file

command = command + " >> " + sum_conservation
os.system(command)
sys.stdout.flush()

command = "cat "
files_overlap = []
for file in os.listdir(OVERLAP):
        files_overlap.append(file) 
files_overlap.sort()
for file in files_overlap:
	command = command + " " +  OVERLAP + file
command = command + " >> " + sum_overlap
os.system(command)
sys.stdout.flush()

command = "cat "
files_overlap_w = []
for file in os.listdir(Weld_OVERLAP):
        files_overlap_w.append(file)
files_overlap_w.sort()
for file in files_overlap_w:
        command = command + " " +  Weld_OVERLAP  +  file
command = command + " >> " + weld_sum_overlap
os.system(command)
sys.stdout.flush()


'''
Open files, go line by line and concatenate results if the line number is correct
'''
conserv = open(sum_conservation, 'r')
overlap = open(sum_overlap, 'r')
weld_overlap = open(weld_sum_overlap, 'r')

inLine = conserv.readline()
line_overlap = overlap.readline()
line_overlap_w = weld_overlap.readline()
while inLine.strip() != '':
	
	inLin = inLine.strip().split('\t')
	l_overlap = line_overlap.strip().split('\t')
	l_overlap_w = line_overlap_w.strip().split('\t')	
	# compare if it's the same variant
	#if (inLin[1] == l_overlap[1] && inLin[2] == l_overla[2] && inLin[3] == l_overlap[3] && inLin[4] == l_overlap[4]: #&& inLin[1] == l_overlap_w[1] && inLin[2] == l_overlap_w[2] && inLin[3] == l_overlap_w[3] && inLin[4] == l_overlap_w[4]):		
	if inLin[1:4] == l_overlap[1:4] and inLin[1:4] == l_overlap_w[1:4]:
		gs = ""
		gs = "\t".join(inLin)
		for i in l_overlap[19:]:
			gs = gs + "\t" + l_overlap[i]
		for j in l_overlap_w[19:]:
			gs = gs + "\t" +l_overlap_w[i]
		gs = gs + "\n"	
		outp.write(gs)
	else:
		print "Wrong line" + inLine
	inLine = conserv.readline()
	line_overlap = overlap.readline()
	line_overlap_w = weld_overlap.readline()

outp.close()
conserv.close()
overlap.close()
weld_overlap.close()

print 'POST PROCESS FINISHED'
print datetime.datetime.now().time()
