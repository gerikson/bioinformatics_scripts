'''
Script that retrieves a retrieves a region from Entrez
author: gerikson
date 14 February 2014
'''

from Bio import Entrez, SeqIO

#open the file with your Entrez gene IDs
input_file = open("path/to/to/the/genelist")
out_handle = open("example.txt", "w")
Entrez.email='gerikson@scripps.edu'

line = input_file.readline()

#this is a loop that goes through every single line of your file
while line != "":
	# Assuming each line it of the format Entrez gene ID\tBegin\tEnd
	line = line.strip().split('\t')
	handle = Entrez.efetch(db="nucleotide", id=, rettype=line[0], strand=1, seq_start=line[1], seq_stop=line[2], retmode='text')
	record = SeqIO.parse (handle, "fasta")
	SeqIO.write(record, out_handle, "fasta")
	line = input_file.readline()
    continue

in_handle.close()
out_handle.close()
input_file.close()