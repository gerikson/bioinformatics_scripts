'''
Script that creates bash file for each chromosoem and starts a process.
author: Galina Erikson
14 February 2014 
'''

import os, sys, datetime


#def main():
QSUB = "qsub -q stsi -M gerikson@scripps.edu -l mem=4G -l cput=24:00:00 -l walltime=24:00:00 "
original_directory = "/gpfs/group/stsi/data/wellderly/CG_data/CG_chromosome/master2014_01_30_unfiltered_chr"
output_directory = "/gpfs/group/stsi/data/gerikson/Wellderly/master2014_01_30_unfiltered_shuffled_chr"
jobdir = "/gpfs/home/gerikson/CNV_pipeline/Wellderly_scripts/"


i = 1
while i < 23:

	input_file = original_directory + str(i) + ".vcf.gz"
	output_file = output_directory + str(i) + ".vcf.gz"
	command = "python " + jobdir + "wellderly_shuffle_multiprocess.py " + input_file + " " + output_file
	jobfile = jobdir + "chr" + str(i) + "_dist_shuffle.job"
	outjob = open(jobfile, 'w')
	outjob.write("#!/bin/csh\n")			
	outjob.write("module load python-addons\n")
	outjob.write(command + "\n")
	outjob.write("touch "+"chr"+str(i)+'_dist_shuffle.job.res\n')				
	outjob.close()		
	execute = QSUB +'-e '+"chr"+str(i)+'_dist_shuffle.job.err -o '+"chr"+str(i)+'_dist_shuffle.job.out ' + jobfile
	#print execute
	sys.stdout.flush()
	clustnum = os.popen(execute, 'r')
	jobnum = clustnum.readline().strip()
	clustnum.close()
	i = i + 1  
	continue

'''
if __name__ == '__main__': 
	main()
'''
