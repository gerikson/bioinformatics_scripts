'''
Script that splits a parsed file by chromosome and starts cnv_genefeat.py on each chromosome 
Author: Galina Erikson
September 12th 2013
'''


import os, sys

varfilename = str(sys.argv[1])
HGTYPE = str(sys.argv[2])

ind1 = varfilename.rindex("/")

ORIGDIR = varfilename[0:ind1+1]
JOBDIR = varfilename[0:ind1+1] + "jobs/"
SUMDIR = varfilename[0:ind1+1] + "summary/"
RESDIR = varfilename[0:ind1+1] + "temp/"
TEMP_CHROMRES = RESDIR + "genefeat_res/"
TEMP_CHROMPATH = RESDIR + "genefeat_temp/"

QSUB = "qsub -q workq -M gerikson@scripps.edu -l mem=4G -l cput=24:00:00 -l walltime=24:00:00 "

#create the directories
allres_genefeat = [JOBDIR, SUMDIR, RESDIR, TEMP_CHROMPATH, TEMP_CHROMRES]
for res in allres_genefeat:
    if os.path.exists(res) != True:
        os.system("mkdir -p " + res)
    else:
        pass

ind = varfilename.rindex('/')
ind2 = varfilename.rindex('.')

outresultlist_genefeat = TEMP_CHROMRES + varfilename[ind+1:ind2] + "_genefeat_reslist.txt"
allres_genefeat = open(outresultlist_genefeat, 'w')

debugfile = varfilename[0:ind1+1] + "debug_file.txt"
debug = open(debugfile, 'w');

# conservation shell script that will copy and execute all of the qsub scripts for cnv_conservation.py 
conservationfile = varfilename[0:ind1+1] + "conservation.sh"
conservation = open(conservationfile, 'w');
conservation.write("#!/bin/csh\n")
                        


JOB_FILES = [] 

cur_chrom = ""
header = ""
invar = open(varfilename)
line = invar.readline()
while line != "":
    	if line.startswith("Haplotype"):
        	header = line.strip()
        	line = invar.readline()
        	continue

	orig_line = line.strip()
    	line = line.strip().split('\t')
	# is this a new chromosome
    	if line[1] != cur_chrom:
		
        	if cur_chrom != "":
            		outvar.close()
            			
            		resfile_genefeat = TEMP_CHROMRES + cur_chrom + "_genefeat.txt"
            		allres_genefeat.write(resfile_genefeat + "\n")

	    		#resfile_var = TEMP_CHROMRES + cur_chrom + "_var.txt"
            		#allres_var.write(resfile_var + "\n")	

			# start the cnv_genefeat.py on this chromosom
			command = "python /gpfs/home/gerikson/CNV_pipeline/scripts/cnv_genefeat.py "  + ORIGDIR + " " + outvarname + " " +  resfile_genefeat + " hg19"
			jobfile = JOBDIR + cur_chrom + "_dist_annotate.job"
			outjob = open(jobfile, 'w')
            		outjob.write("#!/bin/csh\n")			
			outjob.write("module load python-addons\n")
            		outjob.write(command + "\n")
			outjob.write("touch "+SUMDIR+cur_chrom+'_genefeat_init.job.res\n')			
			outjob.close()
			
			execute = QSUB +'-e '+SUMDIR+cur_chrom+'_genefeat_init.job.err -o '+SUMDIR+cur_chrom+'_genefeat_init.job.out ' + jobfile
			print execute
            		sys.stdout.flush()
			
			clustnum = os.popen(execute, 'r')
            		jobnum = clustnum.readline().strip()
            		clustnum.close()
            		ind0 = jobfile.rindex("/")
            		#ind = jobnum.index(".")
			debug.write(jobnum + "\n")
			#debug.close()
            		JOB_FILES.append(SUMDIR+jobfile[ind0+1:] + ".res")
		
			'''
			Create the job file for the next process cnv_conservation.py that will start upon the succesfull completion of the genefeat	
			'''
        		
			resfile_conservation = TEMP_CHROMRES + cur_chrom + "_conservation.txt"
			command_conservation = "python /gpfs/home/gerikson/CNV_pipeline/scripts/cnv_conservation.py " + resfile_genefeat + " " + cur_chrom  + " " + resfile_conservation + " hg19"
			jobfile_conservation = JOBDIR + cur_chrom + "_conservation.job"			
			outjob_conservation = open(jobfile_conservation, 'w')
			outjob_conservation.write("#!/bin/csh\n")                    
                        outjob_conservation.write("module load python-addons\n")
                        outjob_conservation.write(command_conservation + "\n")
			outjob_conservation.write("touch "+SUMDIR+cur_chrom+'_conservation_init.job.res\n')
			outjob_conservation.close()
			
			execute = QSUB + " -W depends=afterok:" + jobnum + ' -e '+SUMDIR+cur_chrom+'_conservation_init.job.err -o '+SUMDIR+cur_chrom+'_conservation_init.job.out ' + jobfile_conservation
                       
			'''
			copy this to the conservation file that will be executed once we escape the while loop
			'''
			conservation.write(execute + "\n")
			#conservation.close()	
			

			
			print execute
                        sys.stdout.flush()

                        clustnum = os.popen(execute, 'r')
                        jobnum_conservation = clustnum.readline().strip()
                        clustnum.close()
                        ind0 = jobfile.rindex("/")
                        debug.write("cnv_conservation.py " + jobnum_conservation + "\n")
                        JOB_FILES.append(SUMDIR+jobfile_conservation[ind0+1:] + ".res")
			

        	#temp_chr = line[1]
		#print "currentc chrom is " + line[1]
		'''
		copying the original parsed file to temp/genefeat_temp 
		'''
		cur_chrom = line[1]
        	outvarname = TEMP_CHROMPATH + cur_chrom + "_temp.txt"
        	
		outvar = open(outvarname, 'w')
        	outvar.write(orig_line +"\n")
	else:
		outvar.write(orig_line + "\n")
	line = invar.readline()
	outvar.flush()


debug.close()
conservation.close()
'''
Execute the conservation script
'''
os.system("sh " + conservationfile)
outvar.close()

	
