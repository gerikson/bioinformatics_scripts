'''
Script that splits a parsed file by chromosome and starts cnv_genefeat.py on each chromosome 
Author: Galina Erikson
September 12th 2013
'''


import os, sys, datetime

varfilename = str(sys.argv[1])
HGTYPE = str(sys.argv[2])
#DIRECTORY = str(sys.argv[3])

ind1 = varfilename.rindex("/")

ORIGDIR = varfilename[0:ind1+1]
JOBDIR = varfilename[0:ind1+1] + "jobs/"
SUMDIR = varfilename[0:ind1+1] + "summary/"
RESDIR = varfilename[0:ind1+1] + "temp/"
TEMP_CHROMRES = RESDIR + "genefeat_res/"
TEMP_CHROMPATH = RESDIR + "genefeat_temp/"
CONSERVATION = RESDIR + "conservation/"
OVERLAP = RESDIR + "overlap/"
OVERLAP_W = RESDIR +"overlap_w/" 

QSUB = "qsub -q stsi -M gerikson@scripps.edu -l mem=4G -l cput=24:00:00 -l walltime=24:00:00 "

#create the directories
allres_genefeat = [JOBDIR, SUMDIR, RESDIR, TEMP_CHROMPATH, TEMP_CHROMRES, CONSERVATION, OVERLAP, OVERLAP_W]
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

'''
New qsub that will start the post process only after all jobs were completed
'''

#qsub_post_process = QSUB

jobarray = ""

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
			#print execute
            		sys.stdout.flush()
		
			clustnum = os.popen(execute, 'r')
            		jobnum = clustnum.readline().strip()
            		clustnum.close()
			'''
			copy all of the jobs to the final qsub that will start post_process.py	
			'''
			#qsub_post_process = qsub_post_process + jobnum + ","
            		ind0 = jobfile.rindex("/")
            		#ind = jobnum.index(".")
			debug.write(jobnum + "\n")
			#debug.close()
            		JOB_FILES.append(SUMDIR+jobfile[ind0+1:] + ".res")
		
			'''
			Create the job file for the next process cnv_conservation.py that will start upon the succesfull completion of the genefeat	
			'''
        		
			resfile_conservation = CONSERVATION + cur_chrom + "_conservation.txt"
			command_conservation = "python /gpfs/home/gerikson/CNV_pipeline/scripts/cnv_conservation.py " + resfile_genefeat + " " + cur_chrom  + " " + resfile_conservation + " hg19"
			jobfile_conservation = JOBDIR + cur_chrom + "_conservation.job"			
			outjob_conservation = open(jobfile_conservation, 'w')
			outjob_conservation.write("#!/bin/csh\n")                    
                        outjob_conservation.write("module load python-addons\n")
                        outjob_conservation.write(command_conservation + "\n")
			outjob_conservation.write("touch "+SUMDIR+cur_chrom+'_conservation_init.job.res\n')
			outjob_conservation.close()
			
			execute = QSUB + " -W depend=afterok:" + jobnum + ' -e '+SUMDIR+cur_chrom+'_conservation_init.job.err -o '+SUMDIR+cur_chrom+'_conservation_init.job.out ' + jobfile_conservation
                        
 
			'''
			copy this to the conservation file that will be executed once we escape the while loop
			'''
			conservation.write(execute + "\n")
			#conservation.close()	
			
			#print execute
                        sys.stdout.flush()

                        clustnum_conservation = os.popen(execute, 'r')
                        jobnum_conservation = clustnum_conservation.readline().strip()
                        clustnum_conservation.close()
                        ind0 = jobfile_conservation.rindex("/")
                        debug.write("cnv_conservation.py " + jobnum_conservation + "\n")
                        JOB_FILES.append(SUMDIR+jobfile_conservation[ind0+1:] + ".res")
			                        
                        #qsub_post_process = qsub_post_process + jobnum_conservation + ","
			jobarray= jobarray + jobnum_conservation + ":"

			'''
			Execute the cnv_overlap.py script
			'''	
			
			resources = "/gpfs/home/gerikson/CNV_pipeline/resources/cnv_reference_db"	
			resfile_overlap = OVERLAP + cur_chrom + "_overlap.txt"
                        command_overlap = "python /gpfs/home/gerikson/CNV_pipeline/scripts/cnv_overlap.py " + outvarname + " " + resfile_overlap + " " + cur_chrom  + " hg19 " + resources +" cnv"
                        jobfile_overlap = JOBDIR + cur_chrom + "_overlap.job"
                        outjob_overlap = open(jobfile_overlap, 'w')
                        outjob_overlap.write("#!/bin/csh\n")
                        outjob_overlap.write("module load python-addons\n")
                        outjob_overlap.write(command_overlap + "\n")
                        outjob_overlap.write("touch "+SUMDIR+cur_chrom+'_overlap_init.job.res\n')
                        outjob_overlap.close()
                        execute = QSUB + ' -e '+SUMDIR+cur_chrom+'_overlap_init.job.err -o '+SUMDIR+cur_chrom+'_overlap_init.job.out ' + jobfile_overlap
                        #print execute
                        sys.stdout.flush()
                        clustnum_overlap = os.popen(execute, 'r')
                        jobnum_overlap = clustnum_overlap.readline().strip()
                        clustnum_overlap.close()
                        ind0 = jobfile_overlap.rindex("/")
                        debug.write("cnv_overlap.py " + jobnum_overlap + "\n")
                        JOB_FILES.append(SUMDIR+jobfile_overlap[ind0+1:] + ".res")
			#qsub_post_process = qsub_post_process + jobnum_overlap + ","
			jobarray = jobarray + jobnum_overlap + ":"
	
			'''
                        Execute the cnv_overlap.py script, wellderly mode
                        '''

                        resources_w = "/gpfs/home/gerikson/CNV_pipeline/resources/wellderly_cnv_reference_db" 
                        resfile_overlap_w = OVERLAP_W + cur_chrom + "_overlapW.txt"
                        command_overlap_w = "python /gpfs/home/gerikson/CNV_pipeline/scripts/cnv_overlap.py " + outvarname + " " + resfile_overlap_w + " " + cur_chrom  + " hg19 " + resources_w +" wellderly"
                        jobfile_overlap_w = JOBDIR + cur_chrom + "_overlapW.job"
                        outjob_overlap_w = open(jobfile_overlap_w, 'w')
                        outjob_overlap_w.write("#!/bin/csh\n")
                        outjob_overlap_w.write("module load python-addons\n")
                        outjob_overlap_w.write(command_overlap_w + "\n")
                        outjob_overlap_w.write("touch "+SUMDIR+cur_chrom+'_overlapW_init.job.res\n')
                        outjob_overlap_w.close()                        
			executeW = QSUB + ' -e '+SUMDIR+cur_chrom+'_overlapW_init.job.err -o '+SUMDIR+cur_chrom+'_overlapW_init.job.out ' + jobfile_overlap_w
                        #print execute
                        sys.stdout.flush()
                        clustnum_overlapW = os.popen(executeW, 'r')
                        jobnum_overlapW = clustnum_overlapW.readline().strip()
                        clustnum_overlapW.close()
                        ind0 = jobfile_overlap_w.rindex("/")
                        debug.write("cnv_overlap.py wellderly mode " + jobnum_overlapW + "\n")
                        JOB_FILES.append(SUMDIR+jobfile_overlap_w[ind0+1:] + ".res")
			#qsub_post_process = qsub_post_process + jobnum_overlapW+ ","
			jobarray = jobarray + jobnum_overlapW + ":"
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


conservation.close()

# remove the last character of the jobarray string (the ':)
jobarray = jobarray[:-1]
qsub_post_process = QSUB + " -W depend=afterok:" + jobarray

print qsub_post_process
print datetime.datetime.now().time()

'''
Create shell script that will start the qsub_post_process 
'''

command_post_process = "python /gpfs/home/gerikson/CNV_pipeline/post_process.py " + ORIGDIR + " " + varfilename[ind+1:ind2] 
jobfile_post_process = ORIGDIR + "post_process.job"
outjob_post_process = open(jobfile_post_process, 'w')
outjob_post_process.write("#!/bin/csh\n")
outjob_post_process.write("module load python-addons\n")
outjob_post_process.write(command_post_process + "\n")
outjob_post_process.close()
execute_post_process = qsub_post_process + ' -e '+ORIGDIR+'post_process.job.err -o '+ORIGDIR+'post_process.job.out ' + jobfile_post_process
print execute_post_process
sys.stdout.flush()
clustnum_post_process = os.popen(execute_post_process, 'r')
jobnum_post_process = clustnum_post_process.readline().strip()
clustnum_post_process.close()
ind0 = jobfile_post_process.rindex("/")
debug.write("Post process " + jobnum_post_process + "\n")
JOB_FILES.append(SUMDIR+jobfile_post_process[ind0+1:] + ".res")
print 'POST PROCESS'
print datetime.datetime.now().time()
debug.close()

'''
Execute the conservation script
'''
#os.system("sh " + conservationfile)
outvar.close()

	
