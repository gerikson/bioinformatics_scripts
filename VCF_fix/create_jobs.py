import os, sys, datetime

inFilename = str(sys.argv[1]) 

infile = open(inFilename, 'rb')


print 'VCF fix start'
print datetime.datetime.now().time()

counter = 0
QSUB = "qsub -q workq -M gerikson@scripps.edu -l mem=4G -l cput=24:00:00 -l walltime=24:00:00 "

for sample in infile:    

    sample = sample.strip()
    counter = counter + 1
    s = sample.split('/')
    resfile = s[7][:-6]
    print "resfile: " + resfile
    outputfile = "/gpfs/group/stsi/data/gerikson/Wellderly_Illumina_vcf_fix/fixed/" + resfile + "fix.vcf.gz"
    command = "python  /gpfs/home/gerikson/scripts/Wellderly_scripts/fix_vcf.py " + sample + " " + outputfile
    print "command " + command 
    jobfile = "/gpfs/group/stsi/data/gerikson/Wellderly_Illumina_vcf_fix/jobs/" + s[7] + ".job"         
    outjob = open(jobfile, 'w')
    outjob.write("#!/bin/csh\n")                    

    outjob.write(command + "\n")
    outjob.close()  
    execute = QSUB + ' -e '+ '/gpfs/group/stsi/data/gerikson/Wellderly_Illumina_vcf_fix/jobs/' + s[7] +'.job.err -o' + '/gpfs/group/stsi/data/gerikson/Wellderly_Illumina_vcf_fix/jobs/'+ s[7] + '.job.out ' + jobfile
    print execute 

    sys.stdout.flush()
    clustnum = os.popen(execute, 'r')
    jobnum = clustnum.readline().strip()
    clustnum.close()

infile.close()