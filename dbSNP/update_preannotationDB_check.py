import os, sys,math, time, pysam, pymongo,json,bson

HG_PATH = os.environ['HG_PATH']
SCHORK_HOME = os.environ['SCHORK_HOME']


sys.path.append(SCHORK_HOME + "biopython-1.54/")
from collections import defaultdict

print "START: "
print time.asctime( time.localtime(time.time()) )
"""
takes to arguments:
1.file to update the db
2.file to extract non ascii variants
"""
varfilename = str(sys.argv[1])
CUR_CHROM ='chr1' 
non_ascii=str(sys.argv[2])

print "Chromosome is: "+CUR_CHROM

print "CHECKING FOR PRE-EXISTING VARIANTS: ", varfilename

inserted_variants = 0
skiped_var = 0
non_ascii_var = 0
invar = open(varfilename)
line = invar.readline()
ascii = open(non_ascii, "w")

'''
Querying mongo db, for now development
'''
con = pymongo.Connection(host='10.128.128.60')#host='10.128.143.255')
#a =eval("con.development."+str(CUR_CHROM))
a = eval("con.gerikson."+str(CUR_CHROM))
PreAnnotationDB = a

while line:
    lin = line.strip().split("\t")
    
    """
    If this is the first line skip
    """    
    if (lin[0]=="Haplotype"):
        continue

    # make sure this is a chromosome string and if it's the same as the current chromosome
    if (lin[1] != CUR_CHROM and lin[1].find('chr') != -1):
        print "new chromosome " + lin[1]
        CUR_CHROM = lin[1]
        a = eval("con.gerikson."+str(CUR_CHROM))
        PreAnnotationDB = a
    
    beg = lin[2]
    end = lin[3]
    vartype = lin[4]
    ref = lin[5]
    obs = lin[6]
    annot = "@@".join(lin[8:])

    '''
    verify is it's a valid annotation, mongo doesn't accept non usii carachters
    '''   
    try:
        annot.decode('ascii')
    except UnicodeDecodeError:
       # print "it was not a ascii-encoded unicode string"
       # print annot
        ascii.write(str(line)+"\n")
        non_ascii_var = non_ascii_var + 1
        line = invar.readline() 
#        continue  
   
    """
    First check if the variant is already present in the database, this step might be skipped if we get a file 
    with only variants that weren't previously found in the preannotationDB
    """
    variant_list = PreAnnotationDB.find_one({'beg':beg, 'end':end,'type':vartype ,'ref':ref,'alt':obs})
    
    """
    If no variant was found, insert the variant into the db
    """
    if str(variant_list) == "None":
        inserted_variants = inserted_variants + 1
        '''
        Variant not present in the preannotation db, update the database with the new variant
        '''
      #  print "Var not found, insert begin = " + beg        
        post = {'chr':CUR_CHROM,'beg':beg, 'end':end, 'type':vartype, 'ref':ref, 'alt':obs, 'annot':annot}
        PreAnnotationDB.insert(post)
        if (inserted_variants%1000 == 0):
            print inserted_variants
            print time.asctime( time.localtime(time.time()) )
       # continue
    else:
         skiped_var = skiped_var + 1
         if (skiped_var%500 == 0):
            print skiped_var
            print time.asctime( time.localtime(time.time()) )
            #continue
    line = invar.readline()

print "END: ", time.asctime( time.localtime(time.time()) )
print "Inserted var = " + str(inserted_variants)
print "skiped var = " + str(skiped_var)
print "non ascii var found = " + str(non_ascii_var)

invar.close()
sys.stdout.flush()

