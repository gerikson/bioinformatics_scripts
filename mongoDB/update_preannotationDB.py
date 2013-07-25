import os, sys,math, time, pysam, pymongo,json,bson

HG_PATH = os.environ['HG_PATH']
SCHORK_HOME = os.environ['SCHORK_HOME']


sys.path.append(SCHORK_HOME + "biopython-1.54/")
from collections import defaultdict

print "START: "
print time.asctime( time.localtime(time.time()) )
"""
Takes 3 arguments: 1.filename_for_DB_update(*.novel_variant) 
2.file to extract the non_ascii variants
3.file to extract the weird chromosome annotation
"""
varfilename = str(sys.argv[1])
CUR_CHROM = 'chr1'
non_ascii = str(sys.argv[2])
weird_chrom = str(sys.argv[3])
print "Chromosome is: "+CUR_CHROM

print "CHECKING FOR PRE-EXISTING VARIANTS: ", varfilename

inserted_variants = 0
skiped_var = 0
non_ascii_var = 0
invar = open(varfilename)
lines = invar.readlines()
ascii = open(non_ascii, "w")
w_chrom = open(weird_chrom, "w")
line = invar.readline()

'''
Querying mongo db, for now development
'''
#con = pymongo.Connection(host='10.128.128.60')#host='10.128.143.255')
con = pymongo.Connection(host='stsia0808.cluster.net')
# a =eval("con.development."+str(CUR_CHROM))
#a = eval("con.gerikson."+str(CUR_CHROM))
a = eval("con.development.test")
PreAnnotationDB = a

while line:
    """
    If this is the first line skip
    """    
    lin = line.strip().split("\t")
    
    if (lin[0]=="Haplotype"):
        line = invar.readline()
        continue     
    # make sure this is a chromosome string and if it's the same as the current chromosome
    if (lin[1].find('chr') == -1):
        w_chrom.write(line+"\t")
        line = invar.readline()
        continue    
    elif (lin[1] != CUR_CHROM):
        print "new chromosome " + lin[1]
        CUR_CHROM = lin[1]
        a = eval("con.gerikson."+str(CUR_CHROM))
        PreAnnotationDB = a
    
    beg = lin[2]
    end = lin[3]
    vartype = lin[4]
    ref = lin[5]
    obs = lin[6]
    annot = "-"
    for int in lin[8:]:
        annot=annot+"@@"+int
 
    '''
    verify is it's a varlid annotation, mongo doesn't accept non usii carachters
    '''   
    try:
        annot.decode('ascii')
    except UnicodeDecodeError:
        #print "it was not a ascii-encoded unicode string"
        #print annot
        ascii.write(str(line)+"\t")
        non_ascii_var = non_ascii_var + 1
        line = invar.readline()
        continue
    
    
    '''
    Update the database with the new variant
    '''
    inserted_variants = inserted_variants + 1
        
    post = {'chr':CUR_CHROM,'beg':beg, 'end':end, 'type':vartype, 'ref':ref, 'alt':obs, 'annot':annot}
 
    #print "inserting brg = " + beg
    if (inserted_variants%100000 == 0):
        print inserted_variants
        print time.asctime( time.localtime(time.time()) )
    PreAnnotationDB.insert(post)
    line = invar.readline()

print "END: ", time.asctime( time.localtime(time.time()) )
print "Inserted var = " + str(inserted_variants)
print "skiped var = " + str(skiped_var)
print "non ascii var found = " + str(non_ascii_var)

w_chrom.close()
ascii.close()
invar.close()
sys.stdout.flush()


