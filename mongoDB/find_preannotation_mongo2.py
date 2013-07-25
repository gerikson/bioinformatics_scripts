import os, sys,math, time, pymongo,json,bson

#GROUP_PATH = os.environ['GROUP_PATH']
#HG_PATH = os.environ['HG_PATH']
#HOME = os.environ['SGHOME'] + "/"
#SCHORK_HOME = os.environ['SCHORK_HOME']
HOME = "/gpfs/group/stsi/methods/annotation/sg-adviser/dev/"
SCHORK_HOME = "/gpfs/home/nschork/"
#HOME = "/ddn/people/sgadvise/"

sys.path.append(SCHORK_HOME + "biopython-1.54/")
from collections import defaultdict

def pr(CUR_CHROM):
    print CUR_CHROM
    return CUR_CHROM    
 
def reorder_list(optfile):
    print "REORDERING"
    
    intemp = open(optfile)
    oldlist = intemp.read().strip().split("\n")
    intemp.close()
    
    newlist = []
    
    inlist = open(KEY_INDEX_FILE)
    for listline in inlist.readlines():
        listline = listline.strip().split("\t")
        if listline[0] in oldlist:
            for num in listline[1:]:
                newlist.append(int(num))
    
    newlist.sort()
    
    return newlist


print "START: "
print time.asctime( time.localtime(time.time()) )
varfilename = str(sys.argv[1])

CUR_CHROM = str(sys.argv[2])
print "Chromosome is: "+CUR_CHROM
dbsnp_resfilename =str(sys.argv[3])
resfilename = str(sys.argv[4])
#ANNOT_OPTIONS = str(sys.argv[5])
HGTYPE = str(sys.argv[5])
SIMPLE_MAP = False

try:
    if str(sys.argv[7]).upper() == 'TRUE':
        SIMPLE_MAP = True
except:
    pass

'''
if HGTYPE == "hg18":
    CHROM_PREFIX = HG_PATH + "hg18/"
elif HGTYPE == "hg19":    
    CHROM_PREFIX = HG_PATH + "hg19/chroms/"
    KEY_INDEX_FILE= HOME + "UPDATE/KEY_INDEX_FILE_hg19.txt"

#HEADER_INDICES = reorder_list(ANNOT_OPTIONS)
'''

print "CHECKING FOR PRE-EXISTING VARIANTS: ", varfilename

outres = open(resfilename, 'w')
outres_dbsnp = open(dbsnp_resfilename, 'w')


all_vars = []
invar = open(varfilename)
line = invar.readline()
while line != "":
    if line.startswith("Haplotype"):
        line = invar.readline()
        continue
    line = line.strip().split("\t")
    temp = []
    temp.append(line[0])
    temp.append(line[1])
    temp.append(int(line[2]))
    temp.append(int(line[3]))
    try: 
        for k in range(4,8):
            temp.append(line[k])
        all_vars.append(temp)
    except:
        temp.append(line[4:])
        all_vars.append(temp)
    line = invar.readline()
invar.close()

var_it = 0
tracker = 0
content_it = 0
var_res = []
var_annot = []
content = []

'''
Querying mongo db, for now development chr1, change later
move the connection out of the while loop otherwise it re-connects for every variant
'''
#conn = pymongo.Connection(host='hadoop00-adm.cluster.net')#10.128.128.60')#host='10.128.143.255')
#conn = pymongo.ReplicaSetConnection('hadoop00-adm.cluster.net',replicaSet='rs0') #10.128.128.60')#host='10.128.143.255')
#conn = pymongo.Connection(host='stsia0808.cluster.net')
conn = pymongo.ReplicaSetConnection(host='stsia0808.cluster.net', replicaSet = 'rs0')
connstr = "conn.development."+str(CUR_CHROM)
a =eval(connstr)#"conn.development."+str(CUR_CHROM))
PreAnnotationDB = a

try:

    while var_it < len(all_vars):
        
        chrom = None #all_vars[var_it][0]
        begpos = None #int(all_vars[var_it][1])
        endpos = None #int(all_vars[var_it][2])
        vartype = None #all_vars[var_it][4]
        ref = None #all_vars[var_it][5]
        obs = None #all_vars[var_it][6]
    
        if SIMPLE_MAP:
            chrom = all_vars[var_it][1]
            begpos = int(all_vars[var_it][2])
            endpos = int(all_vars[var_it][3])
            #vartype = all_vars[var_it][4]
            #ref = all_vars[var_it][5]
            #obs = all_vars[var_it][6]
        else:
            chrom = all_vars[var_it][1]
            #begpos = int(all_vars[var_it][2])
           # endpos = int(all_vars[var_it][3])
            beg = str(all_vars[var_it][2])
            end = str(all_vars[var_it][3])
            vartype = all_vars[var_it][4]
            ref = all_vars[var_it][5]
            obs = all_vars[var_it][6]
         
                
        
        """
        # comparison to check querying only using the indexed keys and doing clientside filtering
        variant_list2 = PreAnnotationDB.find({'chrom':CUR_CHROM,'beg':beg, 'end':end})
        variant_list = None
        for i in variant_list2:
            if i['ref'] == ref and\
                i['alt'] == obs and\
                 i['type'] == vartype:
                 
                variant_list = i
                break
        """
        #variant_list = PreAnnotationDB.find_one({'chrom':CUR_CHROM,'beg':beg, 'end':end}) #,'type':vartype ,'ref':ref,'alt':obs})
    	variant_list = PreAnnotationDB.find_one({'chrom':CUR_CHROM,'beg':beg, 'end':end,'type':vartype ,'ref':ref,'alt':obs})
        """
        If no variant was found, skip to the next element
        """
        if str(variant_list) == "None":
	    #print "Not found " + time.asctime( time.localtime(time.time()) ) 	
            for entry in all_vars[var_it]:
                outres.write(str(entry) + "\t")
            outres.write("\n")
            var_it += 1
            FOUND = False
            continue
        else:
            FOUND = True
    
        if(FOUND):
	    #print "found " + time.asctime( time.localtime(time.time()) )	
            #write stuff at the begin
            outres_dbsnp.write(str(chrom) + "\t" +str(beg) + "\t" + str(end) + "\t" + str(vartype) + "\t" + str(ref) + "\t" + str(obs)+ "\t")    
            for var,val in variant_list.items():
                '''
                if the variable is not the annotation, coppy to file, else split the annotation at @@ and copy to file
                '''
                if (var == 'annot'):
                    val = val.encode('ascii', 'ignore')
                    annot_res = val.split("@@")
                    s = len(annot_res)    
                    outres_dbsnp.write("\t".join(annot_res[0:s]))
   	        '''	 
            for it in HEADER_INDICES:
                outres_dbsnp.write('\t' + annot_res[it])
    	    '''
	
            outres_dbsnp.write("\n")
            outres_dbsnp.flush()
            var_it += 1
            continue   
    
        if FOUND != True:
            for entry in all_vars[var_it]:
                outres.write(str(entry) + "\t")
            outres.write("\n")
            outres.flush()
        
        var_it += 1 
        continue
except:
    # safety to ensure that we can reach the conn.close() statement to ensure we don't have dangling sockets on mongo's end
    pass

conn.close() # must explicitly close the connection
outres.close()
outres_dbsnp.close()

print "END: ", time.time()
print time.asctime( time.localtime(time.time()) )

sys.stdout.flush()
