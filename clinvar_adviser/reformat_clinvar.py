import os, sys, urllib, time, mmap
from SOAPpy import WSDL

class DbSNPConverter:
    
    URL = 'https://mutalyzer.nl/services/?wsdl'
    NCBI_CHROMS = {'NC_000006.11': 'chr6', 'NC_000005.9': 'chr5', 'NC_000023.10': 'chrX', 'NC_000013.10': 'chr13', 'NC_000003.11': 'chr3', 'NC_000024.9': 'chrY', 'NC_000009.11': 'chr9', 'NC_000022.10': 'chr22', 'NC_000011.9': 'chr11', 'NC_000001.10': 'chr1', 'NC_000004.11': 'chr4', 'NC_000008.10': 'chr8', 'NC_000002.11': 'chr2', 'NC_000007.13': 'chr7', 'NC_000016.9': 'chr16', 'NC_000018.9': 'chr18', 'NC_000010.10': 'chr10', 'NC_000019.9': 'chr19', 'NC_000014.8': 'chr14', 'NC_000015.9': 'chr15', 'NC_000012.11': 'chr12', 'NC_000017.10': 'chr17', 'NC_000021.8': 'chr21', 'NC_000020.10': 'chr20'}
    
    def __init__(self):
        self.api = WSDL.Proxy(self.URL)
    
    def convert(self,chrom,begpos,endpos,rsid):
        
        try:
            res = self.api.getdbSNPDescriptions(rs_id=rsid)
        except:
            return None
        #print "converting rsid: ", rsid
        #print "res: ", res
        for recordlist in res:
            #print "recordlist: ", recordlist
            for record in recordlist:
                #print "record: ", record
                id = record.split(":")[0]
                if id in self.NCBI_CHROMS and "g" in record:
                    checkchrom = self.NCBI_CHROMS[id]
                    if chrom == None:
                        chrom = checkchrom
                    elif checkchrom != chrom:
                        continue
                    mut = record.split(":")[1]
                    return extract_varentry(chrom,begpos,endpos,mut)
        return None
    
    def convert_chrom(self,chromid):
        if chromid in self.NCBI_CHROMS:
            return self.NCBI_CHROMS[chromid]
        else:
            return None

def get_chromseq(chrom):
    
    #chromfile = CHROM_PREFIX + chrom + ".mmap"
    chromfile = CHROM_PREFIX + chrom + ".fa"
    inchrom = open(chromfile,'rb')
    CHROMSEQ = mmap.mmap(inchrom.fileno(),0,access=mmap.ACCESS_READ)
    inchrom.close()
    return CHROMSEQ

def check_varentry(varentry):
    global ALLCHROMSEQ
    chrom, begpos, endpos, vartype, ref, obs = varentry.split("_")
    
    if chrom == "chrM":
        return True
    
    if chrom not in ALLCHROMSEQ:
        chromseq = get_chromseq(chrom)
        ALLCHROMSEQ[chrom] = chromseq
        
    if vartype == "snp" or vartype == "del" or vartype == "delins":
        if ref != ALLCHROMSEQ[chrom][int(begpos):int(endpos)]:
            return False
    return True

def extract_varentry(chrom,begpos,endpos,mut,rsid=None):
    
    varentry = None
    """
    if mut == "g.89017961G>A":
        print "PROCESSING g.89017961G>A", mut
    """
    if chrom == "chrM":
        return None
    elif chrom == None:
        if rsid == "N/A" or rsid == None:
            return None
        return dbSNPConverter.convert(chrom,begpos,endpos,rsid)
    try:
        if '>' in mut:
            vartype = "snp"
            mut = mut.split()[0].replace(",","")
            ref = mut[-3]
            alt = mut[-1]
            begpos = str(int(mut.split(">")[0][:-1].split("g.")[-1].replace("(","").replace(")",""))-1)
            endpos = str(int(begpos)+1)
            
        elif "del" in mut and "ins" in mut:
            vartype = "delins"
            ref = mut.split("del")[-1].split("ins")[0]
            alt = mut.split("ins")[-1].split()[0]
            begpos = str(int(mut.split("_")[0].split("g.")[-1].split("del")[0])-1)
            endpos = str(int(begpos)+len(ref))
        elif "del" in mut:
            vartype = "del"
            ref = mut.split("del")[-1]
            alt = "-"
            begpos = str(int(mut.split("_")[0].split("g.")[-1].split("del")[0])-1)
            endpos = str(int(begpos)+len(ref))
        elif "ins" in mut:
            vartype = "ins"
            ref = "-"
            begpos = str(int(mut.split("_")[0].split("g.")[-1]))
            alt = mut.split("ins")[-1]
            endpos = begpos
        else:
            if rsid == "N/A" or rsid == None:
                return None
            else:
                return dbSNPConverter.convert(chrom,begpos,endpos,rsid)
    except Exception, e:
        print "VARENTRY_EXTRACT_ERROR: ", chrom, rsid, e
        sys.stdout.flush()
        return None
    
    CHECKVAR = set(['A','C','G','T','-'])
    
    if len(set(list(ref)).difference(CHECKVAR)) > 0 or len(set(list(alt)).difference(CHECKVAR)) > 0:
        if rsid == "N/A" or rsid == None:
            return None
        else:
            return dbSNPConverter.convert(chrom,begpos,endpos,rsid)
    
    varentry = chrom + "_" + begpos + "_" + endpos + "_" + vartype + "_" + ref + "_" + alt
    
    return varentry

#CHROM_PREFIX = "/home/tools/ccs/Human_Annotation/hg19/chroms/"
CHROM_PREFIX = "/gpfs/group/torkamani/phpham/Human_Annotation/hg19/chroms/"
version = "2013.3"

dbSNPConverter = DbSNPConverter()

infilename= "clinvar_hg19.gff"
#infilename = "temp.gff"
outfilename = "clinvar_hg19." + version + ".varentry_new.txt"
errfilename = "clinvar_hg19." + version + ".varentry_new.errs"


infile = open(infilename)
outfile = open(outfilename,'w')
outfile.write("#varentry\tclinvar_disease~clinvar_prediction~clinvar_id~omim_id\n")

outerr = open(errfilename,'w')

TRACKER = {}
ALLCHROMSEQ = {}

ERR = 0
COUNT = 0
for line in infile:
    COUNT += 1
    if COUNT % 100 == 0:
        print COUNT, time.time()
        sys.stdout.flush()
    if line.startswith("#"):
        continue
    
    line = line.strip().split("\t")
    
    chrom = line[0]
    origbegpos = line[3]
    origendpos = line[4]
    
    entry = line[-1].split("; ")
    
    
    try:
        allentries = dict([(x.split("=")) for x in entry if len(x.split("=")) == 2])
        line[2] = urllib.unquote(allentries['hgvs'])
    except KeyError,e:
        line[2] = None
    
    omim = allentries['omim']
    disease = allentries['clinvar_DiseaseName'].strip()
    pred = allentries['clinvar_ClinicalSignificance']
    clinvar_id = allentries['clinvar_Accession']
    gene = allentries['hgnc'].strip()
    rsid = allentries['rsid']
    
    diseasegene_entry = disease + "(" + gene + ")"
    finalentry = diseasegene_entry + '~' + pred + '~' +clinvar_id + '~' + omim
    
    if rsid == "N/A" and "feature" in allentries and "(rs" in allentries['feature']:
        feature = allentries['feature']
        ind = feature.index("(rs")
        feature = feature[ind+1:]
        ind0 = feature.index(")")
        rsid = feature[:ind0]
    
    #print line[2]
    mut = ":".join(line[2].split(":")[1:]).split(";")[0]
    chromid = line[2].split(":")[0]
    converted_chrom = dbSNPConverter.convert_chrom(chromid)
    if converted_chrom != chrom:
        chrom = None
    
    varentry = extract_varentry(chrom,origbegpos,origendpos,mut,rsid)
    """
    print "PROCESSING in main loop!"
    print "original line[2]: ", line[2]
    print "new mut: ", mut
    print "varentry result: ", varentry
    print "rsid: ", rsid
    """
    if varentry == None:
        outerr.write("\t".join(line) + "\n")
        ERR += 1
        continue
    else:
        isverified = check_varentry(varentry)
        if not isverified:
            print "ERROR IN VARENTRY: ", varentry, line
            sys.stdout.flush()
            outerr.write("VARENTRY_ERROR" + "\t" + "\t".join(line) + "\n")
            ERR += 1
            continue
    
    finalentry = finalentry.decode("ascii","ignore")
    if varentry in TRACKER:
        TRACKER[varentry].append(finalentry)
    else:
        TRACKER[varentry] = [finalentry]

for varentry, allentries in TRACKER.iteritems():
        outfile.write(varentry + '\t' + "///".join(allentries) + "\n")
outfile.flush()
outfile.close()
infile.close()
outerr.flush()
outerr.close()
print "NUMERR: ", ERR





