#!/usr/bin/python
import  os, sys, math, time, pymongo,json,bson
'''
For testing purposes exit program with ctrl+c
'''
import signal
import sys
def signal_handler(signal, frame):
    print 'You pressed Ctrl+C!'
    sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)
#print 'Press Ctrl+C to exit'

HOME = "/gpfs/group/stsi/methods/annotation/sg-adviser/dev/"
SCHORK_HOME = "/gpfs/home/nschork/"

sys.path.append(SCHORK_HOME + "biopython-1.54/")
from collections import defaultdict

''''
Script that extracts annotation from mongoDB given a vcf file, adds genotype in notes column 
and the rest of the info from vcf as a last column
'''

varfilename = str(sys.argv[1])
CUR_CHROM = str(sys.argv[2])
resultfile = str(sys.argv[3])
nfile = str(sys.argv[4])

infile = open(varfilename, 'r')
outfile = open(resultfile, 'w' )
blah = open(nfile, 'w')

outfile.write("Haplotype	Chromosome	Begin	End	VarType	Reference	Allele	Notes	Gene	Gene_Type	Location	Distance	Coding_Impact	Protein_Pos	Original_AA	Allele_AA	Start~Stop_Dist	Prop_Cons_Affected_Upstream	Prop_Cons_Affected_Downstream	Trunc_Prediction	Conserved46way	Conserved46wayPlacental	Conserved46wayPrimates	ASW_minallele	CEU_minallele	CHB_minallele	CHD_minallele	GIH_minallele	JPT_minallele	LWK_minallele	MEX_minallele	MKK_minallele	TSI_minallele	YRI_minallele	1000GENOMES_AF	WELLDERLY_AF325	NHLBI	eQTL_genes	miRNA_BS_influenced	miRNA_BS_impact	miRNA_BS_direct	miRNA_BS_deltaG	miRNA_genomic	miRNA_folding_deltaG	miRNA_binding_deltaG	miRNA_top_targets_changed	Splice_Site_Pred	Splicing_Prediction(MaxENT)	ESE_sites	ESS_sites	Protein_Impact_Prediction(Polyphen)	Protein_Impact_Probability(Polyphen)	Protein_Impact_Prediction(SIFT)	Protein_Impact_Score(SIFT)	Protein_Domains	Protein_Domains_Impact(LogRE)	Protein_Impact_Prediction(Condel)	Protein_Impact_Probability(Condel)	TF_Binding_Sites	TFBS_deltaS	omimGene_ID~omimGene_association	Protein_Domain_Gene_Ontology	dbSNP_ID	HGMD_Variant~PubMedID	HGMD_Gene~disease_association	Genetic_Association_Database~PubMedID	PharmGKB_Database~Drug	Inheritance~Penetrance	Severity~Treatability	COSMIC_Variant~NumSamples	COSMIC_Gene~NumSamples	MSKCC_CancerGenes	Atlas_Oncology	Sanger_CancerGenes	Sanger_Germline_CancerGenes	Sanger_network-informed_CancerGenes~Pval	SegDup_Region	Gene_Symbol	DrugBank	Reactome_Pathway	Gene_Onotology	Disease_Ontology	ACMG_Score_Clinical~Disease_Entry~Explanation	ACMG_Score_Research~Disease_Entry~Explanation")
outfile.write("\n")
print "START: "
print time.asctime( time.localtime(time.time()) )

print "Chromosome is: "+CUR_CHROM

'''
Open db connection
'''
conn = pymongo.ReplicaSetConnection(host='stsia0808.cluster.net', replicaSet = 'rs0')
connstr = "conn.preannotation."+str(CUR_CHROM)
a =eval(connstr)
PreAnnotationDB = a
counter = 0

line = infile.readline()
found_var = 0
not_found = 0
#The name of the genotypes
genotypesHeader = ""
#genotypes usually starts in column 9 but we do the calculation later to make sure this is correct
genStart = 9
while line and True:	
    '''
    If this is are first lines, skip
    '''
    if line.startswith("##"):
	line = infile.readline()
        continue
    elif line.startswith('#CHROM') or line.startswith("#"):
	parts = line.strip().split("\t")
	for i in range(0, len(parts)):
		p = parts[i].lower()
		if p == 'format':
			genStart = i + 1
	#extract the genotypes to a String for future reference if needed
	for j in range(genStart, len(parts)):
		if j == len(parts):
			genotypesHeader = genotypesHeader + parts[j]
		else:
			genotypesHeader = genotypesHeader + parts[j] + "-"
	t = len(parts) - genStart
    	#outfile.write("Haplotype        Chromosome      Begin   End     VarType Reference       Allele  ")
	#outfile.write(genotypesHeader)
	#outfile.write("   Gene    Gene_Type       Location        Distance        Coding_Impact   Protein_Pos     Original_AA     Allele_AA       Start~Stop_Dist Prop_Cons_Affected_Upstream     Prop_Cons_Affected_Downstream   Trunc_Prediction        Conserved46way  Conserved46wayPlacental Conserved46wayPrimates  ASW_minallele   CEU_minallele   CHB_minallele   CHD_minallele   GIH_minallele   JPT_minallele   LWK_minallele   MEX_minallele   MKK_minallele   TSI_minallele   YRI_minallele   1000GENOMES_AF  WELLDERLY_AF325 NHLBI   eQTL_genes      miRNA_BS_influenced     miRNA_BS_impact miRNA_BS_direct miRNA_BS_deltaG miRNA_genomic   miRNA_folding_deltaG    miRNA_binding_deltaG    miRNA_top_targets_changed       Splice_Site_Pred        Splicing_Prediction(MaxENT)     ESE_sites       ESS_sites       Protein_Impact_Prediction(Polyphen)     Protein_Impact_Probability(Polyphen)    Protein_Impact_Prediction(SIFT) Protein_Impact_Score(SIFT)      Protein_Domains Protein_Domains_Impact(LogRE)   Protein_Impact_Prediction(Condel)       Protein_Impact_Probability(Condel)      TF_Binding_Sites        TFBS_deltaS     omimGene_ID~omimGene_association        Protein_Domain_Gene_Ontology    dbSNP_ID        HGMD_Variant~PubMedID   HGMD_Gene~disease_association   Genetic_Association_Database~PubMedID   PharmGKB_Database~Drug  Inheritance~Penetrance  Severity~Treatability   COSMIC_Variant~NumSamples       COSMIC_Gene~NumSamples  MSKCC_CancerGenes       Atlas_Oncology  Sanger_CancerGenes      Sanger_Germline_CancerGenes     Sanger_network-informed_CancerGenes~Pval        SegDup_Region   Gene_Symbol     DrugBank        Reactome_Pathway        Gene_Onotology  Disease_Ontology        ACMG_Score_Clinical~Disease_Entry~Explanation   ACMG_Score_Research~Disease_Entry~Explanation")
	#outfile.write("\n")
	print "Number of genotype is: " + str(t) 
	print "Genotypes are: " + genotypesHeader    
	line = infile.readline()
        continue
   	 
    line = line.strip().split()
    tline = line	
    if line[0].startswith("chr"):
	chrom = line[0]
    else:			
    	chrom = 'chr'+line[0]
    '''
    Make sure it's the same chromosome, if not connect to a different collection	
    '''	 	
    if chrom != CUR_CHROM:
	print "New chromosome!"
	CUR_CHROM = chrom
	connstr = "conn.preannotation."+str(CUR_CHROM)
	a =eval(connstr)
	PreAnnotationDB = a	
    
    '''
    -1 due to one based coordinates
    '''
    origbeg = int(line[1]) - 1
    origref = line[3]
    '''
    Extract the genotypes
    '''
    genotype = []
    som = ""
    for g in range(genStart, len(line)):
        gen = line[g].strip().split(":")
	genotype.append(gen[0])	

    obs = line[4].strip().split(',')
    for index in range(1, len(obs)+1):
	counter = counter + 1
	
	if counter % 100000 == 0:
		print "Variants: " + str(counter)
	'''
	Extract the genotypes
	'''
	
	obsCount = str(index)
	currentGen = []
	for g in genotype:
		som = ""
		if "/" in g:
			tempG = g.strip().split("/")
		else:
			tempG = g.strip().split("|")
		
		for s in tempG:
			if s == obsCount:
				som = som + "1"
			elif s == '0' or s == '.':
				som = som + s
			else:
				som = som + "X"
			continue
		currentGen.append(som)
        
        var = obs[index-1]
        ref = origref
        begpos = origbeg
        start = 0
        end = 0

                #if len(var) == len(ref) and len(var) > 1:
                    # handles messy cases like:
                    # ref =   'TGCCGTGGCGGCA'
                    # vars = ['TGCCGTGGCGGAA','TGCCGTGGTGGCA'] 
                    #  TGCCGTGGCGGCA
                    #             C
                    #  TGCCGTGGCGGAA
                    #             A
                    #  ************
                    #  TGCCGTGGCGGCA
                    #          C
                    #  TGCCGTGGTGGCA
                    #          T
                    #  ************                   
        for i in range(0,min( len(var),len(ref)) ):
                    if var[i] == ref[i]:
                        start += 1
                    else:
                        break
               
                # trim the head 
        ref = ref[start:]
        var = var[start:]
        begpos += start
	endpos = begpos + 1
        vartype = ''            
                # check the tail    
        for i in range(1,min( len(var),len(ref)) ):
                    if var[-i] == ref[-i]:
                        end -= 1
                    else:
                        break
                   
        if end == 0:
                    # if end is > -1 then we have no matching tail
                    pass
        else:
                    var = var[:end]
                    ref = ref[:end]

        if len(var) == 1 and len(ref) == 1:
                    endpos = begpos + 1
                    vartype = 'snp'
                    
        elif len(ref) == 0 and len(var) > 0:
                    begpos += 1
                    endpos = begpos
                    vartype = 'ins'
                    ref = '-'
        elif len(ref) > 0 and len(var) == 0:
                    vartype = 'del'
                    endpos = begpos + len(ref)
                    var = '-'
        elif len(ref) > len(var) and\
                        ( var[0] == '.' or\
                          var[0] == '-' ) and\
                          '<' not in ref: 
                    offset = len(var)
                    ref = ref[offset:]
                    begpos = begpos + offset
                    endpos = begpos + len(ref)
                    vartype = 'del'
                    var = '-'
        elif len(ref) > 0 and  len(var) > 0 and\
                    ref != var:
                    ref = ref.replace("-","").strip()[0:]
                    vartype = "delins"
                    endpos = begpos + len(ref)
        '''
	Query the db
	'''
	'''
	try:
        	variant_list = PreAnnotationDB.find_one({'chrom':CUR_CHROM,'beg':str(begpos), 'end':str(endpos)})
        	
        	if str(variant_list) == "None":
            		print "Variant not found! beg: " + str(begpos) + " beg: " + str(endpos) 
	    		line = infile.readline()
            		continue
        	else:
	    		#print "Found!"
			found_var = found_var + 1
			#gs = "-".join(genotype)
        		#outfile.write(str(var_found) + "\t" + CUR_CHROM + "\t" +str(begpos) + "\t" + str(endpos) + "\t" + str(vartype) + "\t" + str(ref) + "\t" + str(var)+ "\t" + str(gs) + "\t")    
            		outfile.write(str(var_found) + "\t" + CUR_CHROM + "\t" +str(begpos) + "\t" + str(endpos) + "\t" + str(vartype) + "\t" + str(ref) + "\t" + str(var)+ "\t")
		for var,val in variant_list.items():
                	
                	#Split the annotation at @@ and copy to file
              			
                	if (var == 'annot'):
                    		val = val.encode('ascii', 'ignore')
                    		annot_res = val.split("@@")
                    		s = len(annot_res)    
                    		outfile.write("\t".join(annot_res[0:s]))
            	
    	except:
		print "Exception"	
        '''
	variant_list = PreAnnotationDB.find_one({'chrom':CUR_CHROM,'beg':str(begpos), 'end':str(endpos)})
                
        if str(variant_list) == "None":
			not_found = not_found + 1
                        #print "Variant not found! beg: " + str(begpos) + " beg: " + str(endpos) 
                        blah.write(str(not_found) + "\t" + CUR_CHROM + "\t" +str(begpos) + "\t" + str(endpos) + "\t" + str(vartype) + "\t" + str(ref) + "\t" + str(var)) #+ "\t" + str(gs) + "\t")
			line = infile.readline()
                        continue
        else:
                        #print "Found!"
                        found_var = found_var + 1
                        gs = "-".join(currentGen)
                        outfile.write(str(found_var) + "\t" + CUR_CHROM + "\t" +str(begpos) + "\t" + str(endpos) + "\t" + str(vartype) + "\t" + str(ref) + "\t" + str(var)+ "\t" + str(gs) + "\t")    
                        #outfile.write(str(found_var) + "\t" + CUR_CHROM + "\t" +str(begpos) + "\t" + str(endpos) + "\t" + str(vartype) + "\t" + str(ref) + "\t" + str(var)+ "\t")
        		for var,val in variant_list.items():
                        	'''
                        	Split the annotation at @@ and copy to file
                        	'''     
                        	if (var == 'annot'):
                                	val = val.encode('ascii', 'ignore')
                                	annot_res = val.split("@@")
                                	s = len(annot_res)    
                                	outfile.write("\t".join(annot_res[0:s]))
			
			ts = "~".join(tline)
			outfile.write("\t" + ts)
			outfile.write("\n")
        		outfile.flush()
        		line = infile.readline()
        		continue   
        
conn.close() # must explicitly close the connection
infile.close()
outfile.close()

print "Variants not found: " + str(not_found)
print "Variants found: " + str(found_var)
print "END: ", time.time()
print time.asctime( time.localtime(time.time()) )

#sys.stdout.flush()
        
