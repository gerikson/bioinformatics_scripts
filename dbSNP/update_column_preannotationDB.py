import os, sys,math, time, pysam, pymongo,json,bson

HG_PATH = os.environ['HG_PATH']
SCHORK_HOME = os.environ['SCHORK_HOME']


sys.path.append(SCHORK_HOME + "biopython-1.54/")
from collections import defaultdict

print "START: "
print time.asctime( time.localtime(time.time()) )

"""
Takes 3 arguments: 
1.filename_for_DB_update - in format: chr begin end vartype ref alt comment(with rsID) 
2.The column name you wish to update
"""

varfilename = str(sys.argv[1])
CUR_CHROM = 'chr1'
columnName = str(sys.argv[2])

print "Chromosome is: "+CUR_CHROM

print "CHECKING FOR PRE-EXISTING VARIANTS: ", varfilename

inserted_variants = 0
skiped_var = 0
non_ascii_var = 0
invar = open(varfilename)

header = ('Gene',	'Gene_Type',	'Location',	'Distance',	'Coding_Impact', 'Protein_Pos',	'Original_AA',	'Allele_AA',	'Start~Stop_Dist',	'Prop_Cons_Affected_Upstream',	'Prop_Cons_Affected_Downstream',	'Trunc_Prediction',	'Conserved46way',	'Conserved46wayPlacental',	'Conserved46wayPrimates',	'ASW_minallele',	'CEU_minallele	CHB_minallele',	'CHD_minallele',	'GIH_minallele',	'JPT_minallele',	'LWK_minallele',	'MEX_minallele',	'MKK_minallele',	'TSI_minallele',	'YRI_minallele',	'1000GENOMES_AF',	'WELLDERLY_AF325',	'NHLBI',	'eQTL_genes',	'miRNA_BS_influenced',	'miRNA_BS_impact',	'miRNA_BS_direct',	'miRNA_BS_deltaG',	'miRNA_genomic',	'miRNA_folding_deltaG',	'miRNA_binding_deltaG',	'miRNA_top_targets_changed',	'Splice_Site_Pred',	'Splicing_Prediction(MaxENT)',	'ESE_sites',	'ESS_sites',	'Protein_Impact_Prediction(Polyphen)',	'Protein_Impact_Probability(Polyphen)',	'Protein_Impact_Prediction(SIFT)',	'Protein_Impact_Score(SIFT)',	'Protein_Domains',	'Protein_Domains_Impact(LogRE)',	'Protein_Impact_Prediction(Condel)',	'Protein_Impact_Probability(Condel)',	'TF_Binding_Sites',	'TFBS_deltaS',	'omimGene_ID~omimGene_association',	'Protein_Domain_Gene_Ontology',	'dbSNP_ID',	'HGMD_Variant~PubMedID',	'HGMD_Gene~disease_association',	'Genetic_Association_Database~PubMedID',	'PharmGKB_Database~Drug',	'Inheritance~Penetrance',	'Severity~Treatability',	'COSMIC_Variant~NumSamples',	'COSMIC_Gene~NumSamples', 'MSKCC_CancerGenes',	'Atlas_Oncology',	'Sanger_CancerGenes',	'Sanger_Germline_CancerGenes',	'Sanger_network-informed_CancerGenes~Pval',	'SegDup_Region',	'Gene_Symbol',	'DrugBank',	'Reactome_Pathway',	'Gene_Onotology',	'Disease_Ontology')
columnNumber = 0

'''
Extract the column number that needs updating
'''

for index, c in enumerate(header):
	if columnName == c:
		columnNumber = index
		print "Column that needs updating is: " + str(columnNumber)
		


'''
Querying mongo db, for now development
'''

con = pymongo.ReplicaSetConnection(host='stsia0808.cluster.net', replicaSet = 'rs0')
constr = "con.development.test"
#constr ="con.preannotation."+str(CUR_CHROM)
a = eval(constr)
PreAnnotationDB = a

line = invar.readline()
while line:
   # print line	
    """
    If this is the first line skip
    """    
    lin = line.strip().split("\t")
    
    if (lin[1]=="Chromosome"):
        line = invar.readline()
        continue     
    # make sure this is a chromosome string and if it's the same as the current chromosome
    
    if (lin[1].find('chr') == -1):
       # w_chrom.write(line+"\t")
        line = invar.readline()
        continue    
    elif (lin[1] != CUR_CHROM):
        print "new chromosome " + lin[1]
        CUR_CHROM = lin[1]
        a = eval("con.preannotation."+str(CUR_CHROM))
        PreAnnotationDB = a
    
    beg = lin[2]
    end = lin[3]
    vartype = lin[4]
    ref = lin[5]
    obs = lin[6]
    comment = lin[7]
    annot = ""

    '''
    verify is it's a valid annotation, mongo doesn't accept non usii carachters
    '''   
    try:
        annot.decode('ascii')
    except UnicodeDecodeError:
        #print "it was not a ascii-encoded unicode string"
        #print annot
        #ascii.write(str(line)+"\t")
        non_ascii_var = non_ascii_var + 1
        line = invar.readline()
        continue
 	
    variant_list = PreAnnotationDB.find({'beg':beg, 'end':end,'type':vartype ,'ref':ref,'alt':obs})
   
    if str(variant_list) == "None":
        skiped_var = skiped_var + 1
        if (skiped_var%500 == 0):
           print skiped_var
           print time.asctime( time.localtime(time.time()) )
        """
        Insert the variants (for now)
		"""        
        annot = "@@".join(lin[8:])
        
        
        post = {'chr':CUR_CHROM,'beg':beg, 'end':end, 'type':vartype, 'ref':ref, 'alt':obs, 'annot':annot}
        PreAnnotationDB.insert(post)
         
		line = invar.readline()
        continue	 
    else:
    	for var,val in variant_list.items():
            '''
            if the variable is not the annotation, coppy to file, else split the annotation at @@ and copy to file
            '''
            if (var == 'annot'):
                annot_res = val.split("@@")
                '''
                Change the columnNumber with new value 
                '''
                for index, s in enumerate(annot_res):
			if columnNumber == index:
                		annot = annot + comment + "@@"
			else:
				annot = annot + s + "@@"
                	#annot_res[columnNumber] = comment
                	'''
			Temp, insert '-' in Gene
			'''
			#print "Comment is " + str(comment)
			#annot_res[columnNumber] = comment
			#if annot_res.index[x] == columnNumber
                
                """
                Insert new value back into the database
                """	
                '''
                for int in annot_res:
        			annot=annot+"@@"+int
                '''
		print "Annotation is " + str(annot)
        	'''
        	Variant not present in the preannotation db, update the database with the new variant
        	'''      
        	post = {'chr':CUR_CHROM,'beg':beg, 'end':end, 'type':vartype, 'ref':ref, 'alt':obs}
		PreAnnotationDB.update(post, {'chr':CUR_CHROM,'beg':beg, 'end':end, 'type':vartype, 'ref':ref, 'alt':obs, 'annot': annot})
		if (inserted_variants%1000 == 0):
            		print inserted_variants
            		print time.asctime( time.localtime(time.time()) )

        	inserted_variants = inserted_variants + 1
	
    	line = invar.readline()
        continue
print "END: ", time.asctime( time.localtime(time.time()) )
print "Inserted var = " + str(inserted_variants)
print "skiped var = " + str(skiped_var)

con.close()
invar.close()
sys.stdout.flush()


