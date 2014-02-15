'''
Script that will insert clinvar into mongodb abd verifies the integrity of the database
Author: Galina Erikson
Date: February 7th 
'''

# -*- coding: utf-8 -*-
import os, sys, math, time, pymongo, json, bson, datetime, string
from bson.objectid import ObjectId
import traceback

HG_PATH = os.environ['HG_PATH']
SCHORK_HOME = os.environ['SCHORK_HOME']


sys.path.append(SCHORK_HOME + "biopython-1.54/")
from collections import defaultdict

print "START: "
print time.asctime( time.localtime(time.time()) )

'''
varfilename = str(sys.argv[1])
not_found = str(sys.argv[2])
original_line = str(sys.argv[3])
bad_var = str(sys.argv[4])
w_gene = str(sys.argv[5])
dup_entries_db = str(sys.argv[6])
clinvar_d_lenght = str(sys.argv[7])
'''

CUR_CHROM = 'chr1'
invar = open("clinvar_sorted.txt")
n_found = open("not_found.txt", "w")
orig_file = open("original_line.txt", "w")
diffrent_lenght = open("different_lenght.txt", "w")
wrong_gene = open("wrong_gene.txt", 'w')
dup_entries = open("duplicated_entries.txt", 'w')
clinvar_diff_length = open("clinvar_d_lenght.txt", 'w')

print 'START'
print datetime.datetime.now().time()

required = set(['Haplotype','Chromosome','Begin','End','VarType','Reference','Allele','Notes'])
headers = ['Chromosome','Begin','End','VarType','Reference','Allele','Gene','Gene_Type','Location','Distance','Coding_Impact', 'Protein_Pos','Original_AA','Allele_AA','Start~Stop_Dist','Prop_Cons_Affected_Upstream','Prop_Cons_Affected_Downstream','Trunc_Prediction','Conserved46way','Conserved46wayPlacental','Conserved46wayPrimates','ASW_minallele','CEU_minallele','CHB_minallele','CHD_minallele','GIH_minallele','JPT_minallele','LWK_minallele','MEX_minallele','MKK_minallele','TSI_minallele','YRI_minallele','1000GENOMES_AF','WELLDERLY_AF325','NHLBI','eQTL_genes','miRNA_BS_influenced','miRNA_BS_impact','miRNA_BS_direct','miRNA_BS_deltaG','miRNA_genomic','miRNA_folding_deltaG','miRNA_binding_deltaG','miRNA_top_targets_changed','Splice_Site_Pred','Splicing_Prediction(MaxENT)','ESE_sites','ESS_sites','Protein_Impact_Prediction(Polyphen)','Protein_Impact_Probability(Polyphen)','Protein_Impact_Prediction(SIFT)','Protein_Impact_Score(SIFT)','Protein_Domains','Protein_Domains_Impact(LogRE)','Protein_Impact_Prediction(Condel)','Protein_Impact_Probability(Condel)','TF_Binding_Sites','TFBS_deltaS','omimGene_ID~omimGene_association','Protein_Domain_Gene_Ontology','dbSNP_ID','HGMD_Variant~PubMedID','HGMD_Gene~disease_association','Genetic_Association_Database~PubMedID','PharmGKB_Database~Drug','Inheritance~Penetrance','Severity~Treatability','COSMIC_Variant~NumSamples','COSMIC_Gene~NumSamples', 'MSKCC_CancerGenes','Atlas_Oncology','Sanger_CancerGenes','Sanger_Germline_CancerGenes','Sanger_network-informed_CancerGenes~Pval','SegDup_Region','Gene_Symbol','DrugBank','Reactome_Pathway','Gene_Onotology','Disease_Ontology','ACMG_Score_Clinical~Disease_Entry~Explanation','ACMG_Score_Research~Disease_Entry~Explanation']

orig_file.write('\t'.join(headers) + "\n")
diffrent_lenght.write('\t'.join(headers) + "\n")

line = invar.readline().strip()
clinvar_variants = 0
inserted_variants = 0
duplicated_entries = 0
wrong_genes = 0
different_lenght = 0

con = pymongo.ReplicaSetConnection(host="stsia0808.cluster.net", replicaSet = 'rs0')
constr ="con.preannotation."+str(CUR_CHROM)
a = eval(constr)
PreAnnotationDB = a


while line:
    clinvar_variants = clinvar_variants + 1
    #counter that will identify if we find more then one copy of the variant
    count = 0
    lin = line.strip().split("\t")
    beg = lin[1]
    end = lin[2]
    vartype = lin[3]
    ref = lin[4]
    obs = lin[5]
    clinvar = lin[6]
    gene = lin[7].strip()
    #print "gene = " + gene 
    if (lin[0] != CUR_CHROM):
        print "new chromosome " + lin[0]
        CUR_CHROM = lin[0]
        a = eval("con.preannotation."+str(CUR_CHROM))
        PreAnnotationDB = a

    variant_list = PreAnnotationDB.find({'chrom':CUR_CHROM,'beg':beg, 'end':end,'type':vartype ,'ref':ref,'alt':obs}) 

    if str(variant_list) == "None":
        print "SKIP:"+line
        skiped_var = skiped_var + 1
        if (skiped_var%500 == 0):
           print skiped_var
           print time.asctime( time.localtime(time.time()) ) 
        n_found.write(line)
        line = invar.readline().strip()
        continue     
    else:
        for varlist in variant_list:
            count = count + 1   

            #print "Found variants"
            annotation = varlist['annot']
            #annot = annotation.split("@@")
            annotation = filter(lambda x: x in string.printable, annotation)
            annot = annotation.split("@@")
            orig_file.write(CUR_CHROM + "\t" + beg + "\t" + end + "\t" + vartype + "\t" + ref + "\t" + obs + "\t".join(annot) + "\n")
            found = False
            correct_lenght = True 
            #verify that the light is of correct lenght
            if len(annot) == 76:
                #verify that the gene wher ethe clinvar variant is located is the same exact as the gene from mongo
                gene_mongo = annot[69].split("///")
                for g in gene_mongo:
                    #g = g.strip()
                    if g == gene:
                        found = True

            else:
                #print "Annotation found in db has different lenght:  " + str(len(annot))
                #print annotation
                if len(annot) == 77 and clinvar == annot[57]:
                    print "good"
                else:
                    correct_lenght = False
                    different_lenght = different_lenght + 1
                    diffrent_lenght.write(CUR_CHROM + "\t" + beg + "\t" + end + "\t" + vartype + "\t" + ref + "\t" + obs + "\t".join(annot) + "\n")
                    clinvar_diff_length.write(line + "\n")

            modified_array = ""

            if found:
                modified_array = "@@".join(annot[:57]) + "@@" + clinvar + "@@" + "@@".join(annot[57:])
                varlist['annot'] = modified_array
                PreAnnotationDB.save(varlist)
                inserted_variants = inserted_variants + 1
                #print "variant inserted"
                if (inserted_variants%1000 == 0):
                            print inserted_variants
                            print time.asctime( time.localtime(time.time()) )
            elif correct_lenght :
                #print "Not the same gene! DB gene " + annot[69]
                #print "Clinvar gene: " + gene
                wrong_genes = wrong_genes + 1  
                wrong_gene.write(line + "\n")

        #verify to see if we have duplicated values
        if count > 1:
            #print "More then one variant in db" + str(count)
            #print CUR_CHROM + "\t" + beg + "\t" + end + "\t" + vartype + "\t" + ref + "\t" + obs + "\n"
            duplicated_entries = duplicated_entries + 1
            dup_entries.write(line + "\n")
        line = invar.readline().strip()
        continue

con.close()
invar.close()
n_found.close()
orig_file.close()
diffrent_lenght.close()
wrong_gene.close()
dup_entries.close()
clinvar_diff_length.close()

print "total clinvar variants queried " + str(clinvar_variants)
print "total lines updated "+str(inserted_variants)
print "total duplicated lines "+str(duplicated_entries)
print "total number of entries were the genes in db weren't the same as clinvar: " + str(wrong_genes)
print "total number of entries were the genes in db had a different_lenght: " + str(different_lenght)
print 'END'
print datetime.datetime.now().time()
