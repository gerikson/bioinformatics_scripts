'''
Splits Complete Genomics masterVar files into chromosome specific masterVar
files when given an input file path and an output directory path.

e.g. >python masterVar_chr_split.py -i /path/to/masterVar.tsv.bz2 -o /path/to/output_dir/

Python package dependencies:
pandas, numpy
python 2.7 for argparse module
'''




import pandas as pd
import os, sys
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str,help='Specifies the input file, /path/to/CG_data/masterVar.tsv.bz2')
parser.add_argument('-o', type=str,help='Specifies the output directory, e.g. /path/to/CG_data/chromosome/')



def chr_split_mastervar(f_path, target_path):
    
    
    #Get header for masterVar
    header = os.popen('bzcat ' + f_path+ ' | head -100 | grep chromosome -n').readlines()    
    
    #Creating Reader object for iterating through NA12878 CG masterVar file
    skip_rows = int(header[0].split(":")[0]) -1
    mastervar_headings = os.popen('head -' + str(skip_rows) + f_path).readlines()
   
    #Creating pandas dataframe with chunksize 200,000 lines
    chunk = pd.read_table(f_path, chunksize=200000, sep="\t", skiprows=skip_rows,compression='bz2',dtype=object)
    chunk.columns = header[0].rstrip('\n').split(":")[1].split("\t")  #relabeling columns
    
    prev_chr = 'chr1'  #tracking chromosome position
    prev_target_write_file = None
    for mastervar in chunk:  #iterate through mastervar file
        for current_chrom,chr_df in mastervar.groupby(['chromosome']):  #split dataframe by chromosome for writing
            
            #check for increment to new chromosome
            if prev_chr != current_chrom:  
                os.system('bzip2 ' + prev_target_write_file) #compress last chromosome file
                prev_chr = current_chrom
            
            #specifying output file path and chromosome-specific name
            file_name = f_path.split("/")[-1].rstrip(".tsv.bz2") #getting file prefix 
            target_path = target_path.rstrip("/")+"/"  #ensuring target path ends with fwd slash
            write_target_file_path = target_path +file_name + "_" + current_chrom +".tsv"  #specify target directory and chrom file name
           
            
            #print write_target_file_path
            if len(os.popen('find '+ write_target_file_path + '').readlines()) == 0:  #checking for output file
                os.system('bzcat '+ f_path + '| head -' + str(skip_rows) + " > " +write_target_file_path) #writing header if no output file found
                chr_df.to_csv(write_target_file_path, sep="\t", index=False, mode='a')  #writing chromosome specific variants to output file
            else:  #Suppress header if target file found
                chr_df.to_csv(write_target_file_path, sep="\t", index=False, mode='a', header=False)  #writing chromosome specifc variants to output file w/o header
            
            prev_target_write_file = write_target_file_path  #increment to current write_target_file_path
                
            
    return 'complete'

opts = parser.parse_known_args()
f_path, target_path = opts[0].i, opts[0].o 
assert f_path.split(".")[-2:] == ['tsv','bz2'], "expecting masterVar input file suffix .tsv.bz2"

test = chr_split_mastervar(f_path, target_path)
if test == 'complete':
    print 'All chromosomes processed'





