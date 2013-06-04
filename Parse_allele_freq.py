#!/usr/bin/python
import  os, sys 

''''
Script that extracts the allele frequency and parses the file for the wellderly for updating the tabix format SG Adviser.
 Allele frequency files are located /gpfs/home/gerikson/wellderly/alleleFreq
'''
varfilename = str(sys.argv[1])
resultfile = str(sys.argv[2])

infile = open(varfilename, 'r')
outfile = open(resultfile, 'w' )

line = infile.readline()

while line:
    #totalLines += 1
    line = line.strip().replace('"','').split()
    chrom = 'chr'+line[0]
    ''''
    -1 due to one based coordinates
    '''
    origbeg = int(line[1]) - 1
    tempRef = line[4].split(':')
    origref = tempRef[0]
    
    for index in range(5, len(line)):
        TempAltAllele = line[index].split(':')
        var = TempAltAllele[0]
        AlleleFreq = TempAltAllele[1]
        ref = origref
       # var = obs[index-1]
        begpos = origbeg
        ref = origref
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

        outfile.write(chrom + "\t" + str(begpos) + "\t" + str(endpos) + "\t" + vartype + '\t' 
                                + ref.strip() + "\t" + var.strip() + '\t' + AlleleFreq + "\n")
      
	'''	
	else:
                    skips += 1
                    print line    
                    #it += 1
	'''
    line = infile.readline()
    #it += 1
infile.close()
outfile.close()
        
