#!/usr/bin/python

from sga_stdlib import sga_base
from Bio.Seq import Seq
import argparse,traceback,sys,os,subprocess


class validator(sga_base):
    
    def resolveVariant(self,ref,alt,start,chrom):
        """
         ref is the raw data ref allele (less offset / padding etc)
         alt is the raw data alt allele (less offset / padding etc)
         start is the 0-based starting index of the reference 
         chrom is a pointer to the 0-based representation of the string
        
         returns a list containing:
         True/False - match or no match possible
         ref - true reference allele
         alt - true alternate allele
         message - either for the correction OR unknown - ref:TRUE REFERENCE 
         integer -  0, valid
                    1, ref/alt reverse
                    2, negative strand
                    3, ref/alt reverse + negative strand
                    -1, cannot resolve        
         matching order of operations
         1. correct ref on chrom @ start
         2. alt is correct on chrom @ start ( alt / ref reverse in raw data)
         3. reverse complement of ref is correct on chrom @ start (negative strand in raw data)
         4. reverse complment of alt is correct on chrom @ start (alt / reverse AND negative strand in raw data)
         5. cannot match to true chromosome
        """

        ref_end = start + len(ref)
        alt_end = start + len(alt)

        if chrom[start:ref_end] == ref: 
            return [True,ref,alt,'',0]
        elif chrom[start:alt_end] == alt:
            return [True,alt,ref,'ref_alt_rev',1]    
        elif ref != '-' and chrom[start:ref_end] == Seq(ref).reverse_complement().tostring():
            ref = Seq(ref).reverse_complement().tostring()
            alt = Seq(alt).reverse_complement().tostring()
            return [True,ref,alt,'neg_strand',2]    
        elif alt != '-' and chrom[start:alt_end] == Seq(alt).reverse_complement().tostring():
            ref = Seq(ref).reverse_complement().tostring()
            alt = Seq(alt).reverse_complement().tostring()
            return [True,alt,ref,'neg_strand_ref_alt_rev',3]    
        else:
            refStart = start - 5
            refStop = max(ref_end,alt_end) + 6
            print ref+'\t'+alt+'\t'+str(start)+'\t'+ str(ref_end)+'\t'+str(alt_end)+'\t'
            #return [False,ref,alt,'unknown - ref:'+chrom[start:max(ref_end,alt_end)],-1 ]
            return [False,ref,alt,'unknown - ref:'+chrom[start-5:start+6],-1 ]
    
    def prescan(self, inputfile):
        #global chrom_sequences
    
        """
        scan the first 100 "valid" lines of the file ref allele doesn't contain N or '-'
        determine if it is zero-indexed or one-indexed 
        """     
        file = None
        line = None
        try:
            file = open(inputfile)      
            line = file.readline()
        except:
            self.internalerr('unable to prescan '+inputfile,True)
            os._exit(0)

        if len(line) == 0:
            return 999
        
        limit     = 100     # limiting line count
        no_Offset = 0       # counter for zero-index references
        oneOffset = 0       # counter for one-indexed references
        count     = 0       # lines with valid references
        try:
            while line and count < limit:
                # break line into core parts, chrom / start/ end / ref
                parts = line.split('\t')
                line = file.readline()
                chrom = parts[1]
                start = int(parts[2])    
                ref   = parts[5] 
                alt   = parts[6] 
                   
                # IFF we find a new chromosome read it in and get it into a single continuous string 
                if self.chrom_sequences.has_key(chrom) != True:
                    #chromfile = chrom_dir + chrom + ".fa"
                    #inchrom = open(chromfile)
                    #seq = inchrom.read().replace(">","").replace(chrom,"").replace("\n","").upper()
                    #sequences[ chrom ] = seq
                    #inchrom.close()
                    self.load_chrom(chrom)
    
                current_sequence = self.chrom_sequences[ chrom ]
                 
                # skip over insertions or where the ref allele is obviously bad 
                if len(ref) < 1 or 'N' in ref or '-' in ref:
                    continue
    
                # roughly 1/4 of the time we'll end up with something like CGATTA
                # we need to be able to resolve this      
                # count all matches BUT 1 count will equal the limit
                
                noOffSetData = self.resolveVariant(ref,alt,start,current_sequence)
                oneOffSetData = self.resolveVariant(ref,alt,start-1,current_sequence)

                print (ref+'\t'+alt+'\t'+str(start)+'\t'+str(noOffSetData[4])+'\t'+str(oneOffSetData[4]))

                if noOffSetData[0]:
                    no_Offset += 1
                
                if oneOffSetData[0]:
                    oneOffset += 1 
    
                # track how many lines we've been able to compare to the reference genome
                count += 1
        except:
            traceback.print_exc()
            return -999
       
        print "1-Offset final: "+str(oneOffset) 
        print "0-offset final: "+str(no_Offset)
        
        if no_Offset > oneOffset:
            return 0
        elif oneOffset > no_Offset:
            return -1
        else:
            print 'unknown offset no_offset:'+str(no_Offset)+' vs oneOffset:'+str(oneOffset)
            return 1 

        """ 
        if oneOffset == limit:
            zov = -1
        elif no_Offset == limit:
            zov = 0
        elif oneOffset == count:
            # handles files with fewer than 100 lines
            zov = 1
        elif no_Offset == count:
            # handles files with fewer than 100 lines
            zov = 0
        else:
            # probably ran out of lines or hit ref sequence(s) that was crap
            zov = 999
        """
        return zov
        
    def load_chrom(self, chrom):
        #global genome_version,chrom_sequences
        
        print 'loading '+chrom
    
        chrom_dir = os.environ['HG_PATH'] + self.genome_version+"/chroms/" # source path to the chromosome fasta's
        chromfile = chrom_dir + chrom + ".fa"
        inchrom = open(chromfile)
        seq = inchrom.read().replace(">","").replace(chrom,"").replace("\n","").upper()
        self.chrom_sequences[ chrom ] = seq
        inchrom.close()
    
    
    def cleanup(self, inputfile,genome,outputfile=None):
        #global genome_version, chrom_sequences
        self.genome_version = genome
    
        self.chrom_sequences = {}  # dictionary of loaded chromosomes (complete string)
    
        zov = None   # (zeroOffsetValue either 0 or -1)
    
        zov = self.prescan(inputfile)
      
        errors = {}
       
        if zov < -1 or zov > 1:
            print 'failure - either an exception occurred or zero length file'
            self.internalerr('verify ZOV returned '+str(zov)+' check verify_variants.py for that number')
            errors['pass_validation'] = False
            errors['lineCount']       = 0    
            errors['totalCount']      = -1    
            return errors

             
        if zov not in [-1,0]:
            print str(zov)+' warning - pre-scan failed treating as zero-based coordinates'
            zov = 0

        print 'zov is : '+str(zov)
        
        cleanFile = inputfile+'.clean'        # good variants file name 
        excludeFile = inputfile+'.exclude'      # bad variants - will not annotate
        file_out = open(cleanFile,'w')          
        file_exclude = open(excludeFile,'w')   
        file_in = open(inputfile)
    
        line = file_in.readline()
    
        if line[0] == '#':
            # remove the header line
            line = file_in.readline()
    
    
        errors = {}
        errors['excludeFile'] = excludeFile
        errors['cleanFile'] = cleanFile
        errors['pass_validation'] = True # flag for return to signal to calling processes final status 
        errors['short']           = 0    # not enough data in the variant line 
        errors['long']            = 0    # too much data not going to parse it
        errors['chrom']           = 0    # bad chromosome 
        errors['ins']             = 0    # problems with a ins variant     
        errors['del']             = 0    # problems with a del      
        errors['snp']             = 0    # problems with a snp    
                                         # these two errors may overlap, both reverse & negative
        errors['ref_alt_rev']     = 0    # track how many rows had ref/alt columns in reverse order
        errors['strand']          = 0    # track how many rows were from the negative strand  
        errors['delins']          = 0    # problems with a delins     
        errors['type']            = 0    # unknown variant type     
        errors['base']            = 0    # bad base      
        errors['start_num']       = 0    # problems with a starting coordinate      
        errors['end_num']         = 0    # problems with an ending coordinate      
        errors['match']           = 0    # ref matchs alt allele
        #errors['error']           = ''   # empty unless we die early
        errors['lineCount']       = 0    # valid lines that became variants
        errors['totalCount']      = -1   # total lines that were read in 
                                         # needs to be -1 b/c of the initial check at 100 lines
                                         # error line count = totalCount - lineCount
        #errors['oneOffset']       = 0    # track how many times we have had to start-1/end-1
        #errors['no_Offset']       = 0    # track how many times we have had zero-based coords
        errors['refseq_N']        = 0    # how many times we hit an 'N' in the ref seq (snp only?)
        errors['zeroOffset']      = zov  # just for house keeping what was the value of the zero offset
        all_bases = ['A','C','G','T']    # only valid bases
    
        while line:
           
            errors['totalCount'] += 1  
            
            # as soon as there are 100 valid variants 
            """
            if errors['totalCount'] == 100:
                if errors['no_Offset'] == 100:
                    pass
                elif errors['oneOffset'] == 100:
                    zov = -1
                " " "
                else:
                    print 'errors['oneOffset']\t\t'+str(errors['oneOffset'])
                    print 'refSeq unknown (N)\t'+str(errors['refseq_N'])
                    print 'valid lines \t\t\t\t'+str(errors['lineCount'])
                    print 'total lines \t\t\t\t'+str(errors['totalCount'])
                    errors['error'] = 'unable to resolve coordinates between base 0 and 1'
                    for i in errors:
                        print i+'\t'+str(errors[i])
                    return errors
                " " "
                if errors['lineCount'] < errors['totalCount']:   
                    errors['error'] = 'unable to resolve coordinates between base 0 and 1'
                    errors['pass_validation'] = False
                    #print 'oneOffset\t\t\t'+str(errors['oneOffset'])
                    #print 'refSeq unknown (N)\t'+str(errors['refseq_N'])
                    #print 'valid lines \t\t\t'+str(errors['lineCount'])
                    #print 'total lines \t\t\t'+str(errors['totalCount'])
                    #for i in errors:
                    #    print i+'\t'+str(errors[i])
                    return errors
            """
            parts = line.split('\t')
           
            
            # pre-read the next line - variants that aren't acceptable log out
            # the error, write out the bad variant to the exlusion file
            # hit a continue and return to the top of the main while loop 
            line = file_in.readline()  
    
            # try and make a usable data line
            # if <6 elements exclude
            # if 6 elements try and verify column 1 (parts[0] ) as a chromosome
            # if 7 elements check both column 1 & 2 ( parts[0] and [1] ) for a chromosome
            if len(parts) < 8:
                # not enough data!
                if len(parts) < 6:
                    errors['short'] += 1
                    file_exclude.write('\t'.join(parts))
                    continue
                elif len(parts) == 6:
                    # presume we are missing the haplotype and notes column
                    # IFF the first column (parts0) is not chr then exclude
                    if parts[0][:3].lower() == 'chr':
                        # pad a haplotype
                        parts.insert(0,'3')
                    else:
                        errors['short'] += 1
                        file_exclude.write('\t'.join(parts))
                        continue
                    # add a "note" place holder    
                    parts.append('default_note1')    
    
                elif len(parts) == 7:    
                    # try and see if we're missing a haplotype
                    if parts[0][:3].lower() == 'chr':
                        # pad a haplotype
                        parts.insert(0,'3')
                    elif parts[1][:3].lower() == 'chr':
                        # add a "note" place holder    
                        parts.append('default_note2')
                    else:        
                        errors['short'] += 1
                        file_exclude.write('\t'.join(parts))
                        continue
            
            if len(parts) > 8:
                # too much data!
                errors['long'] += 1
                parts[7] = parts[7].strip()+'?0l\n'
                file_exclude.write('\t'.join(parts))
                continue
    
            # start actually checking the fields
            try:
                start = int(parts[2])
            except:
                # starting coordinate isn't a number
                errors['start_num'] += 1   
                parts[7] = parts[7].strip()+'?0s\n'
                file_exclude.write('\t'.join(parts))
                continue
            
            try:
                end = int(parts[3])
            except:
                # ending coordinate isn't a number 
                errors['end_num'] += 1   
                parts[7] = parts[7].strip()+'?0e\n'
                file_exclude.write('\t'.join(parts))
                continue
    
            if end < start:
                # the ending coordinate can't precede the starting
                errors['start_num'] += 1   
                parts[7] = parts[7].strip()+'?0s\n'
                file_exclude.write('\t'.join(parts))
                continue
    
    
            # check and make sure chromosome is legal
            parts[1] = parts[1].lower()
            if len( parts[1] ) != 4 and len( parts[1] ) != 5\
               and parts[1][:3] == 'chr':
                # must be either chr
                errors['chrom'] += 1
                parts[7] = parts[7].strip()+'?0c\n'
                file_exclude.write('\t'.join(parts))
                continue
            number = parts[1][3:]
            try:
                numb = int(number)
                if numb > 0 and numb < 23:
                    pass
                else:
                    # illegal chrom
                    errors['chrom'] += 1
                    parts[7] = parts[7].strip()+'?0c\n'
                    file_exclude.write('\t'.join(parts))
                    continue
            except:
                if number == 'X':
                    pass
                elif number == 'x':
                   parts[1] = 'chrX'
                elif number == 'Y':
                    pass
                elif number == 'y':
                   parts[1] = 'chrY'
                else:
                    #illegal chrom
                    errors['chrom'] += 1
                    parts[7] = parts[7].strip()+'?0c\n'
                    file_exclude.write('\t'.join(parts))
                    continue   
                                 
            # IFF we find a new chromosome read it in and get it into a single continuous string 
            if self.chrom_sequences.has_key(parts[1]) != True:
                self.load_chrom(parts[1])
                #chromfile = chrom_dir + parts[1] + ".fa"
                #inchrom = open(chromfile)
                #seq = inchrom.read().replace(">","").replace(parts[1],"").replace("\n","").upper()
                #chrom_sequences[ parts[1] ] = seq
                #inchrom.close()
    
            # get a pointer to the current chromosome sequence 
            current_sequence = self.chrom_sequences[ parts[1] ]
    
            # ensure ref / alt are upper case
            parts[5] = parts[5].upper()
            parts[6] = parts[6].upper()    
    
            # skip we don't annotate indeterminates
            if 'N' in parts[5]:
                errors['base'] += 1
                parts[7] = parts[7].strip()+'?0b\n'
                file_exclude.write('\t'.join(parts))
                continue
            if 'N' in parts[6]:
                errors['base'] += 1
                parts[7] = parts[7].strip()+'?0b\n'
                file_exclude.write('\t'.join(parts))
                continue
           
            # kick out variants that ref == alt allele
            if parts[5] == parts[6]:
                errors['match'] += 1
                parts[7] = parts[7].strip()+'?0m\n'
                file_exclude.write('\t'.join(parts))
                continue
            
            if len(current_sequence) < start:
                errors['start_num'] += 1   
                parts[7] = parts[7].strip()+'?0s\n'
                file_exclude.write('\t'.join(parts))
                continue
                
            elif len(current_sequence) < end:
                errors['end_num'] += 1   
                parts[7] = parts[7].strip()+'?0e\n'
                file_exclude.write('\t'.join(parts))
                continue
            
            
            # get the variant type right
            parts[4] = parts[4].lower()
    
            # fix known variant type error
            if parts[4] == 'sub' or parts[4] == 'indel' or\
               parts[4] == 'indels' or parts[4] == 'insdel':
                parts[4] = 'delins'
    
    
            if parts[4] == 'snp':
                if len(parts[5]) != 1 and len(parts[6]) != 1:
                    errors['snp'] += 1 
                    parts[7] = parts[7].strip()+'?1s\n'
                    file_exclude.write('\t'.join(parts))
                    continue
                elif parts[5] not in all_bases or\
                     parts[6] not in all_bases:
                    errors['snp'] += 1 
                    parts[7] = parts[7].strip()+'?1b\n'
                    file_exclude.write('\t'.join(parts))
                    continue
                elif current_sequence[start + zov] == 'N':
                    errors['refseq_N'] += 1
                else:
                    resolve = self.resolveVariant(parts[5],parts[6],start+zov,current_sequence) 
                
                    if resolve[0]:
                        # make the coordinates rights
                        parts[2] = str(start+zov)
                        parts[3] = str(start+zov+1)
                         
                        # we have resolved the variant correctly
                        parts[5] = resolve[1]
                        parts[6] = resolve[2]

                        if resolve[4] == 0:
                            pass
                        elif resolve[4] == 1:
                            errors['ref_alt_rev'] += 1                      # ref/alt columns reversed
                            parts[7] = parts[7].strip()+'?'+resolve[3]+'\n'
                        elif resolve[4] == 2:
                            errors['strand'] += 1                           # negative strand read
                            parts[7] = parts[7].strip()+'?'+resolve[3]+'\n'
                        elif resolve[4] == 3:
                            errors['strand'] += 1                           # negative strand read PLUS ref/alt columns reversed
                            errors['ref_alt_rev'] += 1
                            parts[7] = parts[7].strip()+'?'+resolve[3]+'\n'
                        else:
                            # for completeness / sanity
                            errors['snp'] += 1 
                            parts[7] = parts[7].strip()+'?1'+resolve[3]+'?1d:'+current_sequence[ (start + zov)-5:(start+zov)+6]+'\n'
                            file_exclude.write('\t'.join(parts))
                            continue
                    else:
                        errors['snp'] += 1 
                        parts[7] = parts[7].strip()+'?2'+resolve[3]+'?1d:'+current_sequence[ (start + zov)-5:(start+zov)+6]+'\n'
                        file_exclude.write('\t'.join(parts))
                        continue
                
                """
                elif current_sequence[start + zov] == parts[5]:
                    parts[2] = str(start + zov)                     # valid 
                    parts[3] = str(start+1 + zov)
                elif current_sequence[start + zov] == parts[6]:
                    errors['ref_alt_rev'] += 1                      # ref/alt columns reversed
                    parts[7] = parts[7].strip()+'?ref_alt_rev\n'
                    temp = parts[5]
                    parts[5] = parts[6]
                    parts[6] = temp
                    parts[2] = str(start + zov)
                    parts[3] = str(start+1 + zov)
                elif current_sequence[start + zov] == Seq(parts[5]).complement().tostring():
                    errors['strand'] += 1                           # negative strand read
                    parts[7] = parts[7].strip()+'?neg_strand\n'
                    parts[5] = Seq(parts[5]).complement().tostring()
                    parts[6] = Seq(parts[6]).complement().tostring()
                    parts[2] = str(start + zov)
                    parts[3] = str(start+1 + zov)
                elif current_sequence[start + zov] == Seq(parts[6]).complement().tostring():
                    errors['strand'] += 1                           # negative strand read PLUS ref/alt columns reversed
                    errors['ref_alt_rev'] += 1
                    parts[7] = parts[7].strip()+'?neg_strand_ref_alt_rev\n'
                    temp = Seq(parts[5]).complement().tostring()
                    parts[5] = Seq(parts[6]).complement().tostring()
                    parts[6] = temp
                    parts[2] = str(start + zov)
                    parts[3] = str(start+1 + zov)
                """
                #else:
                #    errors['snp'] += 1 
                #    parts[7] = parts[7].strip()+'?1d:'+current_sequence[ (start + zov)-5:(start+zov)+6]+'\n'
                #    file_exclude.write('\t'.join(parts))
                #    continue
            elif parts[4] == 'del':
                if len(parts[5]) != end-start and len(parts[6]) != end - start:
                    errors['del'] += 1
                    parts[7] = parts[7].strip()+'?2s\n'
                    file_exclude.write('\t'.join(parts))
                    continue
                elif current_sequence[ start+zov:end+zov] == parts[5]:
                    parts[2] = str(start + zov)
                    #parts[3] = str(end + zov)
                    parts[3] = str(start + zov + len(parts[5]))
                    parts[6] = '-'                                  # force the alt column to be correct
                elif current_sequence[ start+zov:end+zov] == parts[6]:
                    errors['ref_alt_rev'] += 1                      # ref/alt columns reversed
                    parts[7] = parts[7].strip()+'?ref_alt_rev\n'
                    parts[2] = str(start + zov)
                    parts[3] = str(start + zov + len(parts[6]))
                    #parts[3] = str(end + zov)
                    parts[5] = parts[6]
                    parts[6] = '-'                                  # force the alt column to be correct
                elif current_sequence[ start+zov:end+zov] == Seq(parts[5]).reverse_complement().tostring():
                    errors['strand'] += 1                           # negative strand read
                    parts[7] = parts[7].strip()+'?neg_strand\n'
                    parts[2] = str(start + zov)
                    parts[3] = str(start + zov + len(parts[5]))
                    #parts[3] = str(end + zov)
                    parts[5] = Seq(parts[5]).reverse_complement().tostring()
                    parts[6] = '-'                                  # force the alt column to be correct
                elif current_sequence[ start+zov:end+zov] == Seq(parts[6]).reverse_complement().tostring():
                    errors['strand'] += 1                           # negative strand read
                    errors['ref_alt_rev'] += 1                      # ref/alt columns reversed
                    parts[7] = parts[7].strip()+'?neg_strand_ref_alt_rev\n'
                    parts[2] = str(start + zov)
                    parts[3] = str(start + zov + len(parts[6]))
                    #parts[3] = str(end + zov)
                    parts[5] = Seq(parts[6]).reverse_complement().tostring()
                    parts[6] = '-'                                  # force the alt column to be correct
                else:
                    errors['del'] += 1
                    parts[7] = parts[7].strip()+'?2d\n'
                    file_exclude.write('\t'.join(parts))
                    continue
            elif parts[4] == 'ins':
                for base in parts[6]:
                    if base not in all_bases:
                        errors['ins'] += 1
                        parts[7] = parts[7].strip()+'?3b\n'
                        file_exclude.write('\t'.join(parts))
                        continue
                if len(parts[6]) > 0:
                    # single comparison b/c zov should be predetermined via pre-scan
                    parts[2] = str(start + zov)
                    #
                    #  FIXME write out a note to stdout with a notice about the forced change
                    #
                    parts[3] = parts[2]         # force the ending coordinate to be correct
                    parts[5] = '-'              # force the reference column to be correct
                    #errors['no_Offset'] += 1    # dumby counts but the numbers have to totalup 
                    #errors['oneOffset'] += 1    # dumby counts but the numbers have to totalup 
                else:
                    errors['ins'] += 1
                    parts[7] = parts[7].strip()+'?3d\n'
                    file_exclude.write('\t'.join(parts))
                    continue
            elif parts[4] == 'delins':
                for base in parts[6]:
                    if base not in all_bases:
                        errors['delins'] += 1
                        parts[7] = parts[7].strip()+'?4b\n'
                        file_exclude.write('\t'.join(parts))
                        continue
                
                if len(parts[5]) != end-start:
                    errors['delins'] += 1
                    parts[7] = parts[7].strip()+'?4s\n'
                    file_exclude.write('\t'.join(parts))
                    continue
                elif current_sequence[ start+zov:end+zov] == parts[5]:
                    # single comparison b/c zov should be predetermined via pre-scan
                    parts[2] = str(start + zov)
                    parts[3] = str(start + zov + len(parts[5]))
                    #parts[3] = str(end + zov)
                elif current_sequence[ start+zov:end+zov] == parts[6]:
                    errors['ref_alt_rev'] += 1                      # ref/alt columns reversed
                    parts[7] = parts[7].strip()+'?ref_alt_rev\n'
                    parts[2] = str(start + zov)
                    parts[3] = str(start + zov + len(parts[6]))
                    #parts[3] = str(end + zov)
                    temp = parts[5]
                    parts[5] = parts[6]
                    parts[6] = temp                                  # force the alt column to be correct
                elif current_sequence[ start+zov:end+zov] == Seq(parts[5]).reverse_complement().tostring():
                    errors['strand'] += 1                           # negative strand read
                    parts[7] = parts[7].strip()+'?neg_strand\n'
                    parts[2] = str(start + zov)
                    parts[3] = str(start + zov + len(parts[5]))
                    #parts[3] = str(end + zov)
                    parts[5] = Seq(parts[5]).reverse_complement().tostring()
                    parts[6] = Seq(parts[6]).reverse_complement().tostring() # force the alt column to be correct
                elif current_sequence[ start+zov:end+zov] == Seq(parts[6]).reverse_complement().tostring():
                    errors['strand'] += 1                           # negative strand read
                    errors['ref_alt_rev'] += 1                      # ref/alt columns reversed
                    parts[7] = parts[7].strip()+'?neg_strand_ref_alt_rev\n'
                    parts[2] = str(start + zov)
                    parts[3] = str(start + zov + len(parts[6]))
                    #parts[3] = str(end + zov)
                    temp = Seq(parts[5]).reverse_complement().tostring()
                    parts[5] = Seq(parts[6]).reverse_complement().tostring()
                    parts[6] = temp                                  # force the alt column to be correct
                else:
                    errors['delins'] += 1       
                    parts[7] = parts[7].strip()+'?4d\n'
                    file_exclude.write('\t'.join(parts))
                    continue
            else:
                errors['type'] += 1
                parts[7] = parts[7].strip()+'?5t\n'
                file_exclude.write('\t'.join(parts))
                continue
        
            errors['lineCount'] += 1    
            file_out.write('\t'.join(parts))
        
        file_out.close()
        file_exclude.close()
    
        errors['totalCount'] += 1 # need to add back for the final line
        if errors['lineCount'] != errors['totalCount']:
            errors['pass_validation'] = False
        else:
            os.system('rm '+excludeFile)
            if outputfile is not None:
                os.system(' cp '+cleanFile+' '+outputfile)
            else:
                stop = len(inputfile)
                if '.' in inputfile:
                    stop = inputfile.rindex('.')
                elif '_' in inputfile:
                    stop = inputfile.rindex('_')
                baseDir = inputfile[0:stop]
                pdata = baseDir+'.verify_data'
                os.system(' cp '+cleanFile+' '+pdata)
            os.system('rm '+cleanFile)

            #os.system('rm '+self.internalErr)

        return errors
    
    def processErrors(self,errDict):
        
        output = []
        
        output.append('Total Input Lines:               '+str(errDict['totalCount']))
        output.append('Valid Variants Found:            '+str(errDict['lineCount']))
        output.append('Negative Strand*:                '+str(errDict['strand']))
        output.append('Reference/Alternate Reverse*:    '+str(errDict['ref_alt_rev'])) 
        output.append('Invalid Variants Found:          '+str(errDict['totalCount'] - errDict['lineCount']))
        output.append('* These variants were detected, correct and are noted in the 8th data (notes) column')
        output.append('* Negative stranding and Ref/Alt Reverse are individual occurrance counts')
        output.append("* 'neg_strand' denotes a negative to positive strand correction for a variant")
        output.append("* 'ref_alt_rev' denotes a reversal of the reference and alternate columns for a particular line")
        
        output.append('')  # padding
        
        if errDict['totalCount'] - errDict['lineCount'] == 0:
            # early return to suppress the 
            return output
        
        output.append('')
        output.append(errDict['cleanFile']+'    contains all variants that were valid or corrected')
        output.append(errDict['excludeFile']+'  contains all variants that were invalid and uncorrectable')

        output.append('')  # padding
        output.append('****************************')
        output.append('Variant file errors include:')
        output.append('****************************')
        output.append('')  # padding
    
        if errDict['short'] >  0:
            output.append('There were errors with the number of fields in some lines, attempts to correct failed. ('+str(errDict['short'])+')')
            output.append('\tTry: Verify 8 columns on every line, tab separated, chromosome column (2) is in the format "chrX"')
            output.append('\tTry: If the 8th column contains no data, substitute a dash or a preiod as a place holder')
        if errDict['long'] >  0:
            output.append('There were errors with the number of fields in some lines. Only 8 columns are permitted. ('+str(errDict['long'])+') ')
            output.append('\tTry: Verify 8 columns on every line, tab separated, chromosome column (2) is in the format "chrX"')
        if errDict['chrom'] >  0:
            output.append('Invalid chromosome identifiers were found. ('+str(errDict['chrom'])+')')
            output.append("\tTry: Ensure the chromosomes are formatted as 'chrX' only regular chromosomes 1-22,X,Y are accepted.")
        if errDict['ins'] >  0:
            output.append('Malformed Insert variant found. ('+str(errDict['ins'])+')')
            output.append('\tTry: Ensure the alternate allele only contains valid bases [A,C,G,T]')
        if errDict['del'] >  0:
            output.append('Malformed Deletion variant found. ('+str(errDict['del'])+')')
            output.append('\tTry: Ensure the reference allele contains only valid bases and aligns to the appropriate reference genome')
            output.append("\tTry: Ensure the reference allele is correct with respect to the target genome at the specified coordinates")
        if errDict['snp'] >  0:
            output.append('Malformed SNP variant found. ('+str(errDict['snp'])+')')
            output.append('\tTry: Ensure the coordinates are correctly formatted, end = start + 1')
            output.append('\tTry: Ensure both the reference and alternate allele contain valid bases [A,C,G,T]')
            output.append("\tTry: Ensure the reference allele is correct with respect to the target genome at the specified coordinates")
        if errDict['delins'] >  0:
            output.append(' Malformed delins (indel) found. ('+str(errDict['delins'])+')')
            output.append('\tTry: Ensure the coordinates are correctly formatted end = start + len(ref allele)')
            output.append('\tTry: Ensure both the reference and alternate allele contain valid bases [A,C,G,T]')
            output.append("\tTry: Ensure the reference allele is correct with respect to the target genome at the specified coordinates")
        if errDict['type'] >  0:
            output.append('Unknown variant type found, could not correct. ('+str(errDict['type'])+')')
            output.append('\tTry: Only use valid variant types - snp / ins / del / delins ')
            output.append('\tTry: Remove invalid lines from input file')
        if errDict['base'] >  0:
            output.append('Invalid base found in a reference or alternate allele. ('+str(errDict['base'])+')')
            output.append('\tTry: Ensure both the reference and alternate allele contain valid bases [A,C,G,T]')
            output.append('\tTry: Verify source data')
        if errDict['start_num'] >  0:
            output.append('Illegal starting coordinate found. ('+str(errDict['start_num'])+')')
            output.append('\tTry: Verify all starting coordinates are integers, also the start <= end')
            output.append('\tTry: Verify source data')
            output.append('\tTry: Remove invalid lines from input file')
        if errDict['end_num'] >  0:
            output.append('Illegal ending coordinate found. ('+str(errDict['end_num'])+')')
            output.append('\tTry: Verify all starting coordinates are integers, also the start <= end')
            output.append('\tTry: Verify source data')
            output.append('\tTry: Remove invalid lines from input file')
        if errDict['match'] >  0:
            output.append('Matching reference and alternative allele. ('+str(errDict['match'])+')')
            output.append('\tTry: Verify source data')
            output.append('\tTry: Remove invalid lines from input file')
        #if errDict['error'] != '':
        #    output.append('')
        #    output.append('\tTry:')
    
        return output

    def __init__(self,stdout=None,stderr=None,fail=None):
        # declare parse type symbols
        #for f in os.listdir('.'):
        #   print f
    
        sga_base.__init__(self,stdOUT=stdout,stdERR=stderr,failureFile=fail)
        
        #for f in os.listdir('.'):
        #    self.usermsg(f+'\n')


if __name__ == "__main__":
    
    #if len(sys.argv) < 3:
    #    print '\n1 usage: verify_variants.py [path to variants] [hg19] [opt: stderr_file]\n'
    #    os._exit(-1)
    #if sys.argv[2].lower() != 'hg19':
    #    print '\n2 usage: verify_variants.py [path to variants] [hg19] [opt: stderr_file]\n'
    #    os._exit(-2)

    parser = argparse.ArgumentParser(description='Arguments for Pipeline Driver')
    parser.add_argument('-f','--input_file',required=True, help='required input file for processing')
    parser.add_argument('-o','--output_file',required=False, help='required input file for processing')
    parser.add_argument('-s','--user_out',default=".verify_out",required=False, help='optional file for stdout messages to the user')
    parser.add_argument('-e','--user_err',default=".verify_err",required=False, help='optional file for stderr messages to the user')
    parser.add_argument('-i','--internal',default=".verify_internalERRor",required=False, help='optional file for internal fatal messages')
    parser.add_argument('-g','--genome',required=True, help='what genome to use as reference')
    args = vars(parser.parse_args())

    inx = len(args['input_file'])
    if '.' in args['input_file']:
        inx = args['input_file'].rindex('.')
    elif '_' in args['input_file']:
        inx = args['input_file'].rindex('_')

    handle = args['input_file'][0:inx]

    valid = validator(handle+args['user_out'],handle+args['user_err'],handle+args['internal'])
    output = valid.cleanup(args['input_file'],'hg19',outputfile=args['output_file'])
    
    resultArray = valid.processErrors(output)    

    if not output['pass_validation']:
        valid.usererr('\n'.join(resultArray))
        valid.usererr('\n')
    else:
        valid.usermsg('\n'.join(resultArray))
        valid.usermsg('\n')


