import xml.etree.ElementTree as ET
import os, sys, datetime
outFilename = str(sys.argv[1])
outFile = open(outFilename, 'w')

print "herro"
#tree = ET.parse('ClinVar10k.xml')
#tree = ET.parse('/gpfs/home/gerikson/CNV_pipeline/Clinvar_cnv_02:2014/XML_format/ClinVarFullRelease_00-latest.xml')
tree = ET.parse('/gpfs/home/gerikson/CNV_pipeline/Clinvar_cnv_02:2014/XML_format/ClinVar10k.xml')
root = tree.getroot()

#i = 0
'''
for items in root:
	if items.attrib['MeasureSet'] == "copy number loss": 
		print "copy number loss \n"
'''

'''
for child in root.iter('Measure'):
	print child.attrib

for child in root.iter('SequenceLocation'):
	print child.attrib

for child in root.iter('ClinVarAccession'):
	print child.attrib

'''
'''
for child in root.findall('./ClinVarSet/ReferenceClinVarAssertion'):
	print child.attrib
'''

for child in root.findall('./ClinVarSet/ReferenceClinVarAssertion/MeasureSet/Measure/SequenceLocation'):
	chrom = child.get('Chr')
	start = child.get('start')
	end = child.get('stop')
	#print chrom
	#print start
	#print end
	outFile.write('chr'+str(chrom) + "\t" + str(start) + "\t" + str(end) + "\n")
	#print child.attrib

outFile.close()