import xml.etree.ElementTree as ET
import os, sys, datetime
outFilename = str(sys.argv[1])
outFile = open(outFilename, 'w')

print "herro"

tree = ET.parse('/gpfs/home/gerikson/CNV_pipeline/Clinvar_cnv_02:2014/XML_format/ClinVarFullRelease_00-latest.xml')
#tree = ET.parse('/gpfs/home/gerikson/CNV_pipeline/Clinvar_cnv_02:2014/XML_format/ClinVar10k.xml')
#tree = ET.parse('/gpfs/home/gerikson/CNV_pipeline/Clinvar_cnv_02:2014/XML_format/cnv.xml')
root = tree.getroot()

#Extract Accesion
for child in root.findall('./ClinVarSet/ReferenceClinVarAssertion'):
	#extract varType
	vartype_temp = child.find('MeasureSet/Measure')
	vartype = vartype_temp.get('Type')

	#extract clinical significance
	clin = child.find('ClinicalSignificance/Description').text

	#extract review status
	rev = child.find('ClinicalSignificance/ReviewStatus').text

	#Extract location
	chrom = ""
	start = ""
	end = ""

	#Extract only the copy number gains and losses
	if vartype == 'copy number gain' or vartype == 'copy number loss':
		location = child.find('MeasureSet/Measure/SequenceLocation')
		if location is not None:
			chrom = location.get('Chr')
			start = location.get('innerStart')
			end = location.get('innerStop')
			Acc_temp = child.find('ClinVarAccession')
			Acc = Acc_temp.get('Acc')
			#Extract Phenotype
			Phenotype = child.find('TraitSet/Trait/Name/ElementValue').text
			temp_string = 'chr'+str(chrom) + "\t" + str(start) + "\t" + str(end) + "\t" + str(vartype) + "\t" + clin + "\t" + rev + "\t" + Phenotype + "\t" + str(Acc) + "\n"
			outFile.write(temp_string.encode("UTF-8"))
			print temp_string
			#print 'chr'+str(chrom) + "\t" + str(start) + "\t" + str(end) + "\t" + str(vartype) + "\t" + clin + "\t" + rev + "\t" + Phenotype + "\t" + str(Acc)


outFile.close()