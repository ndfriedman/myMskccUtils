#written by noah friedman
#takes an input maf and returns a mutation spectra file 
#the mutation spectra file is a 97 col tsv: Tumor_Sample_Barcode plus 96 quadnuc (ACTG etc)
#usage python calculate_nucleotide_mutation_spectra_from_maf.py inputMaf path/to/output/filename.tsv 
#versioning: python 2.7.1, pandas 0.19.2

import sys
import pandas as pd

#UTILITY function to merge dicitonaries
def Merge(dict1, dict2):  
	z = dict2.copy()
	z.update(dict1)
	return z

#df['quadNuc'] = df.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
def create_reference_four_nuc(refTri, refAllele, altAllele, variantType):
	#properly invert when needed
	def invert_allele_and_ref_tri(altAllele):
		nucleotideDict = {'A': 'T', 'G': 'C', 'C':'G', 'T':'A'}
		return nucleotideDict[altAllele]

	if variantType != 'SNP': return None  #there is no reference trinuc for non snps
	if not isinstance(refTri, basestring): return None #if the ref tri is not a string we better return none
	if len(refTri) < 3: return None # if the ref tri is less than length 3 (at the end of an exon), we cant do anything
	refTri = str(refTri) #just make sure the ref tri is a string here to avoid funny business
	refAlleleFromReftri = refTri[1]
	alt = altAllele
	if refAlleleFromReftri != refAllele:
		alt = invert_allele_and_ref_tri(altAllele)
	quadNuc = refTri[:2] + alt + refTri[2]
	return quadNuc

#main function that returns a dataframe with mutation spectra
def count_quad_nuc_sprectra_from_maf(maf):

	listOfDataframes = []

	if ['Ref_Tri'] not in maf.columns.values:
		print 'error, the maf we use needs to have a column called ref tri'
		return
	allBases = ['A', 'C', 'G', 'T']
	changes = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG'] #format: 'CA' means a change from C>A
	allQuadNucs = [firstBase + change + lastBase for firstBase in allBases for change in changes for lastBase in allBases] #enumerate all 96 quadnucs for signatures
	maf['quadNuc'] = maf.apply(lambda row: create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
	allCases = set(maf['Tumor_Sample_Barcode'])
	for quadNuc in allQuadNucs:
		qMaf = maf[maf['quadNuc'] == quadNuc]

		#Convluted (but hopefully somewhat more efficient) code to count the number of times a trinucleotide occurs in each case
		valueCountsDict = dict(qMaf['Tumor_Sample_Barcode'].value_counts()) #cases where the trincleotide occurs
		zeroBarcodes = allCases - set(valueCountsDict.keys())
		casesWhereQuadNucDoesNotOccur = {key:value for (key, value) in [(barcode, 0) for barcode in zeroBarcodes]}#create another dictionary to mark other cases as 0
		allCasesDictionary = Merge(valueCountsDict, casesWhereQuadNucDoesNotOccur)
		df = pd.DataFrame(allCasesDictionary.items())
		df = df.rename(columns={0:'Tumor_Sample_Barcode', 1:quadNuc}) #fix the column names of this df we have constructed
		listOfDataframes.append(df) #store all these columns so we can merge them later

	df = listOfDataframes[0]
	for df_ in listOfDataframes[1:]:
		df = df.merge(df_, on='Tumor_Sample_Barcode')
	return df

mafPath = sys.argv[1]
writePath = sys.argv[2]
maf = pd.read_table(mafPath)
spectraDf = count_quad_nuc_sprectra_from_maf(maf)
print 'writing file to ', writePath
spectraDf.to_csv(writePath, index=False, sep='\t')




