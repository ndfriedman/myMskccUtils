#written by Noah Friedman 
#a script intended to be a pipleine for impact data aglommeration and other stuff

import sys
import argparse
import os
import pandas as pd
import numpy as np

import data_compacting_and_cleaning_util

#TODO: make these functions parameter based instead of arguemnt based

#function to combine case based annotations (patient drugs etc)
def combine_case_based_annotations(clinicalInfoTablePath=None, survivalDataPath=None, drugDataPath=None,
 caseDfIdCol='PATIENT_ID', writeCombinedDf=False, outputFilename='combinedCaseBasedAnnotations.tsv', outputDirPath='/ifs/work/taylorlab/friedman/myAdjustedDataFiles'):

	caseBasedDfsToAmalgamate  = [] # a list of case based dfs to amalgamate (clinical dfs etc)
	if clinicalInfoTablePath != None: 
		clinicalDf = pd.read_table(clinicalInfoTablePath)
		clinicalDf = clinicalDf.rename(columns={'Patient_Identifier': caseDfIdCol})
		caseBasedDfsToAmalgamate.append(clinicalDf)
	if survivalDataPath != None: caseBasedDfsToAmalgamate.append(pd.read_table(survivalDataPath))
	if drugDataPath != None: caseBasedDfsToAmalgamate.append(pd.read_table(drugDataPath))
	#ALERT do renames if needed for casedf id col
	mergedCaseBasedDf = data_compacting_and_cleaning_util.amalgamate_dfs_same_col(caseBasedDfsToAmalgamate, caseDfIdCol, 'combinedCaseBasedAnnotations.tsv', mode='return')
	
	if writeCombinedDf:
		writePath = os.path.join(outputDirPath, outputFilename)
		print 'writing data to ', writePath
		dfFinal.to_csv(writePath, index=False, sep='\t')

	return mergedCaseBasedDf

#util function to combine sample specific annotations
def combine_sample_based_annotations(signatureTablePath=None, sampleInfoPath=None, cnaSummaryInfoPath=None,
	sampleDfIdCol='Tumor_Sample_Barcode'):

	sampleBasedDfsToAmalgamate = [] # a list of sample based dfs to amalgamate (signatures etc)
	if signatureTablePath != None: sampleBasedDfsToAmalgamate.append(pd.read_table(signatureTablePath))
	if sampleInfoPath != None:
		sampleDf = pd.read_table(sampleInfoPath)
		sampleDf = sampleDf.rename(columns ={'#Sample Identifier': 'Tumor_Sample_Barcode'})
		sampleBasedDfsToAmalgamate.append(sampleDf)
	if cnaSummaryInfoPath != None:
		sampleBasedDfsToAmalgamate.append(pd.read_table(cnaSummaryInfoPath))
	mergedSampleBasedDf = data_compacting_and_cleaning_util.amalgamate_dfs_same_col(sampleBasedDfsToAmalgamate, sampleDfIdCol, 'combinedSampleBasedAnnotations.tsv', mode='return')
	return mergedSampleBasedDf

def get_specific_cancer_type_case_ids(cancerTypes, pathToCancerTypeLists='/ifs/work/taylorlab/friedman/clinicalData/msk-impact/msk-impact/case_lists'):
	idSet = set()
	for cancerType in cancerTypes:
		path = os.path.join(pathToCancerTypeLists, cancerType)
		with open(path) as f:
			lines = f.readlines()
			ids = lines[4]
			caseIds = ids.strip('\n').split('\t')[1:]
			idSet = idSet | set(caseIds)
		f.close()
	return idSet

def combine_gene_based_dfs(
	somaticMutationsPath=None, signaturesDFPath=None, cnaDfPath=None, makeGenesCols=False,
	writeCombinedDf=False, outputFilename='geneBasedAnnotations.tsv', outputDirPath='/ifs/work/taylorlab/friedman/myAdjustedDataFiles'):

	#util function to rename gene columns with a suffix so there isnt a conflict for gene amplifications/variants etc
	def rename_cols_with_suffix(df, colSuffix):
		cols = list(df.columns.values)
		cols.remove('Tumor_Sample_Barcode')
		renameDict = {key: value for (key, value) in [(string, string + colSuffix) for string in cols]}
		df = df.rename(columns=renameDict)
		return df

	geneBasedDfsToAmalgamage = [] #a list of dfs with information on a per gene basis
	
	if cnaDfPath != None:
		dfCNA = data_compacting_and_cleaning_util.parse_data_format_cases_on_line_one(cnaDfPath)
		dfCNA = rename_cols_with_suffix(dfCNA, '_CNA')
		geneBasedDfsToAmalgamage.append(dfCNA)

	if somaticMutationsPath != None: #TODO qualify this information based on how we are going to treat mutations
		geneMatrix = data_compacting_and_cleaning_util.parse_mut_data_into_gene_matrix(somaticMutationsPath, write=False)
		geneMatrix = rename_cols_with_suffix(geneMatrix, '_MUT')
		geneBasedDfsToAmalgamage.append(geneMatrix)
	
	if signaturesDFPath != None:
		dfSignatures = pd.read_table(signaturesDFPath)
		geneBasedDfsToAmalgamage.append(dfSignatures)

	dfFinal = reduce(lambda left,right: pd.merge(left,right,on='Tumor_Sample_Barcode'), geneBasedDfsToAmalgamage)
	
	if writeCombinedDf:
		writePath = os.path.join(outputDirPath, outputFilename)
		print 'writing data to ', writePath
		dfFinal.to_csv(writePath, index=False, sep='\t')

	return dfFinal

# a utility function for combining multiple dfs with information about variants
#creates a join col based on chromosome, pos alt
def combine_variant_based_dfs(mafAnnoDf=None, triuncDf=None,
 writeCombinedDf=False, outputFilename='mergedAnnotatedMaf.maf', outputDirPath='/ifs/work/taylorlab/friedman/myAdjustedDataFiles'):

	#create a unique identifier for joining the dfs
	def create_chr_pos_join_key(row):
		return str(row['Chromosome']) + ':' + str(row['Start_Position']) + '->'  + str(row['End_Position']) + ';' + str(row['Tumor_Seq_Allele2'])

	dfsToAmalgamate = []
	if mafAnnoDf != None:
		mafAnnoData = pd.read_table(mafAnnoDf)
		mafAnnoData['joinCol'] = mafAnnoData.apply(lambda row: create_chr_pos_join_key(row), axis=1)
		dfsToAmalgamate.append(mafAnnoData)

	if triuncDf != None:
		triuncData = pd.read_table(triuncDf)
		triuncData = triuncData[['Chromosome', 'Start_Position', 'End_Position', 'Tumor_Seq_Allele2', 'Ref_Tri']]
		triuncData['joinCol'] = triuncData.apply(lambda row: create_chr_pos_join_key(row), axis=1)
		dfsToAmalgamate.append(triuncData)

	mergedDf = data_compacting_and_cleaning_util.amalgamate_dfs_same_col(dfsToAmalgamate, 'joinCol')

	if writeCombinedDf:
		writePath = os.path.join(outputDirPath, outputFilename)
		print 'writing data to ', writePath
		mergedDf.to_csv(writePath, index=False, sep='\t')

	return mergedDf

#a function for merging dfs of different sizes (ie a df of tumor sample barcodes and a df of clinical ids)
def merge_dfs_of_different_sizes(tsBarcodeDf, clinicalIdDf,
	writeCombinedDf=False, outputFilename='mergedAnnotations.tsv', outputDirPath='/ifs/work/taylorlab/friedman/myAdjustedDataFiles'):
	
	tsBarcodeDf['PATIENT_ID'] = tsBarcodeDf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
	combinedDf = clinicalIdDf.merge(tsBarcodeDf)

	if writeCombinedDf:
		writePath = os.path.join(outputDirPath, outputFilename)
		print 'writing data to ', writePath
		combinedDf.to_csv(writePath, index=False, sep='\t')

	return combinedDf

#simple merge and write wrapper for merges that require no manipulation of constituent dfs
def simple_df_merge(df1, df2, mergeOn=None,
	writeCombinedDf=False, outputFilename='mergedAnnotations.tsv', outputDirPath='/ifs/work/taylorlab/friedman/myAdjustedDataFiles'):
	mergedDf = df1.merge(df2, on=mergeOn)
	if writeCombinedDf:
		writePath = os.path.join(outputDirPath, outputFilename)
		print 'writing data to ', writePath
		mergedDf.to_csv(writePath, index=False, sep='\t')
	return mergedDf

#returns a list of all the genes in the impact panel
def enumerate_impact_panel_genes(panelFile='/ifs/work/taylorlab/friedman/msk-impact/msk-impact/gene_panels/impact468_gene_panel.txt'):
	with open(panelFile) as f:
		lines = f.readlines()
		panelLine = lines[3]
		return panelLine.strip('\n').split('\t')[1:]

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--amalgamatedFilePath', help='a user specified path to an amalgamated file.  if none, as default, create an amagamated file', default=None)
	parser.add_argument('--signaturesTable', help='Table of signature data', default='/ifs/work/taylorlab/friedman/signatureInvestigation/inputFiles/impactSignaturesAlexGenerated.txt')
	parser.add_argument('--facetsData', help='path to maf annotated with facets data', default='/ifs/work/taylorlab/friedman/signatureInvestigation/inputFiles/data_mutations_extended.mafAnno_fast.txt')	
	
	parser.add_argument('--clinicalInfoTable', help='Table of clinical info table', default='/ifs/work/taylorlab/friedman/signatureInvestigation/inputFiles/clinicalPatientDataReheadered.txt')
	parser.add_argument('--survivalData', help='path to survival data', default='/ifs/work/taylorlab/friedman/signatureInvestigation/inputFiles/survivalData.txt')
	parser.add_argument('--drugData', help='path to data about drug treatment data', default='/ifs/work/taylorlab/friedman/clinicalData/msk-impact/msk-impact/data_timeline.txt')


	#/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impactVarsWithTriunc.maf #VERSION OF MAF WITH REF TRI INFO
	parser.add_argument('--impactSomaticMutations', help='the path to the data containing impact somatic mutations', default='/ifs/work/taylorlab/friedman/clinicalData/msk-impact/msk-impact/data_mutations_extended.txt')
	parser.add_argument('--impactGermlineMutations', help='the path to the data containing impact germline mutations', default=None)

	parser.add_argument('--clonalityData', help='path to clonality info about the sample', default='/ifs/work/taylorlab/friedman/signatureInvestigation/clonalityInfo.tsv')
	parser.add_argument('--cnaData', help='', default='/ifs/work/taylorlab/friedman/clinicalData/msk-impact/msk-impact/data_CNA.txt')


	#arg parser modes
	parser.add_argument('--geneVsVariantLevelAnnotationMode', help='a mode dictating whether the analysis of certain features such as facets or somatic mutations should be performed on the gene or variant level. There are three valid options: geneOnly, variantOnly, or both', default='geneOnly')

	args = parser.parse_args()

	amalgamatedDf  = None
	if args.amalgamatedFilePath != None:
		amalgamatedDf = pd.read_table(args.amalgamatedFilePath)
	else: 


		"""print args.clinicalInfoTable
		mergedCaseBasedDf = combine_case_based_annotations(args, 'PATIENT_ID')
		mergedSampleBasedDf = combine_sample_based_annotations(args, 'Tumor_Sample_Barcode')
		print mergedCaseBasedDf
		print mergedSampleBasedDf"""

		if args.geneVsVariantLevelAnnotationMode == 'geneOnly':
			combine_gene_based_dfs(args.impactSomaticMutations)

		elif args.geneVsVariantLevelAnnotationMode == 'variantOnly':
			combine_variant_based_dfs(args)

		elif args.geneVsVariantLevelAnnotationMode == 'both':
			print 'both variant and gene mode not implemented'
		
		else:
			print 'invalid geneVsVariant annotation option specifed.  valid options are: both, variantOnly and geneOnly'



if __name__ == '__main__':
    main()















