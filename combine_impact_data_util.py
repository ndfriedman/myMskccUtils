#written by Noah Friedman 
#a script intended to be a pipleine for impact data aglommeration and other stuff

import sys
import argparse
import os
import pandas as pd
import numpy as np

import data_compacting_and_cleaning_util


def combine_case_based_annotations(args, caseDfIdCol):
	caseBasedDfsToAmalgamate  = [] # a list of case based dfs to amalgamate (clinical dfs etc)
	if args.clinicalInfoTable != None: 
		clinicalDf = pd.read_table(args.clinicalInfoTable)
		clinicalDf = clinicalDf.rename(columns={'Patient_Identifier': caseDfIdCol})
		caseBasedDfsToAmalgamate.append(clinicalDf)
	if args.survivalData != None: caseBasedDfsToAmalgamate.append(pd.read_table(args.survivalData))
	if args.drugData != None: caseBasedDfsToAmalgamate.append(pd.read_table(args.drugData))
	#ALERT do renames if needed for casedf id col
	mergedCaseBasedDf = data_compacting_and_cleaning_util.amalgamate_dfs_same_col(caseBasedDfsToAmalgamate, caseDfIdCol, 'combinedCaseBasedAnnotations.tsv', mode='return')
	return mergedCaseBasedDf

def combine_sample_based_annotations(args, sampleDfIdCol):
	sampleBasedDfsToAmalgamate = [] # a list of sample based dfs to amalgamate (signatures etc)
	if args.signaturesTable != None: sampleBasedDfsToAmalgamate.append(pd.read_table(args.signaturesTable))
	mergedSampleBasedDf = data_compacting_and_cleaning_util.amalgamate_dfs_same_col(sampleBasedDfsToAmalgamate, sampleDfIdCol, 'combinedSampleBasedAnnotations.tsv', mode='return')
	return mergedSampleBasedDf

def combine_variant_based_dfs(args, somaticMutations, germlineMutations):
	variantBasedDfsToAmalgamate = [] # a list of dfs with information on a per variant basis (ie somatic muts, CNAs etc)
	print 'variant only mode not implemented'

def combine_gene_based_dfs(args, somaticMutations, germlineMutations, makeGenesCols=False):

	variantDf = pd.read_table(args.impactSomaticMutations)
	print variantDf
	sys.exit()

	geneBasedDfsToAmalgamage = [] #a list of dfs with information on a per gene basis
	dfCNA = data_compacting_and_cleaning_util.parse_data_format_cases_on_line_one(args.cnaData)
	print dfCNA.shape
	sys.exit()
	geneBasedDfsToAmalgamage.append(args.cnaData)
	if makeGenesCols:
		print 'TODO implement this'
	print 'gene only mode not implemented'

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--amalgamatedFilePath', help='a user specified path to an amalgamated file.  if none, as default, create an amagamated file', default=None)
	parser.add_argument('--signaturesTable', help='Table of signature data', default='/ifs/work/taylorlab/friedman/signatureInvestigation/inputFiles/impactSignaturesAlexGenerated.txt')
	parser.add_argument('--facetsData', help='path to maf annotated with facets data', default='/ifs/work/taylorlab/friedman/signatureInvestigation/inputFiles/data_mutations_extended.mafAnno_fast.txt')	
	
	parser.add_argument('--clinicalInfoTable', help='Table of clinical info table', default='/ifs/work/taylorlab/friedman/signatureInvestigation/inputFiles/clinicalPatientDataReheadered.txt')
	parser.add_argument('--survivalData', help='path to survival data', default='/ifs/work/taylorlab/friedman/signatureInvestigation/inputFiles/survivalData.txt')
	parser.add_argument('--drugData', help='path to data about drug treatment data', default='/ifs/work/taylorlab/friedman/clinicalData/msk-impact/msk-impact/data_timeline.txt')


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

		print args.clinicalInfoTable
		mergedCaseBasedDf = combine_case_based_annotations(args, 'PATIENT_ID')
		mergedSampleBasedDf = combine_sample_based_annotations(args, 'Tumor_Sample_Barcode')
		print mergedCaseBasedDf
		print mergedSampleBasedDf

		somaticMutations = None
		germlineMutations = None
		if args.impactSomaticMutations: somaticMutations = pd.read_table(args.impactSomaticMutations) 
		if args.impactGermlineMutations: germlineMutations = pd.read_table(args.impactGermlineMutations)
		
		if args.geneVsVariantLevelAnnotationMode == 'geneOnly':
			combine_gene_based_dfs(args, somaticMutations, germlineMutations)

		elif args.geneVsVariantLevelAnnotationMode == 'variantOnly':
			combine_variant_based_dfs(args)

		elif args.geneVsVariantLevelAnnotationMode == 'both':
			print 'both variant and gene mode not implemented'
		
		else:
			print 'invalid geneVsVariant annotation option specifed.  valid options are: both, variantOnly and geneOnly'



if __name__ == '__main__':
    main()















