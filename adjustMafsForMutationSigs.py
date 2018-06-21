#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

#
def map_barcode_to_cancer_type(df):
	d = {}
	for index, row in df.iterrows():
		d[row['Tumor_Sample_Barcode']] = row['Cancer Type']
	return d

def remove_hypermutators():
	return 0

def remove_multiple_mutations(df):
	df['idCol'] =  df.apply(lambda row: '.'.join([str(row['Chromosome']), str(row['Start_Position']), str(row['End_Position'])]), axis=1)
	df = df.drop_duplicates(subset=['idCol'])
	return df

###############################################################

def reheader_with_cancer_types_as_id(df, d, outputDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles', filename='mutationsWithCancerTypeAsId.maf'):
	localDf = df.copy()
	localDf = localDf[localDf['Tumor_Sample_Barcode'].isin(set(d.keys()))]
	localDf['Tumor_Sample_Barcode'] = localDf['Tumor_Sample_Barcode'].apply(lambda x: d[x])
	filepath = os.path.join(outputDir, filename)
	localDf.to_csv(filepath, sep = '\t', index=False)

def reheader_with_hugo_symbol_as_id(df, outputDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles', filename='mutationsWithHugoSymbolAsId.maf'):
	localDf = df.copy()
	localDf['idCol'] =  localDf.apply(lambda row: '.'.join([str(row['Chromosome']), str(row['Start_Position']), str(row['End_Position'])]), axis=1)
	localDf = localDf.drop_duplicates(subset=['idCol'])
	localDf['Tumor_Sample_Barcode'] = localDf.apply(lambda row: row['Hugo_Symbol'], axis=1)
	filepath = os.path.join(outputDir, filename)
	localDf.to_csv(filepath, sep = '\t', index=False)

def reheader_with_hugo_symbol_and_cancer_type_as_id(df, d, outputDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles', filename='mutationsWithHugoSymbolAsId.maf'):
	localDf = df.copy()
	localDf = localDf[localDf['Tumor_Sample_Barcode'].isin(set(d.keys()))]
	localDf['Tumor_Sample_Barcode'] = localDf.apply(lambda row: str(d[row['Tumor_Sample_Barcode']]) + '_' + str(row['Hugo_Symbol']), axis=1)
	filepath = os.path.join(outputDir, filename)
	localDf.to_csv(filepath, sep = '\t', index=False)

def reheader_with_hugo_symbol_pole_specific(df, outputDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles', filename='mutationsWithHugoSymbolAsIdInPoleCases.maf'):
	def get_pole_cases(mutSigsFilepath = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mutationSigFiles/may16unfiltered30sigs.txt'):
		mutSigsDf = pd.read_table(mutSigsFilepath)
		hypermutatedPoles = mutSigsDf[(mutSigsDf['Signature.10'] > .15) & (mutSigsDf['Number of Mutations'] > 30)] #define hypermutatedPoles as cases with 15% pole signature and over 30 mutations
		return set(hypermutatedPoles['Sample Name'])

	poleCaseIds = get_pole_cases()
	localDf = df.copy()
	localDf = localDf[localDf['Tumor_Sample_Barcode'].isin(poleCaseIds)]
	localDf['Tumor_Sample_Barcode'] = localDf.apply(lambda row: row['Hugo_Symbol'], axis=1)
	filepath = os.path.join(outputDir, filename)
	print 'writing to ', filepath
	localDf.to_csv(filepath, sep = '\t', index=False)

def reheader_with_hugo_symbol_mmr_specific(df, outputDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles', filename='mutationsWithHugoSymbolAsIdInMMRCases.maf'):
	def get_mmr_cases(mutSigsFilepath = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mutationSigFiles/may16unfiltered30sigs.txt'):
		mutSigsDf = pd.read_table(mutSigsFilepath)
		hypermutatedMMRs = mutSigsDf[
		(mutSigsDf['Signature.1'] + mutSigsDf['Signature.6'] + mutSigsDf['Signature.15'] + mutSigsDf['Signature.20'] + mutSigsDf['Signature.26'] > .25)
		&
		(mutSigsDf['Number of Mutations'] > 50)] #define hypermutatedMMRs as cases with 15% pole signature and over 30 mutations
		return set(hypermutatedMMRs['Sample Name'])

	mmrCaseIds = get_mmr_cases()
	localDf = df.copy()
	localDf = localDf[localDf['Tumor_Sample_Barcode'].isin(mmrCaseIds)]
	localDf['Tumor_Sample_Barcode'] = localDf.apply(lambda row: row['Hugo_Symbol'], axis=1)
	filepath = os.path.join(outputDir, filename)
	print 'writing to ', filepath
	localDf.to_csv(filepath, sep = '\t', index=False)


def reheader_with_hugo_symbol_tmz_specific(df, outputDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles', filename='mutationsWithHugoSymbolAsIdInTMZCases.maf'):
	def get_tmz_cases(mutSigsFilepath = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mutationSigFiles/may16unfiltered30sigs.txt'):
		mutSigsDf = pd.read_table(mutSigsFilepath)
		hypermutatedTMZs = mutSigsDf[
		(mutSigsDf['Signature.11'] > .15) & (mutSigsDf['Number of Mutations'] > 50)] #define hypermutatedMMRs as cases with 15% pole signature and over 30 mutations
		return set(hypermutatedTMZs['Sample Name'])

	tmzCaseIds = get_tmz_cases()
	localDf = df.copy()
	localDf = localDf[localDf['Tumor_Sample_Barcode'].isin(tmzCaseIds)]
	localDf['Tumor_Sample_Barcode'] = localDf.apply(lambda row: row['Hugo_Symbol'], axis=1)
	filepath = os.path.join(outputDir, filename)
	print 'writing to ', filepath
	localDf.to_csv(filepath, sep = '\t', index=False)

def reheader_samples_with_age(mutationsMafPath='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impactMafs/data_mutations_unfiltered_as_of_may_16.txt'):
	def get_cases(path):
		f = open(path)
		lines = f.readlines()
		return set(lines[4].split('\t'))

	def parse_age_df_into_brackets(ageDf):
		peds = set()
		young = set()
		middle = set()
		old = set()
		for index, row in ageDf.iterrows():
			if int(row['AGE']) < 17: peds.add(row['PATIENT_ID'])
			elif int(row['AGE']) < 40: young.add(row['PATIENT_ID'])
			elif int(row['AGE']) < 50: middle.add(row['PATIENT_ID'])
			else: old.add(row['PATIENT_ID'])
		return peds, young, middle, old

	def logicFunction(x, youngSet, middleSet, oldSet): #to avoid trying to write a non trivial lambda function cause syntax is confusing
		#print len(young), len(middle), len(old)

		if x in oldSet:
			return 'old'
		elif x in middleSet:
			return 'middle'
		elif x in youngSet: return 'young'
		else:
			return 'none'

	def get_pole_cases(mutSigsDf):
		poleSigPositive = mutSigsDf[mutSigsDf['Signature.10'] > .2]
		return set(poleSigPositive['Sample Name'])

	def get_mmr_cases(mutSigsDf):
		mmrPositive = mutSigsDf[(mutSigsDf['Signature.6'] + mutSigsDf['Signature.15'] + mutSigsDf['Signature.20'] + mutSigsDf['Signature.26']) > .2]
		return set(mmrPositive['Sample Name'])

	darwinAgePath = '/home/friedman/friedman/clinicalData/msk-impact/msk-impact/darwin/darwin_age.txt'
	ageDf = pd.read_table(darwinAgePath)
	peds, young, middle, old = parse_age_df_into_brackets(ageDf)
	mutSigsDf = pd.read_table('/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mutationSigFiles/may16unfiltered30sigs.txt')
	mutSigsDf = mutSigsDf[mutSigsDf['Number of Mutations'] < 50] #only look at cases with fewer than 50 mutations
	caseListDir = '/home/friedman/friedman/msk-impact/msk-impact/case_lists'
	mutationsDf = pd.read_table(mutationsMafPath, skiprows=[0])
	mutationsDf['PatientID'] = mutationsDf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
	mutationsDf = mutationsDf[(mutationsDf['PatientID'].isin(young)) | (mutationsDf['PatientID'].isin(middle)) | (mutationsDf['PatientID'].isin(old))] #only keep stuff we have info for
	mutationsDf['Age'] = mutationsDf['PatientID'].apply(lambda x: logicFunction(x, young, middle, old))
	mutationsDf['newCol'] = mutationsDf['Tumor_Sample_Barcode']
	for file in os.listdir(caseListDir)[:1]:
		print file
		filepath = os.path.join(caseListDir, file)
		caseType = file[10:].strip('.txt')
		currentTypeCases = get_cases(filepath)
		#ALERT WHAT TO WORK ON TOMORROW SOMEHOW THIS CASE TYPE INFORMATION DOESNT APPROPRIATELY PROPAGATE
		mutationsDf['newCol'] = mutationsDf.apply(lambda row: row['newCol'] if row['Tumor_Sample_Barcode'] not in currentTypeCases else caseType + '_' + logicFunction(row['Age'], young, middle, old), axis=1)
		
	mutationsDf['Tumor_Sample_Barcode'] = mutationsDf['newCol']
	print 'writing to file'
	mutationsDf.to_csv('allMutationsWithAgeAndCancerSubtypeAsIdentifier.maf', sep='\t', index=False)

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--mafPath', help='path to initial maf', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impactMafs/data_mutations_extended_mafAnno.maf')
	parser.add_argument('--mode', help='mode to run the script in')

	args = parser.parse_args()

	reheader_samples_with_age()
	sys.exit()

	#amalgamatedMaf = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mergedAnnotations.tsv'
	#amalgamatedDf = pd.read_table(amalgamatedMaf)

	#d = map_barcode_to_cancer_type(amalgamatedDf)
	
	df = pd.read_table(args.mafPath, skiprows=[0])

	#reheader_with_hugo_symbol_mmr_specific(df)
	
	#reheader_with_hugo_symbol_pole_specific(df)

	reheader_with_hugo_symbol_tmz_specific(df)

	#df = remove_multiple_mutations(df)




if __name__ == '__main__':
    main()















