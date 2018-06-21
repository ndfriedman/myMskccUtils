#written by Noah Friedman 
#little utility script to do a la carte subsetting of mafs

import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')


def get_cases_for_cancer_type(cancerType, cancerTypeListDir = '/ifs/work/taylorlab/friedman/msk-impact/msk-impact/case_lists'):
	cancerType = 'case_list_' + cancerType + '.txt'
	path = os.path.join(cancerTypeListDir, cancerType)
	f = open(path)
	lines = f.readlines()
	return set(lines[4].split('\t'))

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--mode', help='mode to run the script in', default='subsetByCancerType')
	parser.add_argument('--cancerTypes', default='')
	parser.add_argument('--inputMaf', default='/ifs/work/taylorlab/friedman/myUtils/FilteredMafWithHospotAndSignatures.maf')
	parser.add_argument('--inputMaf2', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impactMafs/data_mutations_extended_mafAnno.maf')
	parser.add_argument('--outputFilename', default='adjustedMaf.maf')
	parser.add_argument('--outputDir', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/adjustedSubsetMafs')
	args = parser.parse_args()

	mafDf = pd.read_table(args.inputMaf, skiprows=[0])
	writePath = os.path.join(args.outputDir, args.outputFilename)

	if args.mode == 'subsetByCancerType':
		cases = get_cases_for_cancer_type('Bladder_Cancer')
		mafDf = mafDf[mafDf['Tumor_Sample_Barcode'].isin(cases)]
		print 'writing file to ', writePath
		mafDf.to_csv(writePath, sep='\t', index=False)

	if args.mode == 'mashMafs': #mode to merge two mafs
		mafDf1 = pd.read_table(args.inputMaf)
		mafDf2 = pd.read_table(args.inputMaf2, skiprows=[0])

		mafDf1['idCol'] = mafDf1.apply(lambda row: str(row['Chromosome']) + '_' + str(row['Start_Position']), axis=1)
		mafDf2['idCol'] = mafDf2.apply(lambda row: str(row['Chromosome']) + '_' + str(row['Start_Position']), axis=1)

		#print mafDf1['idCol']
		#print mafDf2['idCol']
		print len(list(mafDf1['idCol']))
		print len(list(mafDf2['idCol']))
		#mergedDf = mafDf1.merge(mafDf2, how='left', left_on='idCol', right_on='idCol')
		mergedDf = mafDf1.merge(mafDf2, on='idCol')

		#print mergedDf.shape
		#print mergedDf.columns.values
		#print mergedDf['idCol']
		#print mergedDf['Start_Position_x']
		#print mergedDf['Start_Position_y']
		print len(set(mafDf1['idCol']) & set(mafDf2['idCol']))

		mergedDf.to_csv('facetsAndTrinucAndHotspotsAndSignatures.maf', index=False, sep='\t')

if __name__ == '__main__':
    main()















