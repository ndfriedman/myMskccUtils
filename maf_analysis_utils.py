#written by Noah Friedman
#funcions for performing analyses and functions on mafs 
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'

import imp
analysis_utils = imp.load_source('analysis_utils', pathPrefix + '/ifs/work/taylorlab/friedman/myUtils/analysis_utils.py')

#OMNIBUS maf prepping function (adds tid, pid etc)
def prep_maf(df,
			sampleNameCol='Tumor_Sample_Barcode',
			pid=True, tid=True, quadNuc=True, cancerType=True,
			secondSig=None):

	if tid: df['tid'] = df[sampleNameCol].apply(lambda x: x[:13])
	if pid: df['pid'] = df[sampleNameCol].apply(lambda x: x[:9])
	if quadNuc:
		df['quadNuc'] = df.apply(lambda row: 
   			mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'])
    		if row['Variant_Type'] == 'SNP' else None, axis=1)
	if cancerType:
		cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
		df['cancer_type'] = df['pid'].apply(lambda x: cDict[x] if x in cDict else None)
	
	if secondSig != None: #mark the second most common signature
		pass
	return df

#a utility function designed to add information about hotspots (which hotspot, what 4nuc etc) to a df with one entry per case
def add_hotspot_maf_info_to_df(df, mafDf,
 allHotspots, hotspotsToFocusOn, #a set of all hospots and a set of hotspots to focus on (those not to call 'other')
 idCol = 'pid', removeDups=True):
	hotspotsWeCareAbout = mafDf[mafDf['Amino_Acid_Change'].isin(allHotspots)]
	if removeDups:
		hotspotsWeCareAbout = hotspotsWeCareAbout.drop_duplicates(subset=[idCol], keep=False)
	hotspotLabelDict = dict(zip(hotspotsWeCareAbout[idCol], hotspotsWeCareAbout['Amino_Acid_Change']))

	hotspotFourNucDict = dict(zip(hotspotsWeCareAbout[idCol], hotspotsWeCareAbout['quadNuc']))
	df['Amino_Acid_Change'] = df[idCol].apply(lambda x: hotspotLabelDict[x] if x in hotspotLabelDict else None)
	df['Amino_Acid_Change_Adj'] = df[idCol].apply(lambda x:
		hotspotLabelDict[x] if (x in hotspotLabelDict and hotspotLabelDict[x] in hotspotsToFocusOn) else 'otherHotspot' if x in hotspotLabelDict else None)
	df['quadNuc'] = df[idCol].apply(lambda x: hotspotFourNucDict[x] if x in hotspotFourNucDict else None)
	return df


def summarize_per_case_mutation_info_for_mafs(mafDf):
	oncoKbOncogenicAnnotations = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	listOfDicts = []
	cases = set(mafDf['Tumor_Sample_Barcode'])
	cntr = 0
	for case in cases:
		if cntr%20==0: #slow function so I print progress
			print cntr, len(cases)
		localDict = {}
		caseDf = mafDf[mafDf['Tumor_Sample_Barcode'] == case]
		nOncogenicMuts = caseDf[caseDf['oncogenic'].isin(oncoKbOncogenicAnnotations)].shape[0]
		nHotspotMuts = caseDf[caseDf['is-a-hotspot'] == 'Y'].shape[0]
		nActivatingMutations = caseDf[(caseDf['is-a-hotspot'] == 'Y') | (caseDf['oncogenic'].isin(oncoKbOncogenicAnnotations))].shape[0]
		nMut = caseDf.shape[0]

		localDict['Nmut'] = nMut
		localDict['nOncogenicMuts'] = nOncogenicMuts
		localDict['nHotspotMuts'] = nHotspotMuts
		localDict['nActivatingMutations'] = nActivatingMutations

		listOfDicts.append(localDict)
		cntr += 1

	return pd.DataFrame(listOfDicts)

#a function that enumerates all activating mutations in a gene across the cohort
def enumerate_activating_muts_across_cohort(gene, mafDf = None):

	#TODO appropriately allow the user to load the maf df if there is no default
	#if mafDf.sh:
	#	mafDf = pd.read_table('/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')

	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	geneOncogenicOrHotspotMuts = mafDf[(mafDf['Hugo_Symbol'] == gene) &((mafDf['is-a-hotspot'] == 'Y') |(mafDf['oncogenic'].isin(oncogenicMutColNames)))]

	colsToKeep = ['Tumor_Sample_Barcode', 'is-a-hotspot', 'oncogenic', 'Ref_Tri', 'Tumor_Seq_Allele2', 'Reference_Allele']
	return geneOncogenicOrHotspotMuts[[colsToKeep]]

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()

if __name__ == '__main__':
    main()















