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


####little utilities
def fix_mll_genes(maf):
    maf['Hugo_Symbol'] = maf['Hugo_Symbol'].apply(lambda x:
        'KMT2A' if x == 'MLL'
        else 'KMT2B' if x == 'MLL2'
        else 'KMT2C' if x == 'MLL3'
        else x)   
    return maf



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


#a utility function that we use to add a column to a dataframe that summarizes True/False for the column X which is true if its hotspot/oncogenic/truncating
def add_mut_effect_summary_col(mafDf):
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	mafDf['mutKnowOncogenicOrTruncating'] = mafDf.apply(lambda row: 
		True if row['is-a-hotspot'] == 'Y'
		else True if row['oncogenic'] in oncogenicMutColNames
		else True if row['Consequence'] == 'stop_gained'
		else True if row['Consequence'] == 'frameshift_variant' #TODO is this something that should be included
		else False, 
		axis=1)
	return mafDf


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
		snpsOnly = caseDf[caseDf['Variant_Type'] == 'SNP']
		nSnps = snpsOnly.shape[0]
		nOncogenicSnps = snpsOnly[snpsOnly['oncogenic'].isin(oncoKbOncogenicAnnotations)].shape[0]

		localDict['Tumor_Sample_Barcode'] = case
		localDict['Nmut'] = nMut
		localDict['NSnps'] = nSnps
		localDict['nOncogenicSnps'] = nOncogenicSnps
		localDict['nOncogenicMuts'] = nOncogenicMuts
		localDict['nHotspotMuts'] = nHotspotMuts
		localDict['nActivatingMutations'] = nActivatingMutations

		listOfDicts.append(localDict)
		cntr += 1

	return pd.DataFrame(listOfDicts)

#a function that enumerates all activating mutations in a gene across a cohort maf
def enumerate_activating_muts_across_cohort(gene, mafDf = None):

	#TODO appropriately allow the user to load the maf df if there is no default
	#if mafDf.sh:
	#	mafDf = pd.read_table('/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')

	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	geneOncogenicOrHotspotMuts = mafDf[(mafDf['Hugo_Symbol'] == gene) &((mafDf['is-a-hotspot'] == 'Y') |(mafDf['oncogenic'].isin(oncogenicMutColNames)))]

	#colsToKeep = ['Tumor_Sample_Barcode', 'is-a-hotspot', 'oncogenic', 'Ref_Tri', 'Tumor_Seq_Allele2', 'Reference_Allele']
	#return geneOncogenicOrHotspotMuts[[colsToKeep]]
	return geneOncogenicOrHotspotMuts

#a function that given a maf of mutaitons and quadnuc spectra enriched for a signature returns the number of mutations that occur at the enriched spectra
def summarize_signature_attribution_for_case(mafDf, enrichedSigMotifs):
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])

	listOfDicts = []
	cases = set(mafDf['Tumor_Sample_Barcode'])
	for case in cases:
		localD = {}
		caseDf = mafDf[mafDf['Tumor_Sample_Barcode'] == case]

		#get a maf for the current cases hotspot and oncogenic muts 
		hotspotMutDf = caseDf[caseDf['is-a-hotspot'] == 'Y']
		oncogenicMutDf = caseDf[caseDf['oncogenic'].isin(oncogenicMutColNames)]
		nHotspotMutations = hotspotMutDf.shape[0]
		nOncogenicMutations = oncogenicMutDf.shape[0]

		#get data about how many of each occurs at a motif
		nOncogenicMutationsAtEnrichedMotif = oncogenicMutDf[oncogenicMutDf['quadNuc'].isin(enrichedSigMotifs)].shape[0]
		nHotpsotMutationsAtEnrichedMotif = hotspotMutDf[hotspotMutDf['quadNuc'].isin(enrichedSigMotifs)].shape[0]
		nMut = caseDf.shape[0]

		#append info to the local dictionary (row of the future df)
		localD['Tumor_Sample_Barcode'] = case
		localD['nHotspots'] = nHotspotMutations
		localD['Nmut'] = nMut
		localD['nOncogenicMutations'] = nOncogenicMutations
		localD['nOncogenicMutationsAtEnrichedMotif'] = nOncogenicMutationsAtEnrichedMotif
		localD['nHotpsotMutationsAtEnrichedMotif'] = nHotpsotMutationsAtEnrichedMotif  

		listOfDicts.append(localD)

	df = pd.DataFrame(listOfDicts)
	return df      

#returns a dataframe with the cases that have mutations in a gene and the type of mutation in that gene
#note the way we return it is designed for the R long format
def asses_per_case_mut_info_for_gene(mafDf, gene, quadnucSet):

	def classify_mut_residue_data(quadNucs, qNucSet):
		mutType = None
		if len(quadNucs) == 1: #if there is only one quadnuc mutation figure out which type it is
			v = quadNucs.pop()
			if v in qNucSet: return 'favoredMutation'
			else: return 'notFavoredMutation'
		else: #if there are more than one mutations in the gene we need to classify whether the mutations are mixed, from the favored process only or from both
			nFavoredMutations = 0
			nNotFavoredMutations = 0
			for v in quadNucs:
				if v in qNucSet:
					nFavoredMutations += 1
				else:
					nNotFavoredMutations += 1
				if nFavoredMutations >= 1 and nNotFavoredMutations >= 1:
					return 'mixed'
				elif nFavoredMutations >= 1:
					return 'mutlipleFavoredMutation'
				else:
					return 'multipleNotFavoredMutation'

	cases = set(mafDf['Tumor_Sample_Barcode'])
	geneMuts = mafDf[mafDf['Hugo_Symbol'] == gene]
	#we only put information in the dataframe we return if there are snp muts
	listOfDicts = []
	for case in cases:
		caseMuts = geneMuts[geneMuts['Tumor_Sample_Barcode'] == case]
		caseQuadNucs = set(caseMuts['quadNuc'])
		if len(caseQuadNucs) != 0:
			if None not in caseQuadNucs:
				classification = classify_mut_residue_data(caseQuadNucs, quadnucSet)
				localD = dict()
				localD['gene'] = gene
				localD['Tumor_Sample_Barcode'] = case
				localD['mutClassification'] = classification
				localD['nGeneMut'] = len(caseQuadNucs)
				listOfDicts.append(localD)
	return pd.DataFrame(listOfDicts)

#a utility that asseses the SNP burden for each case
def asses_snp_burden_across_cohort(maf):
	cases = set(maf['Tumor_Sample_Barcode'])
	cntr = 0
	listOfDicts = []
	for case in cases:
		if cntr%500 == 0: print cntr, len(cases)
		localDict = dict()
		cntr += 1

		caseMuts = maf[maf['Tumor_Sample_Barcode'] == case]
		caseSnps = caseMuts[caseMuts['Variant_Type'] == 'SNP']
		caseIndels = caseMuts[(caseMuts['Variant_Type'] == 'INS') | (caseMuts['Variant_Type'] == 'DEL')]
		localDict['Tumor_Sample_Barcode'] = case
		localDict['nSnps'] = caseSnps.shape[0]
		localDict['nIndels'] = caseIndels.shape[0]
		localDict['nMuts'] = caseMuts.shape[0]
		oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
		localDict['nOncogenicMutations'] = caseMuts[caseMuts['oncogenic'].isin(oncogenicMutColNames)].shape[0]
		listOfDicts.append(localDict)

	return pd.DataFrame(listOfDicts)

#a utility function that enumerates the top N genes with the most oncogenic mutations in distinct samples
#also enumerates fraction of cases with muts in each case
#returns three things:
#i. A counter of N cases with oncogenic muts
#ii. A dict with gene name: cohort ranking
#iii. A dict mapping gene to fraction of cohort with an oncogenic mutation in that case
def enumerate_top_n_oncogenic_mutated_genes_across_cohort(cohortMaf, n=None):
	oncoKbOncogenicAnnotations = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	oncogenicMutations = cohortMaf[cohortMaf['oncogenic'].isin(oncoKbOncogenicAnnotations)]
	oncogenicMutations['patientGeneMutated'] = oncogenicMutations.apply(lambda row: row['pid'] + '_' + row['Hugo_Symbol'], axis=1)
	oncogenicMutationsSansPatientDuplicates = oncogenicMutations.drop_duplicates(subset=['patientGeneMutated']) #NOTE this analysis isnt perfect as it treats independent primaries as genetically related
	
	nPatients = len(set(cohortMaf['pid']))	
	occurenceCounter = None
	if n!= None:
		occurenceCounter = Counter(oncogenicMutationsSansPatientDuplicates['Hugo_Symbol']).most_common(n)
	else:
		occurenceCounter = Counter(oncogenicMutationsSansPatientDuplicates['Hugo_Symbol']).most_common()

	fractionalDict = {key: value for (key, value) in occurenceCounter}
	for key, value in fractionalDict.items():
		fractionalDict[key] = 1.0*value/nPatients

	rankingDict = dict()
	cntr = 1 #counter is 1 indexed is that a problem
	for gene in occurenceCounter:
		rankingDict[gene[0]] = cntr
		cntr += 1

	return occurenceCounter, rankingDict, fractionalDict

######################LOH/CNA analysis tools

def do_gene_loh_summary(maf, genes=None):
	listOfDicts = []

	if genes == None: 
		genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])

	for case in set(maf['Tumor_Sample_Barcode']):
		localD = {}
		localD = {'Tumor_Sample_Barcode': case}
		for gene in genes:

			mutAtGene = False
			oneMutOncogenic = False
			lohAtGene = False
			caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
			geneMaf = caseMaf[caseMaf['Hugo_Symbol'] == gene]
			if geneMaf.shape[0] > 0:
				mutAtGene = True
			if geneMaf[geneMaf['oncogenic'].isin(oncogenicMutColNames)].shape[0] > 0:
				oneMutOncogenic = True
			if geneMaf[geneMaf['lcn'] == 0].shape[0] > 0:
				lohAtGene = True
			lohPlusOncogenic = False
			if oneMutOncogenic & lohAtGene:
				lohPlusOncogenic = True

			
			localD[gene + '_mut'] = mutAtGene
			localD[gene + '_loh'] = lohAtGene
			localD[gene + '_oncogenicMut'] = oneMutOncogenic
			localD[gene + '_lohPlusOncogenic'] = lohPlusOncogenic
		listOfDicts.append(localD)

	return pd.DataFrame(listOfDicts)

#a function that given a cohort maf enumerates which tumor suppressors and oncogenes are reccurently mutated (mutated in greater than 'thresh' fraction of the cohort)
def enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(cohortMaf, thresh=.1):
    tumorSupressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])
    occurenceCounter, rankingDict, fractionalDict = enumerate_top_n_oncogenic_mutated_genes_across_cohort(cohortMaf, n=50)
    recurrentTumorSupressors = []
    recurrentOncogenes = []
    for key, value in fractionalDict.items():
        if value > thresh:
            if key in tumorSupressors:
                recurrentTumorSupressors.append(key)
            else:
                recurrentOncogenes.append(key) 
    return set(recurrentTumorSupressors), set(recurrentOncogenes)


def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()

if __name__ == '__main__':
    main()















