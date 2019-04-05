#written by Noah Friedman
#a collection of utilities designed to model mutation acquisition in silico 
import sys
import argparse
import os
import pandas as pd
import numpy as np
import math

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')
import variant_consequence_analysis_util
import mutationSigUtils


def init_gene_probabilities():
	return 0

#CURRENTLY DEPRECTED BUT IS A GOOD IDE A FOR WHAT SHIT WOULD LOOK LIKE IF WE RANDOMLY PICK SIGS
def get_trinucleotide_probabilities(trinucEntries):
		rowAsDict = trinucEntries.to_dict()
		trinucChoices = rowAsDict.keys()
		trinucCounts = rowAsDict.values()
		trinucProbabilities = [1.0*i/sum(trinucCounts) for i in trinucCounts]
		return trinucChoices, trinucProbabilities

def normalize_counter(cntrObj, mode='Round', nDigitsRound=2):
	total = sum(cntrObj.values(), 0.0)
	for key in cntrObj:
		cntrObj[key] /= total
		if mode == 'Round': cntrObj[key] = round(cntrObj[key], nDigitsRound)
	return cntrObj

####################################

###################################SAMPLING AREA

#note this could be done a lot faster as a matrix multiplication; todo implement it that way?
def get_quadnuc_fracs_given_decomposition(caseDecomposition, spectraPath = '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt'):
	spectraD = mutationSigUtils.convert_spectrum_file_to_dict_of_dicts(spectrumFile=spectraPath)
	quadNucs = spectraD['Signature.1'].keys()
	quadNucD = {} #initialize the fractions for each quadnuc
	for quadNuc in quadNucs:
		quadNucD[quadNuc] = 0

	for key, value in spectraD.items():
		decompositionFraction = caseDecomposition[key]
		for quadNuc, percent in value.items():
			quadNucD[quadNuc] = quadNucD[quadNuc] + decompositionFraction*percent

	#RETURNS a dictionary for each quadnuc, summing up to one
	return quadNucD
			

###################################OVERALL MODELING/COMPARISSON

def get_average_oncogenic_mut_burdens_for_tcga(caseDecompositionDict, 
	spectraD, 
	quadNucDict,
	quadNucDictAdjusted,
	geneOncogenicOrHotspotProbDict,
	geneMutProbInfo,
	n=1000,
	mutMode='oncogenic'):

	oncogenicMutFracs = []
	for model in ['model1', 'model2', 'model3']:

		nNonImpactMuts, spectraChosen, quadNucsChosen, genesChosen, mutationsChosen = pick_n_tcga_mutations_given_spectra(
		caseDecompositionDict, 
		spectraD, 
		quadNucDict,
		quadNucDictAdjusted,
		geneOncogenicOrHotspotProbDict,
		geneMutProbInfo,
		n=n,
		mutMode=mutMode,
		modelName=model)

		nNonOncogenicMuts = len([i for i in mutationsChosen if i == 'nonOncogenic'])
		nNotOncogenic = nNonImpactMuts + nNonOncogenicMuts
		nOncogenic = n - nNotOncogenic
		fracOncogenic = nOncogenic*1.0/n
		oncogenicMutFracs.append(fracOncogenic)
	averageFrac = np.nanmean(oncogenicMutFracs)
	return averageFrac

def get_simulated_oncogenic_mut_frac_per_case(
	signaturesDf,
	spectraD, 
	quadNucDict,
	quadNucDictAdjusted,
	geneOncogenicOrHotspotProbDict,
	geneMutProbInfo,
	n=1000,
	mutMode='oncogenic',
	idColName='Tumor_Sample_Barcode'):
	
	mutFracDict = {}
	cases = set(signaturesDf[idColName])
	signatureColNames = ['Signature.' + str(i) for i in range(1,31)]
	cntr = 0
	for case in cases:
		if cntr%50 == 0: print cntr
		cntr += 1

		localDf = signaturesDf[signaturesDf[idColName] == case]
		caseSignatures = localDf[signatureColNames].iloc[0].to_dict()
		avgVal = get_average_oncogenic_mut_burdens_for_tcga(
			caseSignatures, 
			spectraD, 
			quadNucDict,
			quadNucDictAdjusted,
			geneOncogenicOrHotspotProbDict,
			geneMutProbInfo,
			n=n,
			mutMode=mutMode)

		mutFracDict[case] = avgVal

	return mutFracDict


def compare_observed_to_expected_mut_burden_data(actualMutInfoDf, simulatedMutInfo, idCol='Tumor_Sample_Barcode'):
	listOfDicts = []
	for caseId, value in simulatedMutInfo.items():
		localD = {}
		acutalMutCaseInfo = actualMutInfoDf[actualMutInfoDf[idCol] == caseId]
		if acutalMutCaseInfo.shape[0] > 0:
			actualNOncogenicSnps = acutalMutCaseInfo.iloc[0]['nOncogenicSnps']
			actualNmut = acutalMutCaseInfo.iloc[0]['Nmut']
			simulatedNOncogenicSNPs = actualNmut*value
			localD['Tumor_Sample_Barcode'] = caseId
			localD['Nmut'] = actualNmut
			localD['nOncogenicSnps'] = actualNOncogenicSnps
			localD['simulatedNOncogenicSNPs'] = simulatedNOncogenicSNPs
			if actualNOncogenicSnps > 0:
				localD['ratioObservedToExpected'] = math.log(1.0*actualNOncogenicSnps/simulatedNOncogenicSNPs, 2)
			else:
				localD['ratioObservedToExpected'] = None
			listOfDicts.append(localD)
	return pd.DataFrame(listOfDicts)

###################################################FUNCTIONS

def calculate_quadnuc_based_oncogenic_susceptibility_dict(df):
	d = {}
	quadNucs = set(filter(lambda x: x == x , set(df['quadNuc']))) #dont allow me to consider nas

	impact368genes = set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1'])
	impact410genes = set(['ABL1', 'ACVR1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VEGFA', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	impact468genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']) 
	
	dfImpact368 = df[df['Hugo_Symbol'].isin(impact368genes)]
	dfImpact410 = df[df['Hugo_Symbol'].isin(impact410genes)]
	dfImpact468 = df[df['Hugo_Symbol'].isin(impact468genes)]
	
	for quadnuc in quadNucs:
		localD = {}
		for impactVersion in ['IMPACT_468', 'IMPACT_410', 'IMPACT_368']:
			if impactVersion == 'IMPACT_368':
				oncogenicMutSum = sum(dfImpact368[dfImpact368['quadNuc'] == quadnuc]['nOncogenicMut'])
				totalMutSum = sum(dfImpact368[dfImpact368['quadNuc'] == quadnuc]['totalNmut'])
				if totalMutSum == 0:
					frac = 0
				else:
					frac = 1.0*oncogenicMutSum/totalMutSum
				localD[impactVersion] = frac
			elif impactVersion == 'IMPACT_410':
				oncogenicMutSum = sum(dfImpact410[dfImpact410['quadNuc'] == quadnuc]['nOncogenicMut'])
				totalMutSum = sum(dfImpact410[dfImpact410['quadNuc'] == quadnuc]['totalNmut'])
				if totalMutSum == 0:
					frac = 0
				else:
					frac = 1.0*oncogenicMutSum/totalMutSum
				localD[impactVersion] = frac
			else:
				oncogenicMutSum = sum(dfImpact468[dfImpact468['quadNuc'] == quadnuc]['nOncogenicMut'])
				totalMutSum = sum(dfImpact468[dfImpact468['quadNuc'] == quadnuc]['totalNmut'])
				if totalMutSum == 0:
					frac = 0
				else:
					frac = 1.0*oncogenicMutSum/totalMutSum
				localD[impactVersion] = frac
		d[quadnuc] = localD
	return d

#calculates the probability that an oncogenic mutation occurs at a gene given a quad nuc
#NOTE MAYBE write differently so its not just a repeat of the above fucntion
def calculate_quadnuc_and_gene_based_oncogenic_susceptibility_dict(df):

	d = {}

	impact368genes = set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1'])
	impact410genes = set(['ABL1', 'ACVR1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VEGFA', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	impact468genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']) 
	
	dfImpact368 = df[df['Hugo_Symbol'].isin(impact368genes)]
	dfImpact410 = df[df['Hugo_Symbol'].isin(impact410genes)]
	dfImpact468 = df[df['Hugo_Symbol'].isin(impact468genes)]

	quadNucs = set([start + change + end for start in ['A', 'C', 'T', 'G'] 
                        for change in ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
                        for end in ['A', 'C', 'T', 'G']])

	cntr = 0
	for gene in impact468genes:
		cntr += 1
		print 'summarizing gene ', gene, ' which is the ', cntr, 'th gene out of 468'
		geneD = {}
		for quadnuc in quadNucs:
			localD = {}
			for impactVersion in ['IMPACT_468', 'IMPACT_410', 'IMPACT_368']:
				if impactVersion == 'IMPACT_368':
					oncogenicMutSum = sum(dfImpact368[(dfImpact368['quadNuc'] == quadnuc) & (dfImpact368['Hugo_Symbol'] == gene)]['nOncogenicMut'])
					totalMutSum = sum(dfImpact368[dfImpact368['quadNuc'] == quadnuc]['totalNmut'])
					if totalMutSum == 0:
						frac = 0
					else:
						frac = 1.0*oncogenicMutSum/totalMutSum
					localD[impactVersion] = frac
				elif impactVersion == 'IMPACT_410':
					oncogenicMutSum = sum(dfImpact410[(dfImpact410['quadNuc'] == quadnuc) & (dfImpact410['Hugo_Symbol'] == gene)]['nOncogenicMut'])
					totalMutSum = sum(dfImpact410[dfImpact410['quadNuc'] == quadnuc]['totalNmut'])
					if totalMutSum == 0:
						frac = 0
					else:
						frac = 1.0*oncogenicMutSum/totalMutSum
					localD[impactVersion] = frac
				else:
					oncogenicMutSum = sum(dfImpact468[(dfImpact468['quadNuc'] == quadnuc) & (dfImpact468['Hugo_Symbol'] == gene)]['nOncogenicMut'])
					totalMutSum = sum(dfImpact468[dfImpact468['quadNuc'] == quadnuc]['totalNmut'])
					if totalMutSum == 0:
						frac = 0
					else:
						frac = 1.0*oncogenicMutSum/totalMutSum
					localD[impactVersion] = frac
			geneD[quadnuc] = localD
		d[gene] = geneD
	return d

	#TODO implement

def get_expected_oncogenic_val_given_quadnuc_fractions(quadNucFractions, susceptibilityDict, impactVersion):
	p = 0
	for key, value in quadNucFractions.items():
		p += value * susceptibilityDict[key][impactVersion]
	return p

#SAME FUNCTION AS ABOVE 
def get_expected_oncogenic_val_given_quadnuc_fractions_and_gene(quadNucFractions, susceptibilityDict, impactVersion, gene):
	p = 0
	for key, value in quadNucFractions.items():
		p += value * susceptibilityDict[gene][key][impactVersion]
	return p

def compare_observed_and_expected_gene_mutation_rates(susceptibilityDf, observedDf):
	for gene in set(susceptibilityDf['Hugo_Symbol']):
		print gene

###&&&&&&&&*************New code 2-18-2019!!!
def get_chance_of_oncogenic_mutation_at_quadnuc(df, quadnuc):
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	allMutsAtQuadNuc = df[df['quadNuc'] == quadnuc]
	nMutsAtQuadNuc = sum(allMutsAtQuadNuc['totalNmut'])
	nOncogenicMutAtQuadNuc = sum(allMutsAtQuadNuc['nOncogenicMut'])
	return 1.0*nOncogenicMutAtQuadNuc/nMutsAtQuadNuc

def get_chance_of_oncogenic_mutation_in_gene_at_quadnuc(df, quadnuc, gene):
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	allMutsAtQuadNuc = df[df['quadNuc'] == quadnuc]
	nMutsAtQuadNuc = sum(allMutsAtQuadNuc['totalNmut'])
	nOncogenicMutAtQuadNuc = sum(allMutsAtQuadNuc[allMutsAtQuadNuc['Hugo_Symbol'] == gene]['nOncogenicMut'])
	return 1.0*nOncogenicMutAtQuadNuc/nMutsAtQuadNuc

#given a decomposition of a Sample expressed as a dict (ie 50%sig1, 25% sig2, 25%sig3) randomly picks a signature to have caused a mutation
def pick_spectra_given_decomposition(decompositionDict):
	return np.random.choice(
			decompositionDict.keys(), p=decompositionDict.values()
		)

#One liner to pick a quadnuc given a signature that theoretically genrted it
def pick_quadnuc_given_spectra(spectraDict):
	return np.random.choice(
			spectraDict.keys(), p=spectraDict.values()
		)

#one liner to return a gene that is mutated given a quadNuc
def pick_gene_given_quadnuc(d, quadNuc): 
	return np.random.choice(
			d[quadNuc][0], p=d[quadNuc][1]
		)

#one liner to pick a mutation at random given a 
def pick_mutation_given_gene_and_quadNuc(gene, quadNuc, d, mode='oncogenic'):
	return np.random.choice(
			d[gene + '_' + quadNuc][mode][0],
			p= d[gene + '_' + quadNuc][mode][1]
		)


def asses_simulated_mutation_rates_by_signatures(sigMagnitudesDict, k, spectraD, n, quadNucDict, geneOncogenicOrHotspotProbDict, mMode='oncogenic'):
    d = dict()
    signatures = ['Signature.' + str(i) for i in range(1,31)]
    for sig in signatures:
        sigMagnitudesDict = {key:value for (key, value) in [('Signature.' + str(i), 0) for i in range(1,31)]}
        sigMagnitudesDict[sig] = 1
        d[sig] = mutation_modeling_util.do_k_simulations_of_mutations(sigMagnitudesDict, k, spectraD, n, quadNucDict, geneOncogenicOrHotspotProbDict, mutMode=mMode)
    return d


#given a maf of all possible mutations in impact summarizes for each gene the number of oncogenic/not oncogenic mutations that occur at each quadnuc and the 
def initiate_gene_mut_mapping(allMutsMaf): #TO improve runtime we precalculate a dict mapping 
															#gene: maf of all mutations and maf of oncogenic mutations
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	listOfDicts = []
	allQuadNucs = set(allMutsMaf['quadNuc'])
	for gene in set(allMutsMaf['Hugo_Symbol']):
		print 'analyzing gene ', gene
		geneMuts = allMutsMaf[allMutsMaf['Hugo_Symbol'] == gene]
		chromosome = geneMuts['Chromosome'].iloc[0]
		impactPanel = None

		#Set the proper impact information
		impact368Panel = set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1'])
		impact410Panel = set(['ABL1', 'ACVR1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VEGFA', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
		impact468Panel = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
		if gene in impact368Panel:
			impactPanel = 'IMPACT_368'
		elif gene in impact410Panel:
			impactPanel = 'IMPACT_410'
		elif gene in impact468Panel:
			impactPanel = 'IMPACT_468'
		else:
			impactPanel = None

		for quadNuc in allQuadNucs:
			geneMutsAtQuadNuc = geneMuts[geneMuts['quadNuc'] == quadNuc]
			nMutsInGeneAtQuadNuc = geneMutsAtQuadNuc.shape[0]
			nOncogenicMutsInGeneAtQuadNuc = geneMutsAtQuadNuc[geneMutsAtQuadNuc['oncogenic'].isin(oncogenicMutColNames)].shape[0]
			localD = {'totalNmut': nMutsInGeneAtQuadNuc, 'nOncogenicMut': nOncogenicMutsInGeneAtQuadNuc,
			'Hugo_Symbol': gene,'quadNuc': quadNuc, 'impactPanel': impactPanel, 'Chromosome': chromosome}
			listOfDicts.append(localD)
			
	return pd.DataFrame(listOfDicts)


#an option intialization function that loads in the simulated mutation mafs and summarizes them one at a time
def load_in_and_summarize_simulated_mut_information(mafDirPath = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs'):
	chromosomes = ['X'] + [str(i) for i in range(1,23)]
	chromosomes.reverse() #reverse the order to do small stuff first for debugging purposes
	listOfDataFrames = []

	for chromosome in chromosomes:
		fileName = 'chr' + chromosome + '_allmut_mafAnno_trinuc.maf'
		path = os.path.join(mafDirPath, fileName)
		print 'loading in file ', path
		maf = pd.read_table(path)
		print 'adding quadNuc to chr', chromosome
		maf['quadNuc'] = maf.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
		print 'sumarizing chr', chromosome
		mutMappingDf = initiate_gene_mut_mapping(maf)
		listOfDataFrames.append(mutMappingDf)
		
	concatDf = pd.concat(listOfDataFrames)
	return concatDf 


#TODO: write a function for reading probabilitiy data from disk and converting it into a dict of the appropriate format


#$#$$$#$#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#OLD CODE THAT IS NO LONGER IN USE
#GIVEN a mutation signature spectra for a case, choose N mutations
#Model 1: the gene mutated depends on the quadnuc mutated
#Model2: The Gene mutated and quadnuc mutated are independent
#model3: the gene mutated depends on the quadnuc mutated with an adjustment for the gene mutation rates

"""def pick_n_mutations_given_spectra(caseDecompositionDict, 
	allSpectraDict, #a dictionary of all spectras from which we get a mutation signature spectra
	nMutations, #the number of mutations we are supposed to simulate
	qNucDict, #probability of picking a specific quadnuc given a spectra 
	qNucDictAdj,#probability of gene mutation given quadnuc weighted by observed mutation rates across impact genes
	geneOncogenicOrHotspotDict, #probability of an oncogenic mutation O given a mutation in a gene G and a quadnuc Q
	geneMutProbInfo, #The probability of mutating a gene (based on the probabilities of genes being mutated across all silent mutations)
	mutMode='oncogenic', #mutMode--do we care about picking oncogenic mutations or picking hotspot mutations
	modelName='model1' #model: specifies the assumptions we use to pick the mutations
	): 

	quadNucsChosen = []
	genesChosen = []
	mutationsChosen = []
	spectraChosen = []
	for i in range(nMutations):
		generatingSignature = pick_spectra_given_decomposition(caseDecompositionDict)
		spectraChosen.append(generatingSignature)
		spectraDict = allSpectraDict[generatingSignature]

		quadNuc = pick_quadnuc_given_spectra(spectraDict)
		quadNucsChosen.append(quadNuc)

		#THIS IS WHERE THE MODEL MATTERs
		if modelName == 'model2':
			gene = np.random.choice(
			geneMutProbInfo[0],
			p= geneMutProbInfo[1]
			)
			genesChosen.append(gene)
		elif modelName == 'model3':
			gene = pick_gene_given_quadnuc(qNucDictAdj, quadNuc)
			genesChosen.append(gene)
		else:
			gene = pick_gene_given_quadnuc(qNucDict, quadNuc)
			genesChosen.append(gene)

		mutation = pick_mutation_given_gene_and_quadNuc(gene, quadNuc, geneOncogenicOrHotspotDict, mode=mutMode)
		mutationsChosen.append(mutation)

	return spectraChosen, quadNucsChosen, genesChosen, mutationsChosen

#A function that picks n mutation in tcga
#first determines if a function occurs in a tcga or non tcga 
def pick_n_tcga_mutations_given_spectra(
	caseDecompositionDict, 
	spectraD, 
	quadNucDict,
	quadNucDictAdjusted,
	geneOncogenicOrHotspotProbDict,
	geneMutProbInfo,
	mutMode='oncogenic',
	modelName='model1',
	n=100,
	mutInImpactProb=1.6/30, 
	):

	#FIRST WE NEED TO SEE HOW MANY MUTATIONS HIT THE IMPACT REGION VS HOW MANY HIT THE OTHER PART OF THE EXOME
	mutAreas = []
	for i in range(n):
		mutArea = np.random.choice(
				['IMPACT', 'NON-IMPACT'],
				p= [mutInImpactProb, 1-mutInImpactProb]
				)
		mutAreas.append(mutArea)
	mutAreaCounts = Counter(mutAreas)
	nImpactMuts = mutAreaCounts['IMPACT']
	nNonImpactMuts = mutAreaCounts['NON-IMPACT']
	
	spectraChosen, quadNucsChosen, genesChosen, mutationsChosen = pick_n_mutations_given_spectra(
			caseDecompositionDict, 
			spectraD,
			nImpactMuts, 
			quadNucDict,
			quadNucDictAdjusted,
			geneOncogenicOrHotspotProbDict,
			geneMutProbInfo,
			mutMode='oncogenic',
			modelName='model1'
			)

	return nNonImpactMuts, spectraChosen, quadNucsChosen, genesChosen, mutationsChosen


#perfroms "pick n mutations given spectra k times and returns the results"
def do_k_simulations_of_mutations(k,
	caseDecompositionDict, 
	allSpectraDict, #a dictionary of all spectras from which we get a mutation signature spectra
	nMutations, #the number of mutations we are supposed to simulate
	qNucDict, #probability of picking a specific quadnuc given a spectra 
	qNucDictAdj,#probability of gene mutation given quadnuc weighted by observed mutation rates across impact genes
	geneOncogenicOrHotspotDict, #probability of an oncogenic mutation O given a mutation in a gene G and a quadnuc Q
	geneMutProbInfo, #The probability of mutating a gene (based on the probabilities of genes being mutated across all silent mutations)
	mutMode='oncogenic', #mutMode--do we care about picking oncogenic mutations or picking hotspot mutations
	modelName='model1' #model: specifies the assumptions we use to pick the mutations
	):

	oncMutationsSimulated = [] # a list of lists of the oncogenic mutations that were simulated each k cases
	oncSummary = []
	for i in range(k):
		spectraChosen, quadNucsChosen, genesChosen, mutationsChosen = pick_n_mutations_given_spectra(
			caseDecompositionDict, 
			allSpectraDict,
			nMutations, 
			qNucDict,
			qNucDictAdj,
			geneOncogenicOrHotspotDict,
			geneMutProbInfo,
			mutMode=mutMode,
			modelName=modelName
			)
		oncogenicMuts = [i for i in mutationsChosen if i != 'nonOncogenic']
		nNonOncogenicMuts = Counter(mutationsChosen)['nonOncogenic']
		nOncogenicMuts = len(mutationsChosen) - nNonOncogenicMuts 
		oncSummary.append(nOncogenicMuts)
		oncMutationsSimulated.append(oncogenicMuts)
	return np.nanmean(oncSummary)/nMutations , oncMutationsSimulated #return the fraction of mutations that are oncogenic

###############################Probability Initiation Area

#TODO weight gene mutations by probability of a mutation occuring in a gene

#given a quadnuc, a dict containing all the maf mutations at that gene&quadnuc that are oncogenic/hotspot and the numbner of times a trinuc appears in a gene
#calculate the probability P(oncogenic/hotspot mut in gene | mut at Quadnuc)
def do_oncogenic_mut_initiation(quadNuc, dictEntry, nTimesTrinucAppearsInGene):

	if nTimesTrinucAppearsInGene == 0: #if we pick a mut here something has gone terribly wrong...
		return {
		'hotspot': (['error this trinuc doesnt appear in gene'], [1]),
		'oncogenic': (['error this trinuc doesnt appear in gene'], [1])
		}

	oncogenicGeneMuts, hotspotGeneMuts = dictEntry
	oncogenicMutsAtQuadNuc = oncogenicGeneMuts[oncogenicGeneMuts['quadNuc'] == quadNuc]
	hotspotMutsAtQuadNuc = hotspotGeneMuts[hotspotGeneMuts['quadNuc'] == quadNuc]
	uniqueOncogenicMutsAtQuadNuc = oncogenicMutsAtQuadNuc.drop_duplicates(subset =['uniqueMutationIdentifier'])
	uniqueHotspotMutsAtQuadNuc = hotspotMutsAtQuadNuc.drop_duplicates(subset =['uniqueMutationIdentifier'])
	#return a dict of the form:
	#{'hotspot': (hotspotMutNames, probs)
	#'oncogenic': (oncogenicMutNames, probs)}
	returnD = dict()
	hotspotNames = list(uniqueHotspotMutsAtQuadNuc['uniqueMutationIdentifier']) 
	oncogenicMutNames = list(uniqueOncogenicMutsAtQuadNuc['uniqueMutationIdentifier']) 

	returnD['hotspot'] = (hotspotNames + ['nonHotspot'], 
		[1.0/nTimesTrinucAppearsInGene for i in hotspotNames] + [1.0*(nTimesTrinucAppearsInGene - len(hotspotNames))/nTimesTrinucAppearsInGene])

	returnD['oncogenic'] = (oncogenicMutNames + ['nonOncogenic'], 
		[1.0/nTimesTrinucAppearsInGene for i in oncogenicMutNames] + [1.0*(nTimesTrinucAppearsInGene - len(oncogenicMutNames))/nTimesTrinucAppearsInGene])
	
	return returnD"""


	

#OLD VERSION
#We call this function with a big code chunk that looks something like this: quadNucDict, quadNucDictAdjusted, geneOncogenicOrHotspotProbDict, geneMutProbInfo = mutation_modeling_util.initiate_models(mafHere, mafBackground, geneDistributionsDf)	
"""def initiate_gene_mut_mapping(mafWithOncogenicAnnotations, geneList): #TO improve runtime we precalculate a dict mapping 
															#gene: maf of all mutations and maf of oncogenic mutations
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	geneToMutDataDict = dict()
	for gene in geneList:
		geneMuts = mafWithOncogenicAnnotations[mafWithOncogenicAnnotations['Hugo_Symbol'] == gene]
		oncogenicGeneMuts = geneMuts[geneMuts['oncogenic'].isin(oncogenicMutColNames)]
		hotspotGeneMuts = geneMuts[geneMuts['is-a-hotspot'] == 'Y']
		geneToMutDataDict[gene] = (oncogenicGeneMuts, hotspotGeneMuts)
	return geneToMutDataDict"""


"""
#
def process_background_maf(maf):
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	nonOncogenicMuts = maf[~maf['oncogenic'].isin(oncogenicMutColNames)]
	nonOncogenicSnps = nonOncogenicMuts[nonOncogenicMuts['Variant_Type'] == 'SNP']
	return nonOncogenicSnps


#geneMutProbInfo, geneRateProbInfo = mutation_modeling_util.init_gene_mutability_profile(mafBackground, geneDistributionsDf)

#Maf: a maf with as many mutations as possible (silent plus etc)
def init_gene_mutability_profile(maf, geneMutDistInfo):
    
    maf = process_background_maf(maf)
    geneLengthD = dict(zip(geneMutDistInfo['Hugo_Symbol'], geneMutDistInfo['cds_length']))
    geneNMutD = dict()
    geneRateD = dict()
    
    impactGenes = ['INPP4A', 'GNA11', 'MEF2B', 'FGF4', 'KEAP1', 'ESR1', 'MAPK1', 'NFKBIA', 'STAG2', 'NTHL1', 'TSC2', 'MGA', 'AGO2', 'CD79A', 'PIK3R2', 'KMT2B', 'ERF', 'HGF', 'PHOX2B', 'CCND1', 'GLI1', 'CDKN1B', 'RAB35', 'MLH1', 'BCL6', 'MSH2', 'MSH6', 'TNFAIP3', 'DUSP4', 'CXCR4', 'FLT3', 'INHBA', 'INHA', 'HIST1H3B', 'CDKN1A', 'SOX9', 'RRAS', 'TRAF2', 'RAC2', 'SMO', 'KNSTRN', 'MYOD1', 'FOXA1', 'RAF1', 'SESN2', 'LATS1', 'RARA', 'DNAJB1', 'H3F3B', 'KRAS', 'RRAS2', 'VHL', 'NOTCH2', 'PDGFRA', 'APC', 'MSI1', 'HNF1A', 'TBX3', 'CDK4', 'TMEM127', 'BRIP1', 'FGFR3', 'BARD1', 'CCND2', 'PALB2', 'CDH1', 'PDGFRB', 'FLT4', 'PIK3C3', 'SMAD2', 'RHEB', 'KMT2C', 'AXIN1', 'CREBBP', 'CCNE1', 'CDKN2C', 'PIK3R3', 'UPF1', 'MAP2K2', 'INPP4B', 'MAPK3', 'SMARCB1', 'EP300', 'EED', 'BRD4', 'NOTCH3', 'BIRC3', 'ACVR1', 'EPAS1', 'EPCAM', 'AKT3', 'KDR', 'PIK3CA', 'CTCF', 'CBL', 'CUL3', 'STAT3', 'DNMT3A', 'TP63', 'SDHA', 'MAP3K13', 'MSH3', 'RAD50', 'NBN', 'CDK6', 'PMS2', 'MAPKAP1', 'PIK3C2G', 'ERBB3', 'RB1', 'RAD51', 'IGF1R', 'ZFHX3', 'NCOR1', 'TP53', 'ERBB2', 'NUF2', 'EPHA5', 'PLK2', 'PIK3R1', 'RASA1', 'EGFR', 'NDUFB11', 'PRDM14', 'CDKN2B', 'NTRK2', 'NOTCH1', 'ATM', 'ARID5B', 'EIF4E', 'MYCN', 'FBXW7', 'CYSLTR2', 'FLT1', 'YAP1', 'MSI2', 'TCEB1', 'FLCN', 'ERCC3', 'CSF1R', 'GNAQ', 'PPARG', 'KIT', 'ERG', 'PREX2', 'BRAF', 'FANCC', 'PIK3CB', 'HIST1H2BD', 'HOXB13', 'U2AF1', 'FGFR4', 'DCUN1D1', 'STAT5B', 'SLX4', 'FGF19', 'REL', 'PRKCI', 'MST1R', 'NPM1', 'SOX17', 'RAD21', 'TSHR', 'INPPL1', 'TSC1', 'SPRED1', 'GREM1', 'RUNX1', 'ANKRD11', 'KMT2D', 'AXL', 'SDHAF2', 'CTLA4', 'INSR', 'IL7R', 'CDKN2A', 'IRS1', 'HIST1H3G', 'PPM1D', 'RPTOR', 'IGF1', 'AXIN2', 'MAP2K1', 'BCL2L1', 'ZRSR2', 'NUP93', 'BTK', 'EGFL7', 'TERT', 'HRAS', 'ERCC4', 'RPS6KB2', 'AURKA', 'YES1', 'CALR', 'GSK3B', 'PMAIP1', 'WHSC1L1', 'CD276', 'ABL1', 'FOXP1', 'ALOX12B', 'EZH2', 'POLE', 'FAM175A', 'PPP2R1A', 'SUZ12', 'MRE11A', 'EIF4A2', 'ARID1A', 'GTF2I', 'SOX2', 'PGR', 'SHQ1', 'TRAF7', 'STK11', 'CARM1', 'SMAD3', 'DNMT3B', 'CHEK2', 'HIST1H3I', 'IDH2', 'AMER1', 'FOXL2', 'SETD8', 'GRIN2A', 'IKZF1', 'HIST2H3D', 'PTCH1', 'PRKD1', 'SOCS1', 'WT1', 'BCL2', 'FGF3', 'RPS6KA4', 'ARID2', 'PDCD1', 'SF3B1', 'CENPA', 'LMO1', 'RAD51D', 'PNRC1', 'ASXL2', 'EPHA3', 'RAD51C', 'FIP1L1', 'MEN1', 'NF2', 'H3F3C', 'DNMT1', 'GATA2', 'SH2B3', 'PDPK1', 'JAK1', 'ERBB4', 'SMAD4', 'DICER1', 'HIST1H1C', 'CDC42', 'DROSHA', 'SMARCA4', 'TCF3', 'IDH1', 'STAT5A', 'ARID1B', 'GATA3', 'E2F3', 'SPOP', 'MALT1', 'AKT1', 'CTNNB1', 'ATR', 'PTPN11', 'MITF', 'PAK7', 'MAP2K4', 'TAP1', 'CRKL', 'NKX2-1', 'BLM', 'XIAP', 'RET', 'TNFRSF14', 'ERCC5', 'RAC1', 'PAK1', 'PTPRD', 'HIST1H3D', 'PPP4R2', 'PTPRS', 'RICTOR', 'HIST1H3A', 'BRCA1', 'NEGR1', 'PAX5', 'SRC', 'NF1', 'CASP8', 'FGFR2', 'RAD52', 'PRKAR1A', 'MAX', 'TGFBR2', 'PIK3CG', 'HIST1H3J', 'XRCC2', 'PLCG2', 'BABAM1', 'ELF3', 'SRSF2', 'HIST1H3E', 'WWTR1', 'NTRK3', 'TOP1', 'MTOR', 'CSF3R', 'SDCCAG8', 'FH', 'HIST3H3', 'PARP1', 'H3F3A', 'PARK2', 'IKBKE', 'MDM4', 'CDC73', 'RFWD2', 'IFNGR1', 'DDR2', 'SDHC', 'INSRR', 'RIT1', 'ROS1', 'FYN', 'MCL1', 'PRDM1', 'HIST1H3H', 'EPHA7', 'FAM46C', 'SHOC2', 'VTCN1', 'NRAS', 'SUFU', 'BCL10', 'PTP4A1', 'FUBP1', 'GNAS', 'SH2D1A', 'JUN', 'PTEN', 'RAD54L', 'NCOA3', 'BMPR1A', 'MUTYH', 'MPL', 'CCND3', 'RRAGC', 'STK40', 'PTPRT', 'ATRX', 'PIM1', 'PPP6C', 'TET1', 'MED12', 'DAXX', 'ID3', 'KLF4', 'AR', 'TAP2', 'TGFBR1', 'NOTCH4', 'STK19', 'KDM5C', 'SDHB', 'SDHD', 'ASXL1', 'SYK', 'SPEN', 'IRS2', 'MDC1', 'GATA1', 'HLA-A', 'ARAF', 'PIK3CD', 'ERRFI1', 'RBM10', 'DIS3', 'KDM6A', 'MYC', 'BCOR', 'FOXO1', 'EIF1AX', 'TET2', 'TEK', 'BRCA2', 'GPS2', 'NKX3-1', 'IRF4', 'CDK8', 'CRLF2', 'CD274', 'JAK2', 'TP53BP1', 'LATS2', 'WHSC1', 'SMYD3', 'ALK', 'FANCA', 'ERCC2', 'AKT2', 'CD79B', 'BCL2L11', 'PBRM1', 'SMARCD1', 'MYD88', 'ETV6', 'CARD11', 'NFE2L2', 'MYCL', 'PDCD1LG2', 'MET', 'EPHB1', 'CYLD', 'TMPRSS2', 'DOT1L', 'AC129492.6', 'MAP3K1', 'KDM5A', 'XPO1', 'SOS1', 'OBSL1', 'ETV1', 'FAM58A', 'ICOSLG', 'RNF43', 'SETD2', 'HLA-B', 'CBFB', 'RHOA', 'RECQL', 'IL10', 'FGFR1', 'RECQL4', 'EZH1', 'CHEK1', 'IGF2', 'SESN1', 'CSDE1', 'NSD1', 'POLD1', 'PMS1', 'FAT1', 'HIST1H3F', 'CDK12', 'BBC3', 'RP11-211G3.3', 'MST1', 'JAK3', 'BAP1', 'MDM2', 'RYBP', 'RXRA', 'RAD51B', 'CEBPA', 'RTEL1', 'LYN', 'VEGFA', 'NTRK1', 'KMT2A', 'SESN3', 'HIST1H3C', 'TIMM8B', 'TCF7L2', 'B2M', 'CIC', 'AURKB', 'MFSD11']
    for gene in impactGenes:
        nmuts = maf[maf['Hugo_Symbol'] == gene].shape[0]
        nmutsPerCDSLen = 1.0*nmuts/geneLengthD[gene]
        geneNMutD[gene] = nmuts
        geneRateD[gene] = nmutsPerCDSLen
        
    normedGeneProbCounter = normalize_counter(Counter(geneNMutD), mode='DontRound')
    genesListGeneMutProb = dict(normedGeneProbCounter).keys()
    probsListGeneMutProb = dict(normedGeneProbCounter).values()
    geneMutProbInfo = (genesListGeneMutProb, probsListGeneMutProb)    
    
    normedGeneRateProbCounter = normalize_counter(Counter(geneRateD), mode='DontRound')
    genesListGeneRateProb = dict(normedGeneRateProbCounter).keys()
    probsListGeneRateProb = dict(normedGeneRateProbCounter).values()
    geneRateProbInfo = (genesListGeneRateProb, probsListGeneRateProb)
    
    return geneMutProbInfo, geneRateProbInfo

#BIG initiation functions that intitate all the background probability information for all muts in impact


#a funciton to adjust the gene probabilities that we have for prob gene G mutated| given quadNuc Q mutated 
def norm_gene_probs(genes, probs, geneRateProbInfo):
	geneRateProbInfoAsDict = dict(zip(geneRateProbInfo[0], geneRateProbInfo[1]))
	cntr = 0
	for gene in genes:
		probs[cntr] *= geneRateProbInfoAsDict[gene]
		cntr += 1 
	probsSum = sum(probs)
	probsNormed = [i/probsSum for i in probs]
	return (genes, probsNormed)


#Function carved out to do a subset of the initiation work for the MAF dict etc
def do_initiation(mafBackground, mafWithOncogenicAnnotations, geneDistributionsDf):

	#create a unique identifer--this is useful later when we simulate picking individual hotspots
	mafWithOncogenicAnnotations['uniqueMutationIdentifier'] = mafWithOncogenicAnnotations.apply(lambda row: str(row['Hugo_Symbol'])+ '_' + str(row['HGVSp_Short']) + '_@' + str(row['Start_Position']), axis=1)

	#get the number of quadnucs in a gene
	quadNucColumns = geneDistributionsDf.columns.values[3:]
	geneDistributionsDf['NQuadNucInGene'] = geneDistributionsDf.apply(lambda row: sum(row[quadNucColumns]),axis=1)

	genes = list(geneDistributionsDf['Hugo_Symbol'])
	print 'Doing initation of mapping from gene to oncogenic mutations' #we precompute this and store it for runtime
	geneMutMappingDict = initiate_gene_mut_mapping(mafWithOncogenicAnnotations, genes)

	geneMutProbInfo, geneRateProbInfo = init_gene_mutability_profile(mafBackground, geneDistributionsDf)
	return geneDistributionsDf, quadNucColumns, genes, geneMutMappingDict, geneMutProbInfo, geneRateProbInfo, mafWithOncogenicAnnotations
"""

"""
#initiates all necessary infomration for a variety of models
#Returns:
    #quadNucDict: a dict mapping quadNucs to p(gene G mutated|quadNuc Q)
    #quadNucDictAdjusted: a dict mapping quadnucs to p(gene G mutated|quadNuc Q & normalized by gene mut freqs)
    #geneMutProbInfo: a list of genes and their probabilities of being mutated (determined from pan impact)
    #geneOncogenicOrHotspotProbDict: a dict mappning p(oncogenic mutation O| gene G and quadnuc Q)
    
def initiate_models(
	mafWithOncogenicAnnotations, #a maf including oncogenic/hotspot annotations for mutations
	mafForBackground, #a large MAF including silent mutations 
    geneDistributionsDf):
    
    #CHECKS FOR NECESSARY COLUMN NAMES
	if 'quadNuc' not in mafWithOncogenicAnnotations.columns.values:
		print 'error quadnuc column not in maf, fix it then feed this function again'
		return 0

	#some of the initial work of this function is decomposed into a "do_initiation" function
	df, quadNucColumns, genes, geneMutMappingDict, geneMutProbInfo, geneRateProbInfo, mafWithOncogenicAnnotations = do_initiation(mafForBackground, mafWithOncogenicAnnotations, geneDistributionsDf)

	#THESE ARE THE THINGS WE RETURN
	quadNucDict = dict()
	quadNucDictAdjusted = dict()
	geneOncogenicOrHotspotProbDict = dict()
	geneMutProbInfo, geneRateProbInfo = init_gene_mutability_profile(mafWithOncogenicAnnotations, df)

	#NOW WE DO REASONING ON THE MAF
	for qCol in quadNucColumns:
		print 'Performing analyses for quad nuc: ', qCol
		counts = df[qCol]
		percentages = [1.0*c/sum(list(counts)) for c in list(counts)]
		quadNucDict[qCol] = (genes, percentages) #store a tuple of names (genes) and counts 

		#need to use [:] to copy the list here otherwise shit gets fucked up 
		quadNucDictAdjusted[qCol] = norm_gene_probs(genes[:], percentages[:], geneRateProbInfo) #Create a dictionary called the quadNucDict adjusted which is like quad nuc dict but with an adjustment for gene mutability

		quadNucCountsDict = dict(zip(list(df['Hugo_Symbol']), list(df[qCol]))) #do this BC local variable counts gets changed
		
		#Now iterate through all the possible genes and calculate the probability P Oncogenic mutatiobn in Gene|mutation at motif
		print 'doing gene analyses'
		for gene in genes:
			probs = do_oncogenic_mut_initiation(qCol, geneMutMappingDict[gene], quadNucCountsDict[gene])
			geneOncogenicOrHotspotProbDict[gene + '_' + qCol] = probs
            
              #WE ALSO NEED P(mutGene&MutOncogenic)
            
	return quadNucDict, quadNucDictAdjusted, geneOncogenicOrHotspotProbDict, geneMutProbInfo

########
#Call this function with: mutation_modeling_util.write_model_information(quadNucDict=quadNucDict, quadNucDictAdjusted=quadNucDictAdjusted, geneOncogenicOrHotspotProbDict=geneOncogenicOrHotspotProbDict, geneMutProbInfo=geneMutProbInfo, write_dir='~/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mutationSimultationData',)
#Writes out the modeling information for later 
def write_model_information(quadNucDict=None, quadNucDictFilename='quadNucProbsFile.tsv',
	quadNucDictAdjusted=None, quadNucDictAdjustedFilename='quadNucProbsAdjFile.tsv',
	geneOncogenicOrHotspotProbDict=None, geneOncogenicOrHotspotProbDictFilename='geneOncogenicOrHotspotProbDictFile.tsv',
	geneMutProbInfo=None, geneMutProbInfoFilename='geneMutProbInfoFile.tsv',
	write_dir='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mutationSimultationData', prefix='IMPACT_'):
	
	def convert_quadNucDict_to_df(qNucDict):
		listOfDicts = []
		for quadnuc in qNucDict:
			genes = quadNucDict[quadnuc][0]
			values = quadNucDict[quadnuc][1]
			i = 0
			for gene in genes:
				listOfDicts.append({'quadNuc': quadnuc, 'gene': gene, 'prob': values[i]})
				i += 1
		df = pd.DataFrame(listOfDicts)
		return df

	def convert_geneOncogenicOrHotspotProbDict_to_df(geneOncogenicOrHotspotProbDict):
		listOfDicts = []
		for key, value in geneOncogenicOrHotspotProbDict.items():
			localD = {}
			hotspotNames = ','.join([str(i) for i in value['hotspot'][0]])
			oncogenicNames = ','.join([str(i) for i in value['oncogenic'][0]])
			hotspotProbs = ','.join([str(i) for i in value['hotspot'][1]])
			oncogenicProbs = ','.join([str(i) for i in value['oncogenic'][1]])
			localD['hotspotNames'] = hotspotNames
			localD['oncogenicNames'] = oncogenicNames
			localD['hotspotProbs'] = hotspotProbs
			localD['oncogenicProbs'] = oncogenicProbs
			listOfDicts.append(localD)
		df = pd.DataFrame(listOfDicts)
		return df

	def convert_gene_mut_prob_info_to_df(geneMutProbInf):
		listOfDicts = []
		for i in range(len(geneMutProbInf[0])):
			listOfDicts.append({'gene': geneMutProbInf[0][i], 'prob': geneMutProbInf[1][i]})
		return pd.DataFrame(listOfDicts)


	quadNucDf = convert_quadNucDict_to_df(quadNucDict)
	quadNucDfAdj = convert_quadNucDict_to_df(quadNucDictAdjusted)
	geneOncogenicOrHotspotProbDf = convert_geneOncogenicOrHotspotProbDict_to_df(geneOncogenicOrHotspotProbDict)
	geneMutProbDf = convert_gene_mut_prob_info_to_df(geneMutProbInfo)

	quadNucDf.to_csv(os.path.join(write_dir, prefix + quadNucDictFilename), sep='\t', index=False)
	quadNucDfAdj.to_csv(os.path.join(write_dir, prefix + quadNucDictAdjustedFilename), sep='\t', index=False)
	geneOncogenicOrHotspotProbDf.to_csv(os.path.join(write_dir, prefix + geneOncogenicOrHotspotProbDictFilename), sep='\t', index=False)
	geneMutProbDf.to_csv(os.path.join(write_dir, prefix + geneMutProbInfoFilename), sep='\t', index=False)

	print 'written to files in the directory ', write_dir
"""

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--fastaPath', help='Path for full FASTA file', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/b37.fasta')
	parser.add_argument('--bedFilePath', help='Path of full Bed file', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/genelist.with_aa.interval_list')   
	parser.add_argument('--scratchBedPath', help='pathForTemporaryBedFile', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/tempScriptFiles/tempBed')
	
	args = parser.parse_args()
	
	enumerate_all_possible_snps(sequenceDfPath = args.bedFilePath,
		genes = ['TP53'],
		scratchBedFilePath = args.scratchBedPath)

if __name__ == '__main__':
    main()






