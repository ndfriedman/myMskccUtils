#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np
import math

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'

import imp
analysis_utils = imp.load_source('analysis_utils', '/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/analysis_utils.py')

def create_facets_dict_key(row):
	return row['Tumor_Sample_Barcode'] + '_' + row['idCol']

def create_facets_clonality_dict(facetsDf):
	#todo make it incorportate patient info too
  	facetsDf['idCol'] = facetsDf.apply(lambda row: str(row['Chromosome']) + '_' + str(row['Start_Position']), axis=1)
  	facetsDf = data_compacting_and_cleaning_util.create_expected_mut_copies_col(facetsDf)
  	facetsDf['naiveClonalStatus'] = facetsDf.apply(lambda row: row['VAF']/row['purity'], axis=1)
  	d = dict()
  	for index, row in facetsDf.iterrows():
  		d[create_facets_dict_key(row)] = (row['naiveClonalStatus'], row['ccf_Mcopies_upper'])
  	return d

def analyze_clonality_across_muts_in_pole_cases(mutsDf, facetsDict, writeDir = ''):
	cases = set(mutsDf['Tumor_Sample_Barcode'])
	for case in cases:
		caseMuts = mutsDf[mutsDf['Tumor_Sample_Barcode'] == case]
		data = []
		if caseMuts.shape[0] > 0:
			for index, row in caseMuts.iterrows():
				facetsDictKey = create_facets_dict_key(row)
				if facetsDictKey in facetsDict:
					ccf, ccf_Mcopies = facetsDict[facetsDictKey]
					if not np.isnan(float(ccf_Mcopies)):
						data.append(ccf_Mcopies)
			if len(data) > 0:
				figureTitle = case + '; NMuts: ' + str(caseMuts.shape[0])
				histogram_util.plot_simple_normed_histogram(data, case, writeDir, figureTitle)
			print '_________'


#a heuristic for testing if a mutation is a double hit of the same location
#assumes that a mutation that occurs in a balanced region with 
#TODO make sure that we dont miscall variants in a small founder clone as doubles

#TODO CHANGE THIS TO CALL DOUBLE MUTATIONS BASED ON SAMPLE PURITY
def is_mut_double_hit(row, 
	flatGenome, #BOOLEAN to tell us if facets said the case's genome is flat
 doubleFactor=2): #factor by which a mutations vaf needs to be bigger than the median to be considered double.  

	if row['Chromosome'] == 'X': return False #DONT CALL A MUTATION ON X AS DOUBLE BY MISTAKE

	if (math.isnan(row['ccf_Mcopies']) or math.isnan(row['ccf_1copy'])) and not flatGenome: #Null values for ccf only matter if the genome isnt flat
		return False
	else:
		if not flatGenome and row['tcn'] - row['lcn'] != row['lcn']: #if the region isnt balanced we wont call it a double hit
			return False
		else:
			if row['t_var_freq'] >= 1.75*row['medianClonalVaf']:
				return True

#MARK mutations that occur at a balanced region of the genome 
#use this for simplified analyses of balanced regions
def mark_mutation_is_balanced(row):
	if row['Chromosome'] == 'X': return False #we dont include X mutations cause its complicated
	
	#if flatGenome: return True #cases with purity = NA, flat genomes are considered to be flat everywhere
	#if (math.isnan(row['ccf_Mcopies']) or math.isnan(row['ccf_1copy'])) #dont analyze regions where Mcopies or 1copy are null (flat genome cases have already been called true)
	#	return False

	if row['tcn'] - row['lcn'] != row['lcn']: #this means the region isnt balanced
		return False
	else: #this means it is
		return True

def mark_maf_with_ccf_for_flat_genomes(maf, gene):
	
	maf = maf[maf['oncogenic'].notnull()]
	print set(maf['oncogenic'])

	for case in set(maf['Tumor_Sample_Barcode']):
		caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
		caseMafPTEN = caseMaf[caseMaf['Hugo_Symbol'] == gene]
		print caseMafPTEN['ccf_Mcopies']


def calculate_delta_vaf_across_mutation_pairs(maf, #MAF to analyze
	doDoubleHitCorrection=True, #does a correction to count all cases where nmut copies = 2 as two mutations at the same spot
	genes = None #genes to look at
	):
	

	#FIRST AND FOREMOST, its only relevant to do this analysis at diploid parts of the genome
	maf = maf[maf['tcn'] == 2]

	dDeltaVaf = {}
	if genes == None: 
		genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	cntr = 0
	nCases = len(set(maf['Tumor_Sample_Barcode']))
	for case in set(maf['Tumor_Sample_Barcode']):
		#COUNTING
		print 'analyzing case number ', cntr, ' out of ', nCases
		cntr += 1

		#SET SOME VARIABLES WE WILL NEED FOR THIS ANALYSIS
		caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
		maxTVarFreq = max(caseMaf['t_var_freq'])
		medianVaf = np.nanmedian(caseMaf['t_var_freq'])
		meanVaf = np.nanmean(caseMaf['t_var_freq'])

		caseMaf['isDouble'] = None #make a dummy column, fill it out if the do double correction is ons
		if doDoubleHitCorrection:
			flatGenome = False
			if caseMaf[caseMaf['ccf_Mcopies'].notnull()].shape[0] == 0: #CATCH ALL THE PURITY=NA flat genomes cases
				flatGenome = True

			caseMaf['isDouble'] = caseMaf.apply(lambda row: is_mut_double_hit(row, medianVaf, flatGenome, 1.75), axis=1)

		for gene in genes:
			geneMaf = caseMaf[caseMaf['Hugo_Symbol'] == gene]

			if geneMaf.shape[0] > 0:

				#DO NOT CONSIDER GENES WHERE THERE IS LOH or the gene is on CHR X
				if geneMaf[geneMaf['lcn'] == 0].shape[0] > 0 or geneMaf[geneMaf['Chromosome'] == 'X'].shape[0] > 0: 
					pass
				else:

					if geneMaf[geneMaf['isDouble'] == True].shape[0] > 0:
						pass
						#ALERT DO WE INCLUDE DOUBLE MUTS AT THE SAME LOCUS OR NOT???
						#diff = max(list(geneMaf['t_var_freq'])))
						#if gene in dDeltaVaf:
						#	dDeltaVaf[gene][0] = dDeltaVaf[gene][0] + [0] #append a zero for double mutations (note that this does not take into account that the double mutations may not occur at the same time!!?>)
						#	dDeltaVaf[gene][1] = dDeltaVaf[gene][1] + [(0.5*max(list(geneMaf['t_var_freq'])))/meanVaf] #take max vaf, divide it by two and compare the ratio to median vaf
						#else:
						#	dDeltaVaf[gene] = [[0],[(0.5*max(list(geneMaf['t_var_freq'])))/meanVaf]] #make it a list of lists
					###TODO add code so we dont do shit like when PTEN IS LOF'D already
					else:
						if geneMaf.shape[0] > 1:
							maxVafToMedianVafRatio = 1.0*max(list(geneMaf['t_var_freq']))/meanVaf

							differences = analysis_utils.calculate_all_pairwise_differences(np.array(list(geneMaf['t_var_freq'])))
							differencesNormed = [1.0*i/maxTVarFreq for i in differences]

							if gene in dDeltaVaf:
								dDeltaVaf[gene][0] = dDeltaVaf[gene][0] + [min(differencesNormed)] # we only take the closest pair if there is more than one pair
								dDeltaVaf[gene][1] = dDeltaVaf[gene][1] + [maxVafToMedianVafRatio]
							else:
								dDeltaVaf[gene] = [[min(differencesNormed)], [maxVafToMedianVafRatio]] # we only take the closest pair if there is more than one pair
	return dDeltaVaf


#TODO add functionality to do percentile for only the biggest mutation as well
def mark_mutations_by_vaf_mut_percentile(maf):

	#only do this at balanced regions of the genome
	maf['isBalanced'] = maf.apply(lambda row: mark_mutation_is_balanced(row), axis=1)
	maf = maf[maf['isBalanced'] == True]

	print maf.shape

	caseDict = {}
	for case in set(maf['Tumor_Sample_Barcode']):
		caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
		percentiles = dict(zip(caseMaf['varUuid'], caseMaf['t_var_freq'].rank(pct=True)))
		caseDict[case] = percentiles
	maf['caseVafPercentileRank'] = maf.apply(lambda row: caseDict[row['Tumor_Sample_Barcode']][row['varUuid']], axis=1)		
	return maf

def enumerate_biggest_and_second_biggest_vaf_for_mutation():
	return 0

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()

if __name__ == '__main__':
    main()















