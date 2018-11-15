#written by Noah Friedman 
#a script containing functions used for analysis task
import sys
import argparse
import os
import pandas as pd
import numpy as np

import scipy.stats

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

#KAPLAN MEIER ANALYSIS

def test_significance(df1, df2):
	results = logrank_test(df1['os_years'], df2['os_years'], df1['CENSOR'], df2['CENSOR'], alpha=.99)
	return results.p_value
	#results.print_summary()

def make_kaplan_meier_plots(dfs):
	def correct_censorship(dfs): #creates a sensorship column which invertys the pt vital status colum (0 is alive, 1 is dead)
		for i in range(len(dfs)):
			df = dfs[i]
			df['CENSOR'] = df['PT_VITAL_STATUS'].apply(lambda x: 1 if x == 0 else 0)
			dfs[i] = df
		return dfs

	dfs = correct_censorship(dfs)
	kmf = KaplanMeierFitter()
	kmf.fit(dfs[0]['os_years'], dfs[0]['CENSOR'], label='cond')
	ax = kmf.plot()
	kmf.fit(dfs[1]['os_years'], dfs[1]['CENSOR'], label='nonCond')
	ax = kmf.plot(ax=ax)
	fig = ax.get_figure()
	#fig.savefig('testFig.pdf')

	print test_significance(dfs[0], dfs[1])


#INFORMATION FUNCTIONS ##################################################

#little utility for normalizing a counter
def normalize_counter(cntrObj, nDigitsRound=2):
	total = sum(cntrObj.values(), 0.0)
	for key in cntrObj:
		cntrObj[key] /= total
		cntrObj[key] = round(cntrObj[key], nDigitsRound)
	return cntrObj

#returns the mean value of a column specified by colname
def get_mean_of_df_col(df, colname, idColumn = 'Tumor_Sample_Barcode'):
	df = df.drop_duplicates(subset=[idColumn])
	return np.nanmean(np.asarray(list(df[colname])))

#returns the median of a df column
def get_median_of_df_col(df, colname, idColumn = 'Tumor_Sample_Barcode'):
	df = df.drop_duplicates(subset=[idColumn])
	return np.nanmedian(np.asarray(list(df[colname])))


def get_n_cases_with_mutation_present(maf, gene, idCol='Tumor_Sample_Barcode'):
	n = len(set(maf[maf['Hugo_Symbol'] == 'TP53'][idCol]))
	print n, 1.0*n/len(set(maf['Tumor_Sample_Barcode']))


def get_age_information(ageInformationPath = '/ifs/work/taylorlab/friedman/msk-impact/msk-impact/darwin/darwin_age.txt'):
	ageDf = pd.read_table(ageInformationPath)
	d = dict(zip(ageDf['PATIENT_ID'], ageDf['AGE']))
	return d

def get_cancer_type_information(cancerTypeDfPath = '/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt', mode='pid'):
	cancerTypeDf = pd.read_table(cancerTypeDfPath)
	if mode == 'pid':
		cancerTypeDf['#Sample Identifier'] = cancerTypeDf['#Sample Identifier'].apply(lambda x: x[:9])
	d = dict(zip(cancerTypeDf['#Sample Identifier'], cancerTypeDf['Cancer Type']))
	return d

#TODO MAKE THIS MORE SUSTAINABLE
def get_gene_length_info(bedFilePath = '/ifs/res/pwg/data/gencode/gencode.v19.all_gene_bounds.bed'):
	impactGenes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	bedDf = pd.read_table(bedFilePath)	
	bedDf = bedDf[bedDf['OR4F5'].isin(impactGenes)]
	bedDf['geneLength'] = bedDf.apply(lambda row: row['70008'] - row['69090'], axis=1)
	return dict(zip(bedDf['OR4F5'], bedDf['geneLength']))

#SIGNATURE SPECIFIC ANALYSIS UTILS #############################################################

#util to give the top N most epxressed signatures:
def get_n_top_signatures(row, n=2):
	#signatureCols = list(row.columns.values)
	row = row[['mean_' + str(i) for i in range(1,31)]]
	l = list(row)
	return list(reversed([str(i + 1) + ':' + str(l[i]) for i in np.argsort(l)[-n:]])) #I plus one to take into account the signatures ordering

#enumerates a set of all 96 possible trinucelotides/changes
def get_all_possible_quadNucs():
	allSpectra = []
	for firstLetter in ['A','C','G','T']:
		for change in ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']:
			for lastLetter in ['A','C','G','T']:
				allSpectra.append(firstLetter+change+lastLetter)
	return(set(allSpectra))

#STATISTICAL TESTS ####################################################

def test_significance(df1, df2):
	results = logrank_test(df1['os_years'], df2['os_years'], df1['CENSOR'], df2['CENSOR'], alpha=.99)
	return results.p_value

#does the mann whitney u test
def do_mann_whitney_test(dist1, dist2):
	return scipy.stats.mannwhitneyu(dist1, dist2).pvalue

#FISHER TEST PLEASE!!
def do_p_val_tests_with_comps(dist1, dist2, comparissonMessage=None):
	print comparissonMessage
	print 'Ns for each distribution: ', len(list(dist1)), len(list(dist2))
	print 'MEANS of distributions: ', np.nanmean(np.asarray(list(dist1))), np.nanmean(np.asarray(list(dist2)))
	print 'MEDIANS of distributions: ', np.nanmedian(np.asarray(list(dist1))), np.nanmedian(np.asarray(list(dist2)))
	print 'P val: ', do_mann_whitney_test(dist1, dist2)

#calcualtes a correlation and pearson p val for two columns of the same df
def do_pearson_correlation_one_df(df, col1, col2):
	df = df[np.isfinite(df[col1])]
	df = df[np.isfinite(df[col2])]
	c1 = np.asarray(df[col1])
	c2 = np.asarray(df[col2])
	#print c1, c2
	return scipy.stats.pearsonr(c1,c2)

#util function to pretty print tumor sample barcodes for cbioportal
def print_for_cbio_portal(s):
	for v in s:
		print v


def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()


if __name__ == '__main__':
    main()















