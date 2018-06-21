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


#SIGNATURE SPECIFIC ANALYSIS UTILS #############################################################

#util to give the top N most epxressed signatures:
def get_n_top_signatures(row, n=2):
	#signatureCols = list(row.columns.values)
	row = row[['mean_' + str(i) for i in range(1,31)]]
	l = list(row)
	return list(reversed([str(i + 1) + ':' + str(l[i]) for i in np.argsort(l)[-n:]])) #I plus one to take into account the signatures ordering


#Conditionals for identifying mutations from specific signatures

#extract
def get_sig_17_pink_peak_muts(df, mode=True):
	if mode:
		return df[
		(df['Ref_Tri'] == 'CTT') #always the ref tri has to be 'CTT' (which implicitly includes 'ACC' as well) for sig 17 pink peak
		& (
		((df['Reference_Allele'] == 'T') & (df['Tumor_Seq_Allele2'] == 'G')) # given the ref tri context, they can be T > G
		| ((df['Reference_Allele'] == 'A') & (df['Tumor_Seq_Allele2'] == 'C')) #or they can be A>C
		)
		]

	else: 
		return df[~ #NOT!
		(
		(df['Ref_Tri'] == 'CTT') #always the ref tri has to be 'CTT' (which implicitly includes 'ACC' as well) for sig 17 pink peak
		& (
		((df['Reference_Allele'] == 'T') & (df['Tumor_Seq_Allele2'] == 'G')) # given the ref tri context, they can be T > G
		| ((df['Reference_Allele'] == 'A') & (df['Tumor_Seq_Allele2'] == 'C')) #or they can be A>C
		)
		)
		]

def get_tcga_sig_17_muts(df, mode=True):
	if mode:
		return df[
			(
			(df['Ref_Tri'] == 'CTT') 
			& (
			((df['Reference_Allele'] == 'T') & (df['Tumor_Seq_Allele2'] == 'G')) # given the ref tri context, they can be T > G
			| ((df['Reference_Allele'] == 'A') & (df['Tumor_Seq_Allele2'] == 'C')) #or they can be A>C
			)
			)

			|

			(
			(df['Ref_Tri'] == 'ATT') 
			& (
			((df['Reference_Allele'] == 'T') & (df['Tumor_Seq_Allele2'] == 'G')) # given the ref tri context, they can be T > G
			| ((df['Reference_Allele'] == 'A') & (df['Tumor_Seq_Allele2'] == 'C')) #or they can be A>C
			)
			)

			|

			(
			(df['Ref_Tri'] == 'GTT') 
			& (
			((df['Reference_Allele'] == 'T') & (df['Tumor_Seq_Allele2'] == 'G')) # given the ref tri context, they can be T > G
			| ((df['Reference_Allele'] == 'A') & (df['Tumor_Seq_Allele2'] == 'C')) #or they can be A>C
			)
			)

			|

			(
			(df['Ref_Tri'] == 'TTT') 
			& (
			((df['Reference_Allele'] == 'T') & (df['Tumor_Seq_Allele2'] == 'G')) # given the ref tri context, they can be T > G
			| ((df['Reference_Allele'] == 'A') & (df['Tumor_Seq_Allele2'] == 'C')) #or they can be A>C
			)
			)
		]

	else:
		return df[~ #NOT!!!
			(
			(
			(df['Ref_Tri'] == 'CTT') 
			& (
			((df['Reference_Allele'] == 'T') & (df['Tumor_Seq_Allele2'] == 'G')) # given the ref tri context, they can be T > G
			| ((df['Reference_Allele'] == 'A') & (df['Tumor_Seq_Allele2'] == 'C')) #or they can be A>C
			)
			)

			|

			(
			(df['Ref_Tri'] == 'ATT') 
			& (
			((df['Reference_Allele'] == 'T') & (df['Tumor_Seq_Allele2'] == 'G')) # given the ref tri context, they can be T > G
			| ((df['Reference_Allele'] == 'A') & (df['Tumor_Seq_Allele2'] == 'C')) #or they can be A>C
			)
			)

			|

			(
			(df['Ref_Tri'] == 'GTT') 
			& (
			((df['Reference_Allele'] == 'T') & (df['Tumor_Seq_Allele2'] == 'G')) # given the ref tri context, they can be T > G
			| ((df['Reference_Allele'] == 'A') & (df['Tumor_Seq_Allele2'] == 'C')) #or they can be A>C
			)
			)

			|

			(
			(df['Ref_Tri'] == 'TTT') 
			& (
			((df['Reference_Allele'] == 'T') & (df['Tumor_Seq_Allele2'] == 'G')) # given the ref tri context, they can be T > G
			| ((df['Reference_Allele'] == 'A') & (df['Tumor_Seq_Allele2'] == 'C')) #or they can be A>C
			)
			)
			)
		]		

def get_pole_mutations(df, mode=True): #note lazy way just gets the big peaks
	
	if mode:
	
		return df[
			(
			(
			(df['Ref_Tri'] == 'TCT') 
			& (
			((df['Reference_Allele'] == 'C') & (df['Tumor_Seq_Allele2'] == 'A')) 
			| ((df['Reference_Allele'] == 'G') & (df['Tumor_Seq_Allele2'] == 'T')) 
			)
			)

			|

			(
			(df['Ref_Tri'] == 'TCG') 
			& (
			((df['Reference_Allele'] == 'C') & (df['Tumor_Seq_Allele2'] == 'T')) 
			| ((df['Reference_Allele'] == 'G') & (df['Tumor_Seq_Allele2'] == 'A')) 
			)
			)
			)
			]

	else:
	
		return df[~
			(
			(
			(df['Ref_Tri'] == 'TCT') 
			& (
			((df['Reference_Allele'] == 'C') & (df['Tumor_Seq_Allele2'] == 'A')) 
			| ((df['Reference_Allele'] == 'G') & (df['Tumor_Seq_Allele2'] == 'T')) 
			)
			)

			|

			(
			(df['Ref_Tri'] == 'TCG') 
			& (
			((df['Reference_Allele'] == 'C') & (df['Tumor_Seq_Allele2'] == 'T')) 
			| ((df['Reference_Allele'] == 'G') & (df['Tumor_Seq_Allele2'] == 'A')) 
			)
			)
			)
			]



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















