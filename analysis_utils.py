#written by Noah Friedman 
#a script containing functions used for analysis task
import sys
import argparse
import os
import pandas as pd
import numpy as np

import scipy.stats


#INFORMATION FUNCTIONS ##################################################

#returns the mean value of a column specified by colname
def get_mean_of_df_col(df, colname, idColumn = 'PATIENT_ID'):
	df = df.drop_duplicates(subset=[idColumn])
	return np.nanmean(np.asarray(list(df[colname])))

#returns the median of a df column
def get_median_of_df_col(df, colname, idColumn = 'PATIENT_ID'):
	df = df.drop_duplicates(subset=[idColumn])
	return np.nanmedian(np.asarray(list(df[colname])))


#SIGNATURE SPECIFIC ANALYSIS UTILS #############################################################

#util to give the top N most epxressed signatures:
def get_n_top_signatures(row, n=2):
	l = list(row)
	return list(reversed([str(i + 1) + ':' + str(l[i]) for i in np.argsort(l)[-n:]])) #I plus one to take into account the signatures ordering


#STATISTICAL TESTS ####################################################

def test_significance(df1, df2):
	results = logrank_test(df1['os_years'], df2['os_years'], df1['CENSOR'], df2['CENSOR'], alpha=.99)
	return results.p_value

#does the mann whitney u test
def do_mann_whitney_test(dist1, dist2):
	return scipy.stats.mannwhitneyu(dist1, dist2).pvalue

#calcualtes a correlation and pearson p val for two columns of the same df
def do_pearson_correlation_one_df(df, col1, col2):
	df = df[np.isfinite(df[col1])]
	df = df[np.isfinite(df[col2])]
	c1 = np.asarray(df[col1])
	c2 = np.asarray(df[col2])
	#print c1, c2
	return scipy.stats.pearsonr(c1,c2)


def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()


if __name__ == '__main__':
    main()















