#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')


def shuffle_arr(arr):
    nameCol = arr[:,0]
    valueCol = arr[:,1]
    np.random.shuffle(nameCol)
    return np.column_stack((nameCol, valueCol))

def convert_df_to_np_array(df):
    dfReduced = df[['Tumor_Sample_Barcode', 'shuffleColumn']]
    #npArr = dfReduced.to_numpy()
    npArr = dfReduced.values
    return npArr

def summarize_shuffle_results_for_double_mut_analysis_gene_dist(shuffledArr, genes):
	zippedList = zip(shuffledArr[:,0], shuffledArr[:,1])

	df = pd.DataFrame([])
	df['Tumor_Sample_Barcode'] = shuffledArr[:,0]
	df['Hugo_Symbol'] = shuffledArr[:,1]
	df['geneCase'] = df.apply(lambda row: row['Tumor_Sample_Barcode'] + '_' + row['Hugo_Symbol'], axis=1)
	return dict(df['geneCase'].value_counts())

	#cntr = Counter(zippedList)
	#return Counter(cntr.values())

def do_n_iterations_of_shuffling(maf, n, genes, mode='geneDist'):
    npArr = convert_df_to_np_array(maf)
    listOfDicts = []
    for i in range(n):
        if i%10 == 0: print 'on iteration ', i
        shuffledArr = shuffle_arr(npArr)
        if mode == 'geneDist':
            shuffleResults = summarize_shuffle_results_gene_dist(shuffledArr, genes)
            listOfDicts.append(shuffleResults)
        else:
            d = summarize_shuffle_results(shuffledArr, genes)
        
        #listOfDicts.append(d)
    return pd.DataFrame(listOfDicts)

#calculates for each observation (all genes, or n mutations in a gene the number of times the permutation test was greater)
def calculate_n_tests_that_differed(testSummaryDf, actualD):
	d = {}
	for val in testSummaryDf.columns.values:
		d[val] = testSummaryDf[testSummaryDf[val] >= actualD[val]].shape[0]
	return d


def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()

if __name__ == '__main__':
    main()















