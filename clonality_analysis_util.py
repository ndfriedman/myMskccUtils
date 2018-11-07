#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

from myUtils import data_compacting_and_cleaning_util

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

def init_for_facets_analysis(): #does a suite of initiation steps for the facets df
	return 0

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()

if __name__ == '__main__':
    main()















