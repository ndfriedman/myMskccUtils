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
	parser.add_argument('inputMaf', default='')
	args = parser.parse_args()

	if args.mode == 'subsetByCancerType':
		cases = get_cases_for_cancer_type('Bladder_Cancer')
		print cases

if __name__ == '__main__':
    main()















