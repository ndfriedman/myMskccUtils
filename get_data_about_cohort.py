#written by Noah Friedman 
#little util that tells us stuff like how many cases have enough mutations etc
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	maf = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impactMafs/data_mutations_unfiltered_mafAnno.maf'
	df = pd.read_table(maf, skiprows=[0])

	d = {}
	for index, row in df.iterrows():
		tCode = row['Tumor_Sample_Barcode']
		if tCode in d:
			d[tCode] = d[tCode] + 1
		else:
			d[tCode] = 1

	cntr = 0
	tTypes = []
	for key, value in d.items():
		if value > 10:
			cntr += 1
			tTypes.append(key)

	f = open('tumorIdsWithEnoughMuts', 'w')
	for tCode in tTypes:
		f.write(tCode + '\n')

	args = parser.parse_args()

if __name__ == '__main__':
    main()















