#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

def move_last_df_col_to_front(df):
	cols = df.columns.tolist()
	cols = cols[-1:] + cols[:-1]
	df = df[cols]
	return df

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--csvPath', help='stub for parser argument', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mutationSigFiles/sigProfiler_exome_SBS_signatures.csv')

	args = parser.parse_args()

	df = pd.read_csv(args.csvPath)
	df['quadNuc'] = df.apply(lambda row: row['SubType'][0] + row['Type'].replace('>', '') + row['SubType'][2], axis=1)
	df = df.drop(['Type', 'SubType'], axis=1)
	df = move_last_df_col_to_front(df)
	df = df.transpose()
	
	print df

	df.columns = df.iloc[0] #make the first row the columns for the df
	df['signature'] = df.index
	df = move_last_df_col_to_front(df)
	

	df = df.drop(df.index[[0]]) #drop the first row which is now redundant
	df.to_csv('newSignatures.txt', sep='\t', index=False)

if __name__ == '__main__':
    main()















