#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np
from plotnine import *
import matplotlib.pyplot as plt


from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')


def test_func(df):


	print 'zeniba'
	#plotnine.options.figure_size = (6.4, 4.8)
	"""
	x = (ggplot(df) +
 	 geom_bar(
 	 	aes(x='hrd_cancer',y='brca_signature'), stat='identity'
 	 	)
 	 	+
 	 	facet_wrap('~biallelic_class') #the facet wrap is the ordering of subplots
 	 )"""

	

	x = (ggplot(df) +
 	 	aes(x='hrd_cancer',fill='cancer_type') +
	geom_bar(position = "fill")
 	 	+
 	 	facet_wrap('~biallelic_class') #the facet wrap is the ordering of subplots
 	 )
	
	display(x)



def plot_df_by_individual_patient(df, manualXLabel):
	"""x = (ggplot(df) +
 	 	aes(x='dmp_sample',fill='brca_signature') +
	geom_bar(position = "stack")
 	+
 	facet_wrap('~biallelic_class') #the facet wrap is the ordering of subplots
 	+
	theme(axis_text_x=element_text(rotation=90, hjust=1)) #rotate x labels so they are less obnoxious
 	)"""

 	x = (ggplot(df) +
 	 	aes(fill='labelSig', y='mean', x='sigLexicalOrder') +
	geom_bar(stat="identity", position="fill")
 	+ 
 	#scale_x_discrete(labels= manualXLabel)
 	#+
 	facet_wrap('~biallelic_class', scales='free_x') #the facet wrap is the ordering of subplots
 	+
	theme(axis_text_x=element_text(rotation=90, hjust=1)) #rotate x labels so they are less obnoxious
 	)

	
	display(x)


def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()

if __name__ == '__main__':
    main()















