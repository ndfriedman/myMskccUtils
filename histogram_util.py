#written by Noah Friedman 
#a script with a lot of different ways to plot histograms basically
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

import matplotlib
matplotlib.use('Agg') #the order of this import crap is cirtical for avoiding a breaking error
import matplotlib.pyplot as plt


def plot_simple_normed_histogram(data, caseName, writeDir, title ='Title', titleSuffix='_clonality_graph.png'):
	indicies = np.array([.1*x for x in range(10)])
	indicies = np.arange(len(data))
	width = 1
	weights = np.ones_like(data)/float(len(data)) #code to make 
	plt.hist(data, weights=weights, bins=25)
	plt.ylim([0,1])
	plt.xlim([0,1])
	plt.title(title)
	saveFileName = caseName + titleSuffix
	savePath = os.path.join(writeDir, saveFileName)
	print 'saving graph to ', savePath
	plt.savefig(savePath)
	plt.clf()

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')
	args = parser.parse_args()

if __name__ == '__main__':
    main()















