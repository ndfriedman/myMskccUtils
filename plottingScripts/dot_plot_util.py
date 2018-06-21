#written by Noah Friedman 
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
import matplotlib.cm as cm

#--------------------

#TAKEN From: http://schubert.atmos.colostate.edu/~cslocum/custom_cmap.html

#-------------------------

def dot_plot_with_words(l): #stub code to make a dot plot where you graph words instead of dots
	fig, ax = plt.subplots()
	for entry in l:
		sense, specif, name = entry
		ax.text(sense, specif, name, ha='center', size=20)
	plt.savefig('noahTest')
	sys.exit()
	plt.clf()

def scatter_plot_with_legend(l, title):

	def make_scatter_plot_legend(legendInfo): #make a legend for a scatter plot
		scatInfo = [x[1] for x in legendInfo]
		names = [x[0] for x in legendInfo]
		plt.legend(scatInfo,
           names,
           scatterpoints=1,
           loc='lower left',
           ncol=3,
           fontsize=8)

	colors = cm.rainbow(np.linspace(0, 1, len(l)))
	cntr = 0
	legendL = []
	for entry in l:
		sense, specif, name = entry
		curColor = colors[cntr]
		scat = plt.scatter(sense, specif, s=100, c=curColor)
		cntr += 1
		legendL.append((name, scat))

	make_scatter_plot_legend(legendL)
	plt.xlabel('Sensitivity')
	plt.ylabel('Specificity')
	plt.xlim([0,1])
	#plt.ylim([.8,1])
	plt.title(title)

	filename = 'noahTest'
	print 'saving figure to ', filename
	plt.savefig(filename)
	sys.exit()

def scatter_plot_line(vals, title):
	#sys.exit()
	senses = [x[0] for x in vals]
	specifs = [x[1] for x in vals]
	colors = [x[2] for x in vals]

	plt.scatter(senses, specifs, c=colors, cmap=plt.cm.get_cmap('RdBu'))
	cbar = plt.colorbar()
	cbar.set_label('1 minus confidence')

	plt.title(title)
	plt.xlim([0,1.1])
	plt.ylim([0,1.1])
	plt.xlabel('Specificity')
	plt.ylabel('Sensitivity')
	#plt.plot(senses, specifs)


	filename = 'signatureConfidenceCapture' + title
	print 'saving figure to ', filename
	plt.savefig(filename)
	plt.clf()

def overlaid_scatterplot_lines(vals, title, lineLabel=None, doSave=False):
	senses = [x[0] for x in vals]
	specifs = [x[1] for x in vals]
	plt.plot(senses, specifs, label=lineLabel)
	plt.xlim(0, 1.1)
	plt.ylim(0, 1.1)
	plt.xlabel('Specificity')
	plt.ylabel('Sensitivity')
	ax = plt.gca()
	handles, labels = ax.get_legend_handles_labels()
	if doSave:
		#create the legend as well
		plt.legend()
		plt.title('Sensitivity/Specificity curves Stratified by Amount of Signature 3 Present in Exome')
		#plt.legend()

		filename = title + '_linePlot'
		print 'saving figure to ', filename
		plt.savefig(filename)
	return 0

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')



	args = parser.parse_args()

if __name__ == '__main__':
    main()















