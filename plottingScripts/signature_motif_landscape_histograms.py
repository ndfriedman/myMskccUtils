#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg') #the order of this import crap is cirtical for avoiding a breaking error
import matplotlib.pyplot as plt

import scipy
from scipy import spatial
import pylab
import scipy.cluster.hierarchy as sch

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

def get_gene_panel_genes(genePanelPath = '/ifs/work/taylorlab/friedman/clinicalData/msk-impact/msk-impact/gene_panels/impact468_gene_panel.txt'):
	f = open(genePanelPath)
	lines = f.readlines()
	geneSet = set(lines[3].strip('\n').split('\t'))
	return geneSet

#from https://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python
def dendogram_function(D, geneNames=None):
	
	def create_dendogram_ticks(axmat, ordering, geneNames):
		sortedGeneNames = [x for _,x in sorted(zip(ordering, geneNames))]
		axmat.set_xticks(range(len(sortedGeneNames)))
		axmat.set_xticklabels(sortedGeneNames, minor=False)
		axmat.xaxis.set_label_position('bottom')
		axmat.xaxis.tick_bottom()

		pylab.xticks(rotation=270, fontsize=1)

		"""axmat.set_yticks(range(40))
		axmat.set_yticklabels(idx2, minor=False)
		axmat.yaxis.set_label_position('right')
		axmat.yaxis.tick_right()"""

		axcolor = fig.add_axes([0.94,0.1,0.02,0.6])


	print geneNames
	# Compute and plot first dendrogram.
	fig = pylab.figure(figsize=(8,8))
	ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
	Y = sch.linkage(D, method='centroid')
	Z1 = sch.dendrogram(Y, orientation='right')
	ax1.set_xticks([])
	ax1.set_yticks([])

	leafOrder = sch.leaves_list(Y)

	# Compute and plot second dendrogram.
	"""ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
	Y = sch.linkage(D, method='single')
	Z2 = sch.dendrogram(Y)
	ax2.set_xticks([])
	ax2.set_yticks([])"""

	# Plot distance matrix.
	axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
	idx1 = Z1['leaves']
	#idx2 = Z2['leaves']
	D = D[idx1,:]
	#D = D[:,idx2]
	im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)

	# Plot colorbar.
	axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
	pylab.colorbar(im, cax=axcolor)
	create_dendogram_ticks(axmatrix, leafOrder, geneNames)
	fig.show()
	saveFileName = 'dendrogram.png'
	print 'saving file to ', saveFileName
	fig.savefig(saveFileName)

def create_distance_matrix(df):

	genes = set(df['hgnc_symbol'])
	dimension = len(genes)
	df = df.drop_duplicates(subset=['hgnc_symbol']) #remove duplicate entries for genes
	orderedGeneNames = list(df['hgnc_symbol'])
	df = df.drop(['ensembl_transcript_id', 'hgnc_symbol'], axis=1)
	D = scipy.zeros([dimension,dimension])

	df = df.div(df.sum(axis=1), axis=0) #normalize df

	dfCopy = df.copy()
	i = 0
	for index, row in df.iterrows():
		j = 0
		l = list(row)
		for index1, row1 in dfCopy.iterrows():
			lAlt = list(row1)
			dist = 1 - spatial.distance.cosine(l, lAlt)
			D[i, j] = dist
			j += 1
		i += 1
	return D, orderedGeneNames

def make_histograms(df, saveDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/graphs'):
	fig, ax = plt.subplots()
	cntr = 0
	for index, row in df.iterrows():
		name = row['hgnc_symbol']
		rowAsDict = row.to_dict()
		del rowAsDict['ensembl_transcript_id']
		del rowAsDict['hgnc_symbol']
		l = []
		for key, value in rowAsDict.items():
			for i in range(value):
				l.append(key)

		labels, values = zip(*Counter(l).items())
		indicies = np.array([i*8 for i in range(len(labels))])
		print indicies
		width = 1

		plt.bar(indicies, values, width)
		plt.xlim([0,1000])
		plt.xticks(indicies, labels)

		#plt.hist(l)
		plt.title(name)
		plt.xlabel('Trinuc')
		plt.xticks(rotation=45)
		plt.ylabel("Frequency")
		saveFileName = name +'_graph.png'
		savePath = os.path.join(saveDir, saveFileName)
		plt.savefig(savePath)
		print 'saving', savePath
		cntr += 1
		sys.exit()

		fig = plt.gcf()

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--trinucsTable', help='stub for parser argument', default='/ifs/work/taylorlab/friedman/scpReceivingStation/gene_reference_signatures_allframe.tsv')

	args = parser.parse_args()

	genePanelGenes = get_gene_panel_genes()
	df = pd.read_table(args.trinucsTable)
	df = df[df['hgnc_symbol'].isin(genePanelGenes)]

	D, geneNames = create_distance_matrix(df)
	dendogram_function(D, geneNames)
	#make_histograms(df)

if __name__ == '__main__':
    main()















