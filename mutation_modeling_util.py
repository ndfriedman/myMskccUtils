#written by Noah Friedman
#a collection of utilities designed to model mutation acquisition in silico 
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')


def init_gene_probabilities():
	return 0

#CURRENTLY DEPRECTED BUT IS A GOOD IDE A FOR WHAT SHIT WOULD LOOK LIKE IF WE RANDOMLY PICK SIGS
def get_trinucleotide_probabilities(trinucEntries):
		rowAsDict = trinucEntries.to_dict()
		trinucChoices = rowAsDict.keys()
		trinucCounts = rowAsDict.values()
		trinucProbabilities = [1.0*i/sum(trinucCounts) for i in trinucCounts]
		return trinucChoices, trinucProbabilities

#DEPRECATED???
"""def do_initiation(geneDistributions='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_trinuc_distributions.tsv'):
	
	df = pd.read_table(geneDistributions)
	quadNucColumns = df.columns.values[3:]
	df['NQuadNucInGene'] = df.apply(lambda row: sum(row[quadNucColumns]),axis=1)

	genes = set(df['Hugo_Symbol'])

	geneChoices = [] #the choices of gene mutations
	geneCounts = [] #the counts for each gene

	for gene in genes:
		geneChoices.append(gene)
		curGeneRow = df[df['Hugo_Symbol'] == gene].iloc[0]
		geneCounts.append(curGeneRow['NQuadNucInGene'])
		
	geneProbabilities = [1.0*i/sum(geneCounts) for i in geneCounts]
	return geneChoices, geneProbabilities"""


#One liner to pick a quadnuc given a signature that theoretically genrted it
def pick_quadnuc_given_spectra(spectraDict):
	return np.random.choice(
			spectraDict.keys(), p=spectraDict.values()
		)

#one liner to return a gene that is mutated given a quadNuc
def pick_gene_given_quadnuc(d, quadNuc):
	return np.random.choice(
			d[quadNuc][0], p=d[quadNuc][1]
		)


#BIG initiation function that intitates all the background probability information for all muts in impact
def do_initiation(geneDistributions='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_trinuc_distributions.tsv', mafWithOncogenicAnnotations):
	
	#given a maf of all muts in a gene, oncogenic muts and hotspot muts calculates the total fractions that occur at each gene
	#TODO in the future---? draw a specific oncogenic or hotspot mutation
	def do_oncogenic_mut_initiation(quadNuc, dictEntry):
		oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
		oncogenicMuts = mafWithOncogenicAnnotations[mafWithOncogenicAnnotations['oncogenic'].isin(oncogenicMutColNames)]
		print list(oncogenicMuts[oncogenicMuts['Hugo_Symbol'] == 'JAK1']['quadNuc'])
		#TODO use the background QUAD nuc distribution to quantify--what fraction of possible mutations in a quadnuc in a gene are ONCOGENIC

	def initiate_gene_mut_mapping(mafWithOncogenicAnnotations, geneList): #TO improve runtime we precalculate a dict mapping 
																#gene: maf of all mutations and maf of oncogenic mutations
		oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
		geneToMutDataDict = dict()
		for gene in geneList:
			geneMuts = mafWithOncogenicAnnotations[mafWithOncogenicAnnotations['Hugo_Symbol'] == gene]
			oncogenicGeneMuts = geneMuts[geneMuts['oncogenic'].isin(oncogenicMutColNames)]
			hotspotGeneMuts = geneMuts[geneMuts['is-a-hotspot'] == 'Y']
			geneToMutDataDict[gene] = (geneMuts, oncogenicGeneMuts, hotspotGeneMuts)
		return geneToMutDataDict


	df = pd.read_table(geneDistributions)
	quadNucColumns = df.columns.values[3:]
	df['NQuadNucInGene'] = df.apply(lambda row: sum(row[quadNucColumns]),axis=1)
	genes = list(df['Hugo_Symbol'])

	#FOR TESTING PURPOSES

	geneMutMappingDict = initiate_gene_mut_mapping(maf, genes)

	quadNucDict = dict()
	for qCol in quadNucColumns:
		counts = df[qCol]
		percentages = [1.0*c/sum(list(counts)) for c in list(counts)]
		quadNucDict[qCol] = (genes, percentages) #store a tuple of names (genes) and counts 
		
		#Now iterate through all the possible genes and calculate the probability P Oncogenic mutatiobn in Gene|mutation at motif
		for gene in genes:
			

	return quadNucDict


np.random.choice(
  ['pooh', 'rabbit', 'piglet', 'Christopher'], 
  5,
  p=[0.5, 0.1, 0.1, 0.3]
)












