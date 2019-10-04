#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np
import imp

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

analysis_utils = imp.load_source('analysis_utils', '/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/analysis_utils.py')


#creates the reference four nucleotide context for signatures
def create_reference_four_nuc(refTri, refAllele, altAllele):
	
	#properly invert when needed
	def invert_allele_and_ref_tri(altAllele):
		nucleotideDict = {'A': 'T', 'G': 'C', 'C':'G', 'T':'A'}
		return nucleotideDict[altAllele]

	refAlleleFromReftri = refTri[1]
	alt = altAllele
	if refAlleleFromReftri != refAllele:
		alt = invert_allele_and_ref_tri(altAllele)
	quadNuc = refTri[:2] + alt + refTri[2]
	return quadNuc

def normalize_probs(counterOfProbs):
	total = sum(counterOfProbs.values(), 0.0)
	for key in counterOfProbs:
		counterOfProbs[key] /= total
	return counterOfProbs


######NEW FUNCTION 6/13
#THIS marks cases whose motifs are consistent with the dominant signature pattern in the case (hypermutator motif)
#and mutations that are not (not hypermutator motif)
def mark_mutations_in_case_as_motif_outliers(maf):

	#TWO different ways to enumerate related quad nucs: quadnucs above a theshold
	#quad nucs representing ~Threshold amount of probability mass of quad nucs

	def enumerate_related_quadnucs(normedQuadNucFracs, rawQuadNucCounts):
		#quadNucsEnrichedForSignature = {x : normedCounter[x] for x in normedCounter if normedCounter[x] >= thresh}
		
		RELATED_MAX_FRAC = .75 #THE FRACTION ABOVE WHICH WE WILL START TO CONSIDER QUAD NUCS unrelated
		MIN_QUAD_NUC_THESH = .03 #the minimum fraction of quad nucs that can be done here
		MIN_QUAD_NUC_COUNT = 3 #the minimum number of mutations at a quad nuc tgat can occur
		quadNucFracSum = 0

		relatedQuadnucs = []
		for quadNuc, val in normedCounter.most_common(len(normedCounter)):
			quadNucFracSum += val
			if quadNucFracSum < RELATED_MAX_FRAC and val >= MIN_QUAD_NUC_THESH and rawQuadNucCounts[quadNuc] >= MIN_QUAD_NUC_COUNT:
				relatedQuadnucs.append(quadNuc)
			else:
				return relatedQuadnucs
	
		return relatedQuadnucs

	quadNucAetiologyDict = {}
	for case in set(maf['Tumor_Sample_Barcode']):
		caseMafWithQuadNucs = maf[(maf['Tumor_Sample_Barcode'] == case) & (maf['quadNuc'].notnull())]
		quadNucs = caseMafWithQuadNucs['quadNuc']
		
		#TRY A NON theshold based method?
		quadNucCounts = Counter(quadNucs)
		normedCounter = analysis_utils.normalize_counter(quadNucCounts, nDigitsRound=3)
		relatedQuadnucs = enumerate_related_quadnucs(normedCounter, quadNucCounts)
		quadNucAetiologyDict[case] = relatedQuadnucs
	maf['isRelatedQuadnuc'] = maf.apply(lambda row:
		None if row['quadNuc'] == None
		else True if row['quadNuc'] in quadNucAetiologyDict[row['Tumor_Sample_Barcode']]
		else False, axis=1)
	return maf
		
#given a row delineating signatures in a case goes bit by bit and assigns which signature likely caused a mutation
def assign_mutation_signature_probabilities(row, spectrumDict, quadNuc, signatureConsiderationMode = 'all', mode='singleAnswer'):
	signatures = ['Signature.' + str(i) for i in range(1,31)]

	#signatures = signatures[:18] #code to get rid of the "garbage" signatures
	#signatures += ['Signature.20', 'Signature.26', 'Signature.29']
	#signatures.remove('Signature.10')
	#signatures.remove('Signature.11')
	#signatures.remove('Signature.8')
	
	d = {}
	for signature in signatures:
		spectrumValue = spectrumDict[signature][quadNuc]
		caseValueForSignature = row[signature]
		
		if caseValueForSignature < 0.01: caseValueForSignature = 0

		d[signature] = spectrumValue * caseValueForSignature
	normalizedProbs = normalize_probs(Counter(d))
	if mode == 'singleAnswer':
		return normalizedProbs.most_common(1)[0][0]
	elif mode == 'fullDict':
		return normalizedProbs
	else:
		print 'error imporoper mode specified'
		sys.exit()

#function to annotate a mutation as a possible snp driver (a snp at a hotspot or oncogenic)
#Note requires a maf with variant type, oncogenic and is_a_hotspot columns
def annotate_possible_snp_driver(row): 
	#if its a hotspot snp
	if row['Variant_Type'] != 'SNP':
		return 'Not SNP'
	elif row['is-a-hotspot'] == 'Y':
		return True
	elif row['oncogenic'] == 'Likely Oncogenic' or row['oncogenic'] == 'Oncogenic' or row['oncogenic'] == 'Predicted Oncogenic':
		return True
	else:
		return False


def calculate_driver_fraction(currentSampleSigMotifs, signatureName, signatureMotifDict):
    curSignatureMotifs = signatureMotifDict[signatureName]
    nDriversAtMotif = 0.0 #counters for number of signatures
    nDrivers = len(currentSampleSigMotifs)
    for motif in currentSampleSigMotifs:
    	if motif in curSignatureMotifs:
    		nDriversAtMotif += 1
    	nDrivers += 1
    return nDriversAtMotif/nDrivers


#POSSIBLY deprecated functions
def get_cytosine_deanimation_pct(mafPath): #returns percentage of snps in the classic aging cytosine deanimation bins
	
	def getAgingSigFraction(trinucCntr): #gets the fraction that happens at the four peaks of the aging signature
		v = 0
		if 'ACTG' in trinucCntr:
			v += trinucCntr['ACTG']
		if 'CCTG' in trinucCntr:
			v += trinucCntr['CCTG']
		if 'GCTG' in trinucCntr:
			v += trinucCntr['GCTG']
		if 'TCTG' in trinucCntr:
			v += trinucCntr['TCTG']
		return v

	def get_c_t_fraction(trinucCntr): #get the fraction of ct transversions among all mutations
		ctMotifs = ['ACTA', 'ACTC', 'ACTG', 'ACTT', 'CCTA', 'CCTC', 'CCTG', 'CCTT', 'GCTA', 'GCTC', 'GCTG', 'GCTT', 'TCTA', 'TCTC', 'TCTG', 'TCTT']
		v = 0
		for motif in ctMotifs:
			v += trinucCntr[motif]
		return v

	maf = pd.read_table(mafPath)
	mafSnpsOnly = maf[maf['Ref_Tri'].notnull()]
	mafSnpsOnly['fourNuc'] = mafSnpsOnly.apply(lambda row: create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2']), axis=1)
	nMuts = mafSnpsOnly.shape[0]
	normedCntr = normalize_counter(Counter(mafSnpsOnly['fourNuc']))
	return nMuts, getAgingSigFraction(normedCntr), get_c_t_fraction(normedCntr)


def init_spectrum_dict(spectrumFile='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt'): #one liner to create a spectrum dict for later use
	spectrumDict = mutationSigUtils.convert_spectrum_file_to_dict_of_dicts(spectrumFile)
	return spectrumDict

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()

if __name__ == '__main__':
    main()















