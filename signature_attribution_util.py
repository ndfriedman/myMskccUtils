#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

#from myUtils import mutationSigUtils



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

#given a row delineating signatures in a case goes bit by bit and assigns which signature likely caused a mutation
def assign_mutation_signature_probabilities(row, spectrumDict, signatureConsiderationMode = 'all', mode='singleAnswer'):
	quadNuc = create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'])
	signatures = ['Signature.' + str(i) for i in range(1,31)]

	signatures = signatures[:18] #code to get rid of the "garbage" signatures
	signatures += ['Signature.20', 'Signature.26', 'Signature.29']
	signatures.remove('Signature.10')
	signatures.remove('Signature.11')
	signatures.remove('Signature.8')
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















