#written by Noah Friedman 
#USAGE:

#This script takes in a tsv file with signatures and confidences and adds a column with "signature interpretations"

import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

#a function that will return all the names of signatures that we believe are present
def get_signatures_we_believe_are_present(row, threshold= 0.85, nSignatures = 30):
	confidenceCols = ['confidence_' + str(i) for i in range(1,nSignatures + 1)]
	row = row[confidenceCols]
	detectedSigs = [str(i) for i in range(1, nSignatures + 1) if (row[i - 1] > threshold)] #recover the indicies of signatures we are confident are present
	return detectedSigs

#uses a buttload of if statements to generate the message for the signature
#uses if statements because some signatures have extra processing
def get_interpretation_message_for_signature(signature, nmut, signaturesMode = 'Stratton'):
	
	if signaturesMode == 'Stratton':
		if signature == '1':
			if nmut < 100:
				return 'Signature 1, the aging signature, is detected in this case.'
			else:
				return 'Signature 1, the aging signature, is detected in this case. However, because the high mutation burden in this case, we believe this signature may actually be caused by MSI/MMR dsyfunction (signature 6)'

		if signature == '2':
			return 'Signature 2, the APOBEC signature, is detected in this case.  This signature often coccurs with signature 13, the other APOBEC signature'

		if signature == '3':
			return 'Signature 3, the signature of Homologous Recombination Repair deficiency is detected in this case.  This signature is most commonly associated with BRCA mutations'

		if signature == '4':
			return 'Signature 4, the smoking signature is detected in this case'

		if signature == '5':
			return 'Signature 5 is detected in this case.  We are not confident that we are able to detect signature 5 in the IMPACT cohort.  It is a "flat" signature--when it is detected it is more likely to be an artefact. In the literature it is associated with age'

		if signature == '6':
			return 'Signature 6, a MMR signature, is detected in this case.  It is usually associated with high mutational burden.  This signature often co-occurs with other MMR signatures (14, 15, 20, 21 26)'

		if signature == '7':
			return 'Signature 7, the UV light signature, is detected in this case.'

		if signature == '8':
			return 'Signature 8 is detected in this case. We are not confident that we are able to detect signature 8 in the IMPACT cohort. It is a "flat" signature--when it is detected it is more likely to be an artefact. In the literature it is associated with HRD defects'

		if signature == '9':
			return 'Signature 9 is detected in this case.  We are not confident that we are able to detect signature 9 in the IMPACT cohort.  In the literature it is associated with POLH'

		if signature == '10':
			return 'Signature 10, the POLE signature, is detected in this case.  It is associated with functions to the exonucleus domain of the POLE gene and enormous mutational burden.  Oftentimes MMR signatures 6, 14,16, 20,21 and 26 co-occur with the POLE signature'

		if signature == '11':
			return 'Signature 11, the Temozolomide (TMZ) signature, is detected in this case'

		if signature == '12':
			return 'Signature 12 is detected in this case.  We are not confident that we are able to detect signature 9 in the IMPACT cohort.  In the literature it is found in liver cancer'

		if signature == '13':
			return 'Signature 13, the APOBEC signature, is detected in this case.  This signature often coccurs with signature 2, the other APOBEC signature'

		if signature == '14':
			return 'Signature 14, the signature of simultaneous MMR and POLE dysfunction is detected in this case.  This signature usually occurs in cases with the POLE signature (signature 10) and other MMR signatures (6, 15, 20, 21 26)'

		if signature == '15':
			return 'Signature 15, a MMR signature, is detected in this case.  It is usually associated with high mutational burden'

		if signature == '16':
			return 'Signature 16 is detected in this case. We are not confident that we are able to detect signature 16 in the IMPACT cohort.  In the literature it is associated with Liver cancer and alcohol consumption'

		if signature == '17':
			return 'Signature 17 is detected in this case.  The aetiology of this signature is unknown.  It is predominantly found in gastric cancers'

		if signature == '18':
			return 'Signature 18 is detected in this case.  This signature is associated with MUTYH dysfunction and neuroblastoma'

		if signature == '19':
			return 'Signature 19 is detected in this case. We are not confident that we are able to detect signature 19 in the IMPACT cohort'

		if signature == '20':
			return 'Signature 20 is detected in this case. This signature is associated with MMR and usually occurs in cases with the POLE signature (signature 10) and other MMR signatures (6, 14, 15, 21, 26)'

		if signature == '21':
			return 'Signature 21 is detected in this case. This signature is associated with MMR and usually co-occurs with other MMR signatures (6, 14, 15, 21, 26)'

		if signature == '22':
			return 'Signature 22 is detected in this case. We are not confident that we are able to detect signature 22 in the IMPACT cohort. In the literature it is associated with exposure to Aristolochic Acid'

		if signature == '23':
			return 'Signature 23 is detected in this case. We are not confident that we are able to detect signature 23 in the IMPACT cohort'

		if signature == '24':
			return 'Signature 24 is detected in this case. We are not confident that we are able to detect signature 24 in the IMPACT cohort.  In the literature it is associated with aflatoxin exposure.  In our cohort we believe it is detected by accident in cases with the smoking signature (signature 4)'

		if signature == '25':
			return 'Signature 25 is detected in this case. We are not confident that we are able to detect signature 25 in the IMPACT cohort'

		if signature == '26':
			return 'Signature 26 is detected in this case.  This signature is associated with MMR and usually co-occurs with other MMR signatures (6, 14, 15, 20, 21)'

		if signature == '27':
			return 'Signature 27 is detected in this case. We are not confident that we are able to detect signature 27 in the IMPACT cohort'

		if signature == '28':
			return 'Signature 28 is detected in this case. We are not confident that we are able to detect signature 28 in the IMPACT cohort.  It often co-occurs with signature 28'

		if signature == '29':
			return 'Signature 29, the mutational signature of chewing tobacco is detected in this case'

		if signature == '30':
			return 'Signature 30 is detected in this case. We are not confident that we are able to detect signature 30 in the IMPACT cohort'


def add_signature_interpretation(row, 
	signaturesMode = 'Stratton',
	minNumberOfMutationsThreshold = 5, #below this number we will not even accept signatures
	dubiousNumberOfMutationsThreshold = 15, #below this number the user gets a disclaimer that there aren't very many mutations

	):



	"""
	This function is used to add a "signature interpretation" to a dataframe of cases and their respective signature decompositions
	We slowly build up an interpretation string which we eventually return to the client

	"""



	interpretationString = ''
	nmut = row['Nmut']

	if signaturesMode == 'Stratton':

		#FIRST add a note about mutation number:
		if nmut < minNumberOfMutationsThreshold:
			return interpretationString + 'There are not enough mutations to call mutation signatures for this case.'
		elif nmut < dubiousNumberOfMutationsThreshold:
			interpretationString += 'This case has a low number of mutations, signatures are highly uncertain. '
		
		#NOW get whichever signatures we feel confident about:
		presentSignatures = get_signatures_we_believe_are_present(row, threshold= 0.85)
		if len(presentSignatures) > 0:
			for signature in presentSignatures:
				interpretationString += get_interpretation_message_for_signature(signature, nmut)
		else:
			interpretationString += 'No signatures confidently detected'

		return interpretationString

	elif mode == 'SBS':
		print 'alert ERROR SBS signature mode has not been implemented'
		sys.exit()

	else:
		print 'alert ERROR improper signatures mode specified'
		sys.exit()

def main():

	parser = argparse.ArgumentParser(description='his script takes in a tsv file with signatures and confidences and adds a column with "signature interpretations"')
	parser.add_argument('--signaturesDf', help='Path to the dataframe with mutational signatures and confidence information', default='/ifs/work/taylorlab/pensona/dmp_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
	parser.add_argument('--mode', help='If you specify mode=test it will only run on 100 cases for quicker runtime', default='Normal')

	args = parser.parse_args()

	signaturesDf = None
	if args.mode == 'test':
		signaturesDf = pd.read_table(args.signaturesDf, nrows=100)
	else:
		signaturesDf = pd.read_table(args.signaturesDf)

	signaturesDf['signatureInterpretation'] = signaturesDf.apply(lambda row: add_signature_interpretation(row), axis=1)

	signaturesDf.to_csv('signaturesInformationWithInterpretation.tsv', sep='\t', index=False)

if __name__ == '__main__':
    main()















