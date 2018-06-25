#written by noah friedman friedman@mskcc.org

import sys
import os
import subprocess
import pandas as pd
import argparse


#run code to get mutational signatures
def run_mutational_signatures_code(scriptDir, outputFilePath, sMaf, signaturesFilePath, triuncOnly=False, tMaf=None):
	def run_triunc_command(scriptDir, sMaf, tMaf):
		triuncScriptPath = os.path.join(scriptDir, 'make_trinuc_maf.py')
		cmd = 'python {sPath} {sourceMafPath} {TargetMafPath}'.format(sPath = triuncScriptPath, sourceMafPath = sMaf, TargetMafPath = tMaf)
		print 'executing Triunc command: ', cmd
		subprocess.Popen(cmd, shell=True).wait()
		print 'Triunc command completed'
		return tMaf

	def run_signatures_code(scriptDir, sMaf, outputFile, signaturesFilePath):
		signaturesScriptPath = os.path.join(scriptDir, 'main.py')
		cmd = 'python {sPath} {sigFilePath} {sourceMafPath} {targetFilePath} --spectrum_output spectrumNoah.txt'.format(sPath = signaturesScriptPath, sigFilePath = signaturesFilePath, sourceMafPath = sMaf, targetFilePath = outputFile)
		print 'executing Signatures command: ', cmd
		subprocess.Popen(cmd, shell=True).wait()
		print 'signatures command completed'

	triuncMaf = ''
	print sMaf
	print 'Sirt'
	if 'triunc' in sMaf: triuncMaf = sMaf
	else:
		triuncMaf = run_triunc_command(scriptDir, sMaf, tMaf)
		if triuncOnly:
			print 'triunc only mode, returning'
			return
	print triuncMaf
	print 'RUMi'
	targetSignaturesFile = outputFilePath
	run_signatures_code(scriptDir, triuncMaf, targetSignaturesFile, signaturesFilePath)

#sanity check that the code works by summing up the counts
def sanity_check_mutation_counts(signaturesFilePath):
	#util function for getting sample names from the key name
	#takes the sample name, and a key, 0 or 1 referring to which value we want
	def extract_sample_name(sampleNameComparison, index):
		if '!' in sampleNameComparison:
			sampleName = sampleNameComparison.split('!')[index]
		elif '&' in sampleNameComparison:
			sampleName = sampleNameComparison.split('&')[index]
		else:
			print 'Major error! sample name improperly formatted'
			sys.exit()
		return sampleName


	def create_sanity_check_dictionary(signaturesDf): #create a dictionary that we will use to do the sanity check
		sanityCheckDict = dict()
		for index, row in signaturesDf.iterrows(): 
			sampleNameComparison = row['Sample Name'] 
			sampleName = extract_sample_name(sampleNameComparison, 0)
			if sampleName in sanityCheckDict:
				curListOfRows = sanityCheckDict[sampleName]
				curListOfRows.append(row)
				sanityCheckDict[sampleName] = curListOfRows
			else:
				sanityCheckDict[sampleName] = [row]

		return sanityCheckDict

	#for each entry in the dictionary we create another dictionary and sanity check it
	def perform_sanity_check_on_dict(sanityCheckDict):
		for key, value in sanityCheckDict.items():
			subDict = dict()
			for v in value:
				sampleNameComparison = v['Tumor_Sample_Barcode']
				sampleName = extract_sample_name(sampleNameComparison, 1)
				if sampleName in subDict:
					l = subDict[sampleName]
					l.append(v['Number of Mutations'])
					subDict[sampleName] = l
				else:
					subDict[sampleName] = [v['Number of Mutations']]
			prevS = -1
			for key1, v in subDict.items():
				v1, v2 = v
				s = v1 + v2
				if prevS != -1 and s != prevS:
					print 'error on key:', key
					print v1, ' + ', v2, '!=', prevS
					sys.exit()
				prevS = s
			print 'sums validated for ', key
	
	#makes logic easier by duplicating and row data--ie if we have 1&2 we also create 2&1 which is identical but makes logic and code easier to manage
	def add_extra_and_rows(df):
		for index, row in df.iterrows():
			if '&' in row['Sample Name']:
				v1, v2 = row['Sample Name'].split('&')
				reversedSampleName = v2 + '&' + v1
				row['Sample Name'] = reversedSampleName
				df = df.append(row)
		return df

	print 'performing sanity check of mutation counts'
	signaturesDf = pd.read_table(signaturesFilePath)
	signaturesDf = signaturesDf[['Sample Name', 'Number of Mutations']]
	signaturesDf = add_extra_and_rows(signaturesDf)
	sanityCheckDict = create_sanity_check_dictionary(signaturesDf)
	perform_sanity_check_on_dict(sanityCheckDict)


def subset_triunc_signature_fractions(signature, nucleotides, spectrumFile='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt'):
	mutationSigTable = pd.read_table(spectrumFile)
	mutationSigTable['signature'] = mutationSigTable.index
	curSignature = mutationSigTable[mutationSigTable['signature'] == 'Signature.' + signature]
	s = 0
	for nuc in nucleotides:
		s += float(curSignature[nuc])
	return s

#UTILITIES FOR ASSIGING mutations to the signature that most likely caused them
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


def assign_most_likely_mutation(spectrumDict, row, n = 1):
	fourNuc = create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'])
	signatures = ['Signature.' + str(i) for i in range(1,31)]
	row = row[signatures]
	rowAsDict = row.to_dict()
	l = []
	for key, value in rowAsDict.items():
		l.append((key, value*spectrumDict[key][fourNuc]))
	signatures = [i[0] for i in l]
	values = [i[1] for i in l]
	sortedL = [x for _,x in sorted(zip(values, signatures), reverse=True)]
	return sortedL[:n] #returns the n most ocmmon signatures

# a little utility function to convert the mutations spectrum file to a dictionary of dictionaries for quicker lookup and access
def convert_spectrum_file_to_dict_of_dicts(spectrumFile='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt'):
	d = {}
	df = pd.read_table(spectrumFile)
	for index, row in df.iterrows():
		localD = row.to_dict()
		d[str(index)] = localD
	return d

def annotate_mutations_with_signatures_in_case(outputFilename, outputDir, mutationsFileToAnnotate='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/maf2mafAnnotatedMay16filteredMafWithIsHotspot.maf',
	signaturesFile='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mutationSigFiles/may16unfiltered30sigs.txt'
	):

	def convert_signatures_file_to_dict(signaturesFile):
		dictToDict = {}
		df = pd.read_table(signaturesFile)
		dfDictList = df.to_dict(orient='records')
		for row in dfDictList:
			barcode = row['Sample Name']
			del row['Number of Mutations']
			del row['Sample Name']
			dictToDict[barcode] = row
		return dictToDict

	barcodeToSignaturesDict = convert_signatures_file_to_dict(signaturesFile)
	mutationsDf = pd.read_table(mutationsFileToAnnotate)
	dfMutationsList = mutationsDf.to_dict(orient='records')
	listOfDicts = []
	for row in dfMutationsList:
		barcode = row['Tumor_Sample_Barcode']
		if barcode in barcodeToSignaturesDict:
			signaturesInfo = barcodeToSignaturesDict[barcode]
			rowCopy = row.copy()
			rowCopy.update(signaturesInfo)
			listOfDicts.append(rowCopy)
	df = pd.DataFrame(listOfDicts)
	writePath = os.path.join(outputDir, outputFilename)
	print 'writing file to ', writePath
	df.to_csv(writePath, sep='\t', index=False)


def create_limited_spectrum_file(signaturesToInclude, oldSpectrumFile='/ifs/work/taylorlab/friedman/myUtils/newSignatures.txt', outputDir='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/spectrumFiles'):
	spectrumDf = pd.read_table(oldSpectrumFile)
	spectrumDf = spectrumDf.ix[signaturesToInclude]
	writePath = os.path.join(outputDir, 'bladderSignatures.txt')
	print 'writing file to ', writePath
	spectrumDf.to_csv(writePath, index=True, sep='\t')

#####################UTILITIES FOR MERGING MUTATIONAL SIGNATURE COLUMNS

def merge_signature_columns(df, mode='Stratton'):
	if mode == 'Stratton':
		df['confidence_APOBEC'] = df.apply(lambda row: max(row['confidence_2'], row['confidence_13']), axis=1)
		df['mean_APOBEC'] = df.apply(lambda row: row['mean_2'] + row['mean_13'], axis=1)
		df['confidence_MMR'] = df.apply(lambda row: max(row['confidence_6'], row['confidence_15'], row['confidence_20'], row['confidence_26']), axis=1)
		df['mean_MMR'] = df.apply(lambda row: row['mean_6'] + row['mean_15'] + row['confidence_20'] + row['confidence_26'], axis=1)	
		df = df.drop(['confidence_2', 'confidence_13', 'confidence_6', 'confidence_15', 'confidence_20', 'confidence_26',
			'mean_2', 'mean_13', 'mean_6', 'mean_15', 'mean_20', 'mean_26'
			], axis=1)
		return df

	#TODO implement mean merges for SBS mode
	elif mode == 'SBS':
		df['confidence_APOBEC'] = df.apply(lambda row: max(row['confidence_SBS2'], row['confidence_SBS13']), axis=1)
		df['confidence_MMR'] = df.apply(lambda row: max(row['confidence_SBS6'], row['confidence_SBS15'], row['confidence_SBS20'], row['confidence_SBS20'], row['confidence_SBS26'], row['confidence_SBS44']), axis=1)
		df['confidence_BRCA'] = df.apply(lambda row: max(row['confidence_SBS3'], row['confidence_SBS39']), axis=1)
		df['confidence_UV'] = df.apply(lambda row: max(row['confidence_SBS7a'], row['confidence_SBS7b'], row['confidence_SBS7c'], row['confidence_SBS7d']), axis=1)
		df['confidence_POLE'] = df.apply(lambda row: max(row['confidence_SBS10a'], row['confidence_SBS10b']), axis=1)
		df['confidence_Sig17'] = df.apply(lambda row: max(row['confidence_SBS17a'], row['confidence_SBS17b']), axis=1)
		df = df.drop(['confidence_SBS2', 'confidence_SBS13', 
			'confidence_SBS6', 'confidence_SBS15', 'confidence_SBS20', 'confidence_SBS21', 'confidence_SBS26', 'confidence_SBS44', 
			'confidence_SBS3', 'confidence_SBS39',
			'confidence_SBS7a', 'confidence_SBS7b', 'confidence_SBS7c', 'confidence_SBS7d',
			'confidence_SBS10a', 'confidence_SBS10b',
			'confidence_SBS17a', 'confidence_SBS17b'], axis=1)
		return df
	else:
		print 'error improper mode specified'

def main():

	parser = argparse.ArgumentParser(description='Noahs script!')
	parser.add_argument('--inputMaf', help='maf to run signatures/triunc on', default='/ifs/work/taylorlab/friedman/clinicalData/msk-impact/msk-impact/data_mutations_extended.txt')
	parser.add_argument('--outputDir', help='output directory', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles')
	parser.add_argument('--outputFilename', help='output filename', default=None)
	parser.add_argument('--mutationalSignaturesScriptPath', help='path to the mutational signatures script', default='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures')
	parser.add_argument('--trinucOnly', help='mode for whether we just generate the triunc only file', default=False)
	parser.add_argument('--mode', help='mode for whether we just generate the triunc only file', default='triuncOnly')
	parser.add_argument('--spectrumFilePath', help='path to the spectrum file', default='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')  #signaturesFilePath = '/ifs/work/taylorlab/friedman/myUtils/newSignatures.txt'


	args = parser.parse_args()

	trinucOnly = False
	if args.mode == 'trinucOnly':
		args.mode = 'runMutSigs' #change mode to move with logic (alert confusing logic and dsylexic variable naming)
		trinucOnly = True

	if args.mode == 'annotateMutations':
		annotate_mutations_with_signatures_in_case(args.outputFilename, args.outputDir, mutationsFileToAnnotate=args.inputMaf)

	elif args.mode == 'createSpectrumFile':
		signaturesToInclude = ['SBS1', 'SBS2', 'SBS5', 'SBS10a', 'SBS10b', 'SBS13', 'SBS31', 'SBS35']
		create_limited_spectrum_file(signaturesToInclude)

	elif args.mode == 'runMutSigs':
		args.mutationalSignaturesOutputPath = os.path.join(args.outputDir, ''.join([args.outputFilename, 'mutationalSignatuesOutput.txt']))
		outputPath = os.path.join(args.outputDir, args.outputFilename)
		run_mutational_signatures_code(args.mutationalSignaturesScriptPath, args.mutationalSignaturesOutputPath, args.inputMaf, args.spectrumFilePath, triuncOnly=trinucOnly, tMaf=outputPath)

	else: print 'invalid mode specified', args.mode


if __name__ == '__main__':
    main()



