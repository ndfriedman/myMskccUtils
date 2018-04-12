#written by noah friedman friedman@mskcc.org

import sys
import os
import subprocess
import pandas as pd



#run code to get mutational signatures
def run_mutational_signatures_code(scriptDir, outputFilePath, sMaf, tMaf=None):
	def run_triunc_command(scriptDir, sMaf, tMaf):
		triuncScriptPath = os.path.join(scriptDir, 'make_trinuc_maf.py')
		cmd = 'python {sPath} {sourceMafPath} {TargetMafPath}'.format(sPath = triuncScriptPath, sourceMafPath = sMaf, TargetMafPath = tMaf)
		print 'executing Triunc command: ', cmd
		subprocess.Popen(cmd, shell=True).wait()
		print 'Triunc command completed'
		return tMaf

	def run_signatures_code(scriptDir, sMaf, outputFile):
		signaturesScriptPath = os.path.join(scriptDir, 'main.py')
		signaturesFilePath = '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt'
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
	print triuncMaf
	print 'RUMi'
	targetSignaturesFile = outputFilePath
	run_signatures_code(scriptDir, triuncMaf, targetSignaturesFile)

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





