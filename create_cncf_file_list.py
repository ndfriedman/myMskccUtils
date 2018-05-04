#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np
import re

from multiprocessing import Process, Queue

def write_cncf_file(lines, headerLine=['Tumor_Sample_Barcode',  'Rdata_filename'], 
	writeFileDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles',
	writeFileName='cncf_filenames.txt',
	mode='writeNewFile',
	alreadyWrittenFile=None):
	
	if mode == 'writeNewFile': #only writes the header
		writePath = os.path.join(writeFileDir, writeFileName) 
		f = open(writePath, 'w')
		f.write('\t'.join(headerLine) + '\n')
		for line in lines:
			slashIndicies = [m.start() for m in re.finditer('\\/', line)]
			secondToLastSlashIndex = slashIndicies[len(slashIndicies) - 2]
			tBarcode = line[secondToLastSlashIndex + 1: secondToLastSlashIndex + 18] #magic indexing to extract the tumor sample barcode (alert not robust)
			f.write('\t'.join([tBarcode, line]) + '\n')

	elif mode == 'appendToExistingFile':
		if alreadyWrittenFile == None:
			print 'error if mode is append to existing file you need to specify the existing file'
			sys.exit()
		else:
			f = open(alreadyWrittenFile, 'a')

import gc

def find_all_matching_files(rootPath, dirs, showStats=True):

	#We have to perform a lot of byzantine memory gymnastics to keep this program running
	def enumerate_paths(rootPath, dirs, writeNLines=500):
		listOfLists = []
		initialPaths = np.array(['*'.join(['+' for i in range(50)]) for i in range(writeNLines)], dtype="S100") #we do this to placate luna's fussiness about memory and python's sloppiness with memory
		paths = initialPaths
		index = 0 #index for writing to our paths array

		nDirectories = 0
		nFacetsDirsConsidered = 0
		cntr = 0 #help to understand runtime
		iterCntr = 0 #help to understand runtime in case some dirs have. bazillions of files
		for d in dirs:
			cntr += 1 
			if cntr %100 == 0:
				print 'considering the ', cntr, 'th entry in the facets directory out of ', len(dirs), ' total entries in the facets dir'
				print iterCntr, ' iterations involved over the last hundred files'
				iterCntr = 0
				gc.collect()

				if cntr==writeNLines:
					paths = paths[:index]
					return paths

			dirPath = os.path.join(rootPath, d)
			if os.path.isdir(dirPath):
				nDirectories += 1
				#facets = [x for x in os.listdir(dirPath) if 'facets_R0.5.6s' in x][0] #list comprehension (commented bc theoretically slower)
				facets = None
				for x in os.listdir(dirPath):
					if 'facets_R0.5.6s' in x:
						facets = x
						break #save time

				if facets != None: 
					facetsPath = os.path.join(dirPath, facets)
					if os.path.isdir(facetsPath):
						nFacetsDirsConsidered += 1
						for x in os.listdir(facetsPath):
							iterCntr += 1
							if '_hisens.cncf.txt' in x:
								cncfFilePath = os.path.join(facetsPath, x)
								paths[index] = cncfFilePath
								index += 1

								break #save time
		return paths

	paths = enumerate_paths(rootPath, dirs)
	write_cncf_file(paths, mode='writeNewFile')

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--facetsDirPath', help='path to the dir for facets output', default='/ifs/res/taylorlab/impact_facets/all')
	parser.add_argument('--keyFilePath', help='path to file with keys', default='/ifs/res/taylorlab/dmp_mirror_key/2018_05_01/dmp_key_paired.txt')



	args = parser.parse_args()

	dirs = os.listdir(args.facetsDirPath)
	find_all_matching_files(args.facetsDirPath, dirs)

if __name__ == '__main__':
    main()















