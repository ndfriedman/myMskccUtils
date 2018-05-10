#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np
import re

import gc

from multiprocessing import Process, Queue

def write_cncf_file(lines, headerLine=['Tumor_Sample_Barcode',  'Rdata_filename'], 
	writeFileDir = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles',
	writeFileName='cncf_filenames.txt',
	mode='writeNewFile',
	alreadyWrittenFile=None):
	
	def write_lines(f, lines):
		for line in lines:
			pDashIndicies = [m.start() for m in re.finditer('P-', line)]
			if len(pDashIndicies) > 0: #we can only write things with properly formed ts barcodes
				tumorSampleBarcode = line[pDashIndicies[0] : pDashIndicies[0] + 17] #magic indexing; may be fragile
				f.write('\t'.join([tumorSampleBarcode, line]) + '\n')

	if mode == 'writeNewFile': #only writes the header
		print 'writing new file'
		writePath = os.path.join(writeFileDir, writeFileName) 
		f = open(writePath, 'w')
		f.write('\t'.join(headerLine) + '\n')
		write_lines(f, lines)
		f.close()
		return writePath

	elif mode == 'appendToExistingFile':
		print 'appending to existing file'
		if alreadyWrittenFile == None:
			print 'error if mode is append to existing file you need to specify the existing file'
			sys.exit()
		else:
			f = open(alreadyWrittenFile, 'a')
			write_lines(f, lines)
			f.close()
			return alreadyWrittenFile


def find_all_matching_files(rootPath, dirs, mode, showStats=True):

	#We have to perform a lot of byzantine memory gymnastics to keep this program running
	def enumerate_paths(rootPath, dirs, writeNLines=1000, runtimeCntr=0):
		listOfLists = []
		initialPaths = np.array(['*'.join(['+' for i in range(50)]) for i in range(writeNLines)], dtype="S200") #we do this to placate luna's fussiness about memory and python's sloppiness with memory
		paths = initialPaths
		index = 0 #index for writing to our paths array

		nDirectories = 0
		nFacetsDirsConsidered = 0
		iterCntr = 0 #help to understand runtime in case some dirs have. bazillions of files
		for d in dirs:
			runtimeCntr += 1 
			if runtimeCntr %500 == 0:
				print 'considering the ', runtimeCntr, 'th entry in the facets directory out of ', len(dirs), ' total entries in the facets dir'
				print iterCntr, ' iterations involved over the last hundred files'
				iterCntr = 0
				gc.collect()

				if runtimeCntr%writeNLines == 0: #once we've considered n lines, write to the file and keep going:
					paths = paths[:index]
					return paths, runtimeCntr

			dirPath = os.path.join(rootPath, d)
			if os.path.isdir(dirPath):
				nDirectories += 1
				#facets = [x for x in os.listdir(dirPath) if 'facets_R0.5.6s' in x][0] #list comprehension (commented bc theoretically slower)
				facets = None
				for x in os.listdir(dirPath):
					if 'facets_R0.5.6s' in x:
						facets = x
						break #save time

					#Add in functionality for if the s version doesnt exist go and get the one without 

				if facets != None: 
					facetsPath = os.path.join(dirPath, facets)
					if os.path.isdir(facetsPath):
						nFacetsDirsConsidered += 1
						for x in os.listdir(facetsPath):
							iterCntr += 1
							if mode == 'createNewCncfFile':
								if '_hisens.cncf.txt' in x:
									cncfFilePath = os.path.join(facetsPath, x)
									paths[index] = cncfFilePath
									index += 1
									break #save time
							else:
								if 'hisens.Rdata' in x:
									rDataFilePath = os.path.join(facetsPath, x)
									paths[index] = rDataFilePath
									index += 1
									break #save time
		return paths, runtimeCntr

	writeNLines=1000
	#we do it once outside the while loop to write a new file than work inside the while loop to append to the existing file
	paths, cntr = enumerate_paths(rootPath, dirs, writeNLines)
	file = None
	if mode == 'createNewCncfFile':
		file = write_cncf_file(paths, mode='writeNewFile', writeFileName='cncf_filenames.txt')
	else:
		file = write_cncf_file(paths, mode='writeNewFile', writeFileName='rdata_filenames.txt')
	dirs = dirs[writeNLines:]
	while(len(dirs) > 0):
		paths, cntr = enumerate_paths(rootPath, dirs, writeNLines, runtimeCntr=cntr)
		write_cncf_file(paths, mode='appendToExistingFile', alreadyWrittenFile=file)
		dirs = dirs[writeNLines:]


def add_new_ids_to_cncf_file(preexistingCncfFilePath= ''):

	def get_ids_from_exisiting_cncf():
		return 0


def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--facetsDirPath', help='path to the dir for facets output', default='/ifs/res/taylorlab/impact_facets/all')
	parser.add_argument('--keyFilePath', help='path to file with keys', default='/ifs/res/taylorlab/dmp_mirror_key/2018_05_01/dmp_key_paired.txt')
	parser.add_argument('--mode', help='mode to run this script in', default='createNewCncfFile')



	args = parser.parse_args()

	dirs = os.listdir(args.facetsDirPath)

	if args.mode == 'createNewCncfFile':
		find_all_matching_files(args.facetsDirPath, dirs, args.mode)

	if args.mode == 'createNewRdataFile':
		find_all_matching_files(args.facetsDirPath, dirs, args.mode)

if __name__ == '__main__':
    main()















