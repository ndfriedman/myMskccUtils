#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np


def enumerate_paths(rootPath, dirs, writeToDiskEveryNIterations=5000):
		listOfLists = []
		initialPaths = np.array(['*'.join(['+' for i in range(50)]) for i in range(writeToDiskEveryNIterations)], dtype="S100")
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

				if cntr%writeToDiskEveryNIterations ==0:
					return

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

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')



	args = parser.parse_args()

if __name__ == '__main__':
    main()















