#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np

import subprocess


def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--facetsDirPath', help='path to the dir for facets output', default='/ifs/res/taylorlab/impact_facets/all')



	args = parser.parse_args()

	print os.listdir(args.facetsDirPath)

if __name__ == '__main__':
    main()















