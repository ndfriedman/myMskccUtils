#written by Noah Friedman 
#a suite of processing functions to prepare myriad data into the formats I need to make signature landscape plots


import sys
import argparse
import os
import pandas as pd
import numpy as np
import math

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

#FUNCTIONALITY for interactive python on my desktop
pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'
    

def find_second_most_common(row, primarySig, returnMode, 
                            sigNamesToSpecify = set(['mean_1', 'mean_3', 'mean_4', 'mean_7', 'mean_10', 'mean_11','mean_14', 'mean_17', 'mean_MMR', 'mean_APOBEC']) #a set of signatures we actually mark on the chart
                            ):
    colNames = row.to_dict().keys()
    signatureColumns = [i for i in list(row.keys()) if 'mean' in i]
    rowSigsOnly = row[signatureColumns]
    rowAsDict = rowSigsOnly.to_dict()
    items = rowAsDict.items()
    sortedItems = sorted(items, key=lambda x: x[1], reverse=True)
    if sortedItems[0][0] == primarySig:
        if returnMode == 'name':
            sigName = sortedItems[1][0]
            if sigName in sigNamesToSpecify:
                return sigName
            else:
                return 'other'
        else:
            return sortedItems[1][1]
    else:
        if returnMode == 'name':
            sigName = sortedItems[0][0]
            if sigName in sigNamesToSpecify:
                return sigName
            else:
                return 'other' 
        else:
            return sortedItems[0][1]


#define an ordering for the ggplot plot based on the first signature to hit a clipping thresh
#add functionality for low mut burden
def ordering_function_clip_mode(row, clipThresh = .15, sigsToOrderBy = ['mean_MMR', 'mean_1', 'mean_APOBEC', 'mean_3', 'mean_4', 'mean_7', 'mean_10']):
    orderingNum = len(sigsToOrderBy)
    for i in range(len(sigsToOrderBy)):
        curSigToConsider = sigsToOrderBy[i]    
        if row[curSigToConsider] > clipThresh and (curSigToConsider == row['otherPredominantSigName'] or curSigToConsider == 'mean_MMR'):
            return orderingNum + row[curSigToConsider]
        else:
            orderingNum -=1
    return 0


#an ordering function for a df     
def ordering_function_dom_sig_mode(row, domSig, #the signature to have as the primary ordering
	ageSigReOrderMode = False, sigsToOrderBy = ['mean_APOBEC', 'mean_3', 'mean_4', 'mean_7', 'mean_10', 'mean_14', 'mean_11', 'mean_17']):
    orderingNum = len(sigsToOrderBy) + 2
    if row['Nmut'] < 10 or math.isnan(row['Nmut']): return -1 #cases with less than 10 mutations go at the far side
    
    if row[domSig] > row['otherPredominantSigMagnitude']:
        return orderingNum + row[domSig]
    orderingNum -= 1
    if row['otherPredominantSigName'] == 'mean_1' or row['otherPredominantSigName'] == 'Age' or row['otherPredominantSigName'] == '1': 
    	if ageSigReOrderMode:
        	return orderingNum + row[domSig] #if age sig reorder mode reorer the age sig column by the dominant signature
        else:
        	return orderingNum + row['mean_1']
    orderingNum -= 1
    for i in range(len(sigsToOrderBy)):
        curSigToConsider = sigsToOrderBy[i]
        if curSigToConsider == row['otherPredominantSigName']:
        	return orderingNum + row[curSigToConsider]
        else:
            orderingNum -=1
    return 0

def ordering_function_two_signatures_mode(row, ordering=None):
    orderingNum = len(ordering)
    if row['Nmut'] < 10: return -1
    for i in range(len(ordering)):
        curSigToConsider = ordering[i]
        if curSigToConsider == row['dominantSignatureName']:
            return orderingNum + row[curSigToConsider]
        else:
            orderingNum -=1
    return 0

def temp_func():
    return 0

#functions for data prep

def rename_column_for_aesthetic_purposes():
	return 0

def rename_columns_from_philip_data():
	return 0











