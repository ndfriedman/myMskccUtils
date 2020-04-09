import pandas as pd
import sys
import os

rDataFilenamesPath = '/juno/work/taylorlab/friedman/myAdjustedDataFiles/rdata_filenames.txt'
rDataWritePath = '/juno/work/taylorlab/friedman/myAdjustedDataFiles/mafAnnoRuns/rDataFilenames'
rDataFilenames = pd.read_table(rDataFilenamesPath)

print 'loading unannotatedMaf'
unannotatedMaf = pd.read_table('/juno/work/taylorlab/friedman/myAdjustedDataFiles/dataMutationsUnfilteredNov19.txt', skiprows=[0])
fileSize = 1000
dfDict = {}
for i in range(0, rDataFilenames.shape[0]/fileSize):
	dfDict[i] = rDataFilenames[(rDataFilenames.index >= i*fileSize) & (rDataFilenames.index < (i + 1)*fileSize)]

writeLines = []
commandPrefix = 'bsub -We 59 -o /home/friedman/LSF/ -e /home/friedman/Err/ -W 8:00 -n 4 -R "rusage[mem=4]" Rscript /juno/work/taylorlab/friedman/myUtils/runAnnotateMaf.R'
cntr = 0
for key, values in dfDict.items():
	print cntr
	cntr += 1
	#Write the R data file
	rDataFilename = os.path.join('/juno/work/taylorlab/friedman/myAdjustedDataFiles/mafAnnoRuns/rDataFilenames', str(key) + '_th_group_rDataFilenames.txt')
	values.to_csv(rDataFilename, index=False, sep='\t')

	#Write a smaller maf
	mafPath = os.path.join('/juno/work/taylorlab/friedman/myAdjustedDataFiles/mafAnnoRuns/unannotatedMafs', str(key) + '_th_group.maf')
	tids = values['Tumor_Sample_Barcode']
	maf = unannotatedMaf[unannotatedMaf['Tumor_Sample_Barcode'].isin(tids)]
	mafPath = os.path.join('/juno/work/taylorlab/friedman/myAdjustedDataFiles/mafAnnoRuns/unannotatedMafs', str(key) + '_th_group.maf')
	maf.to_csv(mafPath, index=False, sep='\t')

	writePath = os.path.join('/juno/work/taylorlab/friedman/myAdjustedDataFiles/mafAnnoRuns/annotatedMafs', str(key) + '_th_group_cncf_annotated.maf')
	writeLines.append(' '.join([commandPrefix,rDataFilename, mafPath, writePath]))

shFilePath = '/juno/work/taylorlab/friedman/myUtils/runMafAnno.sh'
shFile = open(shFilePath, 'w')
shFile.write('#!/bin/bash' + '\n')
for line in writeLines:
	shFile.write(line + '\n')
print 'instructions written to ', shFilePath
	