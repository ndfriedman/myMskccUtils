#written by Noah Friedman
#bsub -We 59 -o LSF/ -e Err/ -W 24 -n 4 -R "rusage[mem=4]" Rscript runAnnotateMaf.R
#bsub -We 59 -o LSF/ -e Err/ -We 8:00 -n 8 -R "rusage[mem=4]" Rscript runAnnotateMaf.R

library(facetsSuite, lib.loc = '/juno/work/ccs/bandlamc/software/R_libs/facetsSuite/2.0.1-beta/')
library(data.table)

args = commandArgs(trailingOnly=TRUE)
rDataPath <- args[1]
mafPath <- args[2]
writePath <- args[3]

#rDataPath <- '/juno/work/taylorlab/friedman/myAdjustedDataFiles/rdata_filenames.txt'
#mafPath <- '/juno/work/taylorlab/friedman/myAdjustedDataFiles/hypermutantOnlyUnfilteredMaf.txt'
#writePath <- '/juno/work/taylorlab/friedman/myAdjustedDataFiles/unfilteredMafHypermutatedOnlu_withCNCFAnnotation.maf'


#mafPath <- '' #small maf for testing?
maf <- read.csv(mafPath, sep = '\t', header=TRUE, fill = TRUE)

#NOW WE GO AND ITERATE OVER ALL THE CASES AND ONE BY ONE ADD CLONALITY INFO
rDataLines <- readLines(rDataPath)

cntr <- 1
listOfMafs <- list()
for(line in rDataLines[2:length(rDataLines)]){ #skip the first line of the file which is the column headers
	vals <- strsplit(line, '\t')
	tsb <- vals[1][[1]][1] #IM BAD AT R
	rDataPath <- vals[1][[1]][2]

	if(tsb %in% unique(maf$Tumor_Sample_Barcode)){

		caseMaf <- maf[maf$Tumor_Sample_Barcode == tsb,]

		print(paste('annotating case ', tsb, 'which is the ', cntr, 'th case out of ',
		 	min(length(rDataLines), length(unique(maf$Tumor_Sample_Barcode))), ' cases'))
		
		facets_output = facetsSuite::load_facets_output(rDataPath)
		if(!is.na(facets_output$purity) #THE CCF annotate maf program does not work on cases with purity equals NA
		 & nrow(caseMaf) > 0){ #THE CCF annotate maf program does not work on cases with no mutations

			annotatedMaf <- facetsSuite::ccf_annotate_maf(caseMaf, facets_output$segs, facets_output$purity, algorithm = 'em')

			listOfMafs[[cntr]] = annotatedMaf 
		}
		else{
			listOfMafs[[cntr]] = caseMaf
		}
		cntr <- cntr + 1
		}	
}
print('combining all mafs together')
combinedMafs <- rbindlist(listOfMafs, fill=TRUE)

print(paste('writing file to', writePath))
write.table(combinedMafs, file=writePath, quote=FALSE, sep='\t')

