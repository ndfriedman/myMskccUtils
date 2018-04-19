import sys
import os
from collections import Counter
import pandas as pd

#fixes the df to make the tumor sample barcode the patient id
def extract_patient_name(df, patientNameColumn='Tumor_Sample_Barcode'):
	df[patientNameColumn] = df[patientNameColumn].apply(lambda x: x[:9])
	return df

#util function to expand a column of a df
def expand_col(df, col):

	#utils lambda style function for making a column lowercase
	def to_lower_col(v):
		if not isinstance(v, float):
			return v.lower()
		else:
			return v

	def process_and_clean_entry(entry):
		s = []
		if not isinstance(entry, float): #AVOID Wasting our time with floats
			entry = entry.strip(' ')
			entry = entry.strip('\\+')
			if '/' in entry:
				s = entry.split('/')
			elif ',' in entry:
				s = entry.split(',')
			else:
				s = [entry]
		return s


	def identify_cols(df, col):
		l = list(df[col])
		colsList = []
		for v in l:
			for value in process_and_clean_entry(v):
				colsList.append(value)
		counts = Counter(colsList)
		counterThreshold = 50
		columnsSet = set()
		for key, count in counts.items():
			if count > counterThreshold:
				columnsSet.add(key)
		return columnsSet

	def assign_cols(df, newCols, originalCol):
		for col in newCols:
			print col
			df[col] = 0
			df[col] = df.apply(
				lambda row:
					0 if isinstance(row[originalCol], float) else 1 if col in row[originalCol] else 0 #convoluted lambda if else syntax
				,axis = 1)
		return df

	df[col] = df[col].apply(lambda x: to_lower_col(x))
	cols = identify_cols(df, col)
	df = assign_cols(df, cols, col)
	return df

#creates a PATIENT ID col for the df by doing reasoning on the tumor sample barcode col
def create_patient_id_col(df):
	df['PATIENT_ID'] = df['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
	return df

#master function that merges all dfs (NOTE ALERT IT DOES NOT WORK well for multiple samples, it aribtrarily drops them)
#takes a list of tuples of a df in the following format: (df, idCol, [colsToExpand])
#NOTE DEPRECATED BY COMBINE IMPACT DATA UTIL DOT PY
def amalgamate_all_dfs(dfTuples, masterIdColumnName, mode='return'):
	firstDfFlag = True
	amalgamatedDf = None
	for v in dfTuples:
		df, idCol, expandCols = v
		for col in expandCols:
			expand_col(df, col)
		print df.shape
		print '--------'
		df = extract_patient_name(df, patientNameColumn=idCol)
		df = df.drop_duplicates(subset=[idCol])
		print df.shape
		df = df.rename(columns={idCol: masterIdColumnName})
		if firstDfFlag:
			amalgamatedDf = df
			firstDfFlag = False
		else:
			amalgamatedDf = amalgamatedDf.merge(df, how='left', left_on=masterIdColumnName, right_on=masterIdColumnName) #Do a left join so we dont lost unmatched values
	if mode == 'return':
		return amalgamatedDf
	elif mode == 'write':
		print 'writing file to amalgamatedDf.tsv'
		amalgamatedDf.to_csv('amalgamatedDf.tsv', index=False, sep='\t')
	else:
		print 'uh-oh, improper mode specified'

#simple df amalgamation
#if mode = return returns the df, if mode equals write, writes the df
#note breaks if ID columns dont match (reuqires client to fix id columns beforehand) 
def amalgamate_dfs_same_col(dfs, idColName, writePath=None, mode='return'):
	finalDf = dfs[0]
	for df in dfs[1:]:
		finalDf = finalDf.merge(df, how='left', left_on=idColName, right_on=idColName)
	if mode == 'return':
		return finalDf
	if mode == 'write':
		print 'writing file to ', writePath
		finalDf.to_csv(writePath, index=False, sep='\t')		
		return finalDf

#parses data into a sample based dataframe with one tumor per row
#assumes data comes in a format where the first line is a list of patients, and other lines are zero matricies
def parse_data_format_cases_on_line_one(dataPath, write=False):
	dfList = []
	df = pd.read_table(dataPath)
	patients = list(df.columns)[1:]
	patientDf = pd.DataFrame(patients, columns=['Tumor_Sample_Barcode'])
	dfList.append(patientDf)
	for index, row in df.iterrows():
		if index%100 == 0: print index
		l = list(row)
		rowAsDf = pd.DataFrame(l[1:], columns=[l[0]])
		dfList.append(rowAsDf)
	df_final = reduce(lambda left,right: pd.merge(left,right, left_index=True, right_index=True), dfList)
	if write:
		df_final.to_csv('dataCNAReshaped.tsv', index=False, sep='\t')
	return df_final

#parses maf mutation data into a matrix of tumors mapped to whether genes are mutated or not
#considerPathogenicity: consider features of the mutations (missense silent etc)
#if countNMutationsPerGene write counts of mutations to columns, otherwise write boolean there/not there
#write: write the output to a file
def parse_mut_data_into_gene_matrix(dataPath, countNMutationsPerGene=False, considerPathogenicity=False, write=False):
	df = pd.read_table(dataPath, skiprows=[0])
	print df
	l = list(df['Tumor_Sample_Barcode'])
	d = {}
	for index, row in df.iterrows():
		barcode = row['Tumor_Sample_Barcode']
		if barcode in d:
			curVal = d[barcode]
			curVal.append(row['Hugo_Symbol'])
			d[barcode] = curVal
		else:
			d[barcode] = [row['Hugo_Symbol']]
	l = [] #create a list of dictionaries so we can force this into a dataframe
	for key, val in d.items():
		if countNMutationsPerGene: val = set(val)
		counts = Counter(val)
		d = dict(counts)
		d['Tumor_Sample_Barcode'] = key
		l.append(d)
	df = pd.DataFrame(l)
	df = df.fillna(value=0)
	if write:
		df.to_csv('mutationDataGenesAsCols.tsv', index=False, sep='\t')
	return df





