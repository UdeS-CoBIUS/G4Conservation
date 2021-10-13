#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import os
import argparse
import pandas as pd
from pprint import pprint

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement of computation.

Date:
	Jully 2020

Description:
	This script reads all ouput files from G4RNA Screener under the name
	'Sequences_WT_xxxxx.csv'. Overlapping windows will be merged. Here are the
	columns in the output : Strand, Chromosome, locStart, locEnd, GeneID,
	Location, TranscriptID, meancGcC, meanG4H, meanG4NN, pG4Start, pG4End,
	G4Sequence.

"""

def mergeOverlappingSequences(dfTmp):
	"""Merge the sequences of overlaping windows.

	:param dfTmp: contain overlaping windows.
	:type dfTmp: dataFrame

	:returns: seq, sequence merged.
	:rtype: string
	"""
	dfTmp = dfTmp.sort_values(by=['wStart'])
	seq = str(dfTmp.seqG4.iloc[0])
	for w in range(1,len(dfTmp)):
		stepTmp = int(dfTmp.wStart.iloc[w] - dfTmp.wEnd.iloc[w-1])-1
		# convert to int elsewise it's a float
		wSeq = dfTmp.seqG4.iloc[w]
		seq += wSeq[-stepTmp:]
	return seq

def getInfo(df):
	"""Retrieves informations of a windows and parse it into a dictionary.

	As gene windows and junction windows are not formated the same way, this
	function aims to parse them into the same type of dictionary.

	:param df: contain all overlaping windows.
	:type df: dataFrame

	:returns: dico, contains all infromation for one window.
	:rtype: dictionary
	"""
	geneDesc = df.geneDesc.iloc[0]
	geneDescSplit = geneDesc.split(':')
	dico = {'Gene' : [geneDescSplit[0]],
			'meancGcC' : [df.cGcC.mean()],
			'meanG4H' : [df.G4H.mean()],
			'meanG4NN' : [df.G4NN.mean()],
			'pG4Start' : [min(df.wStart)],
			'pG4End' : [max(df.wEnd)]}
	if len(geneDescSplit) > 5:
		dico['Chromosome'] = [geneDescSplit[2]]
		dico['Strand'] = [geneDescSplit[4]]
	else:
		dico['Chromosome'] = [geneDescSplit[1]]
		dico['Strand'] = [geneDescSplit[3]]
	return dico

def mergeWindows(df):
	"""Merge overlaping windows.

	:param df: contain overlaping windows.
	:type df: dataFrame

	:returns: pG4, contains the pG4 which is the merge of overlaping windows.
	:rtype: dictionary
	"""
	pG4rSeq = mergeOverlappingSequences(df)
	if len(pG4rSeq) >= 20:
		dicoInfo = getInfo(df)
		pG4Start = dicoInfo['pG4Start'][0]
		pG4End = dicoInfo['pG4End'][0]
		pG4 = {}
		pG4 = dicoInfo
		pG4['pG4Start'] = [min(pG4Start, pG4End)]
		pG4['pG4End'] = [max(pG4Start, pG4End)]
		pG4['Sequence'] = [pG4rSeq]
		pG4['Description'] = [df.geneDesc.iloc[0]]

	else:
		pG4 = {}
	return pG4

def filterOnScores(dicoParam, dfWindows):
	"""Filter the windows based on thresholds.

	:param dicoParam: contains all parameters that were given to g4rna screener.
	:type dicoParam: dictionnary
	:param dfWindows: contains all windows of all genes from one specie.
	:type dfWindows: dataframe

	:returns: dfWindows, with only windows upper thresholds.
	:rtype: dataFrame
	"""
	dfWindows = dfWindows[ dfWindows.cGcC >= dicoParam["cGcC"] ].dropna()
	dfWindows = dfWindows[ dfWindows.G4H >= dicoParam["G4H"] ].dropna()
	dfWindows = dfWindows[ dfWindows.G4NN >= dicoParam["G4NN"] ].dropna()
	return dfWindows

def mergeG4(df, dicoParam):
	"""Browses all junction window to find those that are overlapping.

	Here we browse all junctions windows. We will only kept those that overlap
	the 100 nucleotid. Indeed, if the window over thresholds don't overlap this
	position, it only in a gene and not a junction.

	:param df: contains all windows.
	:type df: dataFrame
	:param dicoParam: contains all parameters that were given to g4rna screener.
	:type dicoParam: dictionnary

	:returns: dfpG4, contain all pG4 for that strand.
	:rtype: dataFrame
	"""
	dfTmp = pd.DataFrame()
	dfpG4 = pd.DataFrame()
	dfTmp = dfTmp.append(df[0:1]) # store the first window
	if len(df) == 1:
		pG4 = mergeWindows(dfTmp)
		dfTmp = pd.DataFrame.from_dict(pG4)
		dfpG4 = dfpG4.append(dfTmp)
	else:
		for w in range(1,len(df)): # w for window
			pG4 = mergeWindows(dfTmp)
			# browses all windows over thresholds, exept the first one
			if (df.geneDesc.iloc[w] == df.geneDesc.iloc[w-1] and
			  (df.wStart.iloc[w] >= df.wStart.iloc[w-1] and \
			  df.wStart.iloc[w] <= df.wEnd.iloc[w-1])):
				# if window overlap, add window at the current pG4
				dfTmp = dfTmp.append(df[w:w+1])
				if w == len(df)-1:
					pG4 = mergeWindows(dfTmp)
					dfTmp = pd.DataFrame.from_dict(pG4)
					dfpG4 = dfpG4.append(dfTmp)
			else: # new pG4
				pG4 = mergeWindows(dfTmp)
				dfTmp = pd.DataFrame.from_dict(pG4)
				dfpG4 = dfpG4.append(dfTmp)
				dfTmp = df.iloc[w:w+1]
				if w == len(df)-1 :
					pG4 = mergeWindows(dfTmp)
					dfTmp = pd.DataFrame.from_dict(pG4)
					dfpG4 = dfpG4.append(dfTmp)
	return dfpG4

def merge(filename, dicoParam, repro):
	dfpG42 = pd.DataFrame()
	try:
		dfWindows = pd.read_csv(filename, sep='\t', index_col=0)
	except:
		print("This file couldn't be converted in data frame : " + filename)
	else:
		# dataFrame with all windows from G4RNA Screener
		dfWindows.columns = ['geneDesc','cGcC',
							'G4H','seqG4','wStart',
							'wEnd', 'G4NN']
		dfWindows = filterOnScores(dicoParam, dfWindows)
		dfpG42 = dfpG42.append(mergeG4(dfWindows, dicoParam))
		dfpG42['Repro'] = repro
		return dfpG42

def main(dicoParam, directory, repro):
	dfpG4MonoGene = pd.DataFrame()
	dfpG4DiGene = pd.DataFrame()
	dfpG4TriGene = pd.DataFrame()

	for path, dirs, files in os.walk(directory):
		for file in files:
			if '_00' in file and '.csv' in file:
				inputfile = directory+'/CSVFile/'+file
				if  '_Mono_' in file and '_Gene_' in file:
					dfpG4MonoGene = dfpG4MonoGene.append(merge(inputfile, dicoParam, repro))
					dfpG4MonoGene = dfpG4MonoGene.reset_index(drop=True)
				elif '_Di_' in file and '_Gene_' in file:
					dfpG4DiGene = dfpG4DiGene.append(merge(inputfile, dicoParam, repro))
					dfpG4DiGene = dfpG4DiGene.reset_index(drop=True)
				elif '_Tri_' in file and '_Gene_' in file:
					dfpG4TriGene = dfpG4TriGene.append(merge(inputfile, dicoParam, repro))
					dfpG4TriGene = dfpG4TriGene.reset_index(drop=True)
	if len(dfpG4MonoGene) > 0:
		dfpG4MonoGene = dfpG4MonoGene.drop_duplicates(subset=None, keep='first', inplace=False)
		dfpG4MonoGene = dfpG4MonoGene.reset_index(drop=True)
	if len(dfpG4DiGene) > 0:
		dfpG4DiGene = dfpG4DiGene.drop_duplicates(subset=None, keep='first', inplace=False)
		dfpG4DiGene = dfpG4DiGene.reset_index(drop=True)
	if len(dfpG4TriGene) > 0:
		dfpG4TriGene = dfpG4TriGene.drop_duplicates(subset=None, keep='first', inplace=False)
		dfpG4TriGene = dfpG4TriGene.reset_index(drop=True)

	dfpG4MonoGene.to_csv(path_or_buf=directory+'/pG4_Shuffle_Mono_Gene.csv', header=True, index=None, sep='\t')
	dfpG4DiGene.to_csv(path_or_buf=directory+'/pG4_Shuffle_Di_Gene.csv', header=True, index=None, sep='\t')
	dfpG4TriGene.to_csv(path_or_buf=directory+'/pG4_Shuffle_Tri_Gene.csv', header=True, index=None, sep='\t')


def createDicoParam(arg):
	"""Retrieves arguments and put them in a dictionary.

	:param arg: contains all arguments given to the script, those are principaly
		parameters from G4RNA Screener.
	:type arg: arg_parser

	:returns: dicoParam, contains all arguments given to the script.
	:rtype: dictionary
	"""
	dicoParam = {"G4H" : float(arg.THRESHOLD_G4H),
				"cGcC" : float(arg.THRESHOLD_CGCC),
				"G4NN" : float(arg.THRESHOLD_G4NN),
				"windowLength" : int(arg.WINDOW),
				"step" : int(arg.STEP)}
	return dicoParam

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	parser.add_argument ('-p', '--path', default = '/home/vana2406/scratch/'+\
		'G4Conservation/reviewShuffle/')
	parser.add_argument ('-sp', '--specie', default = \
		'escherichia_coli_str_k_12_substr_mg1655')
	parser.add_argument ('-r', '--repro', default = '1')
	parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
	parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
	parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
	parser.add_argument ('-W', '--WINDOW', default = 60)
	parser.add_argument ('-S', '--STEP', default = 10)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	sp = arg.specie
	repro = arg.repro
	path = arg.path+sp+'/Repro'+repro+'/'
	print("specie : " + sp)
	dicoParam = createDicoParam(arg)
	main(dicoParam, path, repro)
	print("\tDone")
