#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	September 2019

Description:
	With the output file of getMainDensities.py, this script compute all other
		statistics needed. There is one tsv output by figure.

Command Line:
	For one chr : python ~/PATH/getDataFig.py
"""

import os
import argparse
import pandas as pd
from pprint import pprint

def removeBiotype(df):
	"""Removes from a data frame, informations of some biotype we don't want.

	A lot of biotypes are removed for mainly 2 reasons : the biotype doesn't
	have enough transcript (less than 50) or the biotype is present in 
	less than 10 different species.

	:param df: origin dataframe with all biotypes.
	:type df: DataFrame

	:returns: updated dataframe with only biotype of interest
	:rtype: DataFrame
	"""
	df = df[ df.Biotype != 'IG_C_gene']
	df = df[ df.Biotype != 'IG_D_gene']
	df = df[ df.Biotype != 'IG_J_gene']
	df = df[ df.Biotype != 'IG_V_gene']
	df = df[ df.Biotype != 'TR_C_gene']
	df = df[ df.Biotype != 'TR_D_gene']
	df = df[ df.Biotype != 'TR_J_gene']
	df = df[ df.Biotype != 'TR_V_gene']
	df = df[ df.Biotype != 'TR_J_pseudogene']
	df = df[ df.Biotype != 'IG_C_pseudogene']
	df = df[ df.Biotype != 'IG_J_pseudogene']
	df = df[ df.Biotype != 'IG_pseudogene']
	df = df[ df.Biotype != 'TR_V_pseudogene']
	df = df[ df.Biotype != 'IG_LV_gene']
	df = df[ df.Biotype != 'non_stop_decay']
	df = df[ df.Biotype != 'nontranslating_CDS']
	df = df[ df.Biotype != 'RNase_MRP_RNA']
	df = df[ df.Biotype != 'pre_miRNA']
	df = df[ df.Biotype != 'SRP_RNA']
	df = df[ df.Biotype != 'vaultRNA']
	df = df[ df.Biotype != 'scaRNA']
	df = df[ df.Biotype != 'lncRNA']
	df = df[ df.Biotype != 'lincRNA']
	df = df[ df.Biotype != 'sense_overlapping']
	df = df[ df.Biotype != 'sense_intronic']
	df = df[ df.Biotype != 'rRNA_pseudogene']
	df = df[ df.Biotype != 'retained_intron']
	df = df[ df.Biotype != 'scRNA']
	df = df[ df.Biotype != 'antisense']
	df = df[ df.Biotype != 'ribozyme']
	df = df[ df.Biotype != 'sRNA']
	df = df[ df.Biotype != 'antisense_RNA']
	df = df[ df.Biotype != 'tmRNA']
	df = df[ df.Biotype != 'RNase_P_RNA']
	df = df[ df.Biotype != 'processed_transcript']
	df = df[ df.Biotype != 'piRNA']
	df = df[ df.Biotype != 'translated_unprocessed_pseudogene']
	df = df[ df.Biotype != 'transcribed_unitary_pseudogene']
	df = df[ df.Biotype != 'IG_V_pseudogene']
	df = df[ df.Biotype != 'tRNA_pseudogene']
	df = df[ df.Biotype != 'translated_processed_pseudogene']
	df = df[ df.Biotype != 'polymorphic_pseudogene']
	df = df[ df.Biotype != 'transcribed_unprocessed_pseudogene']
	df = df[ df.Biotype != 'transcribed_processed_pseudogene']
	df = df[ df.Biotype != 'unitary_pseudogene']
	df = df[ df.Biotype != 'processed_pseudogene']
	df = df[ df.Biotype != 'unprocessed_pseudogene']
	return df

def sumSubTable(group, name):
	"""sum columns a groupBy.

	:param group: contains a groupBy information.
	:type group: DataFrame
	:param name: name of the groupBy.
	:type name: string

	:returns: row, contains global information of a group.
	:rtype: dictionnary
	"""
	NbpG4rShuf = sum(group.NbpG4rShuf)
	NbpG4rWt = sum(group.NbpG4rWt)
	NbTrpG4Wt = sum(group.NbTrpG4Wt)
	NbTrpG4Shuf = sum(group.NbTrpG4Shuf)
	nbTr = sum(group.nbTr)
	Tot = sum(group.Tot)
	G = sum(group.nuclG)
	C = sum(group.nuclC)
	T = sum(group.nuclT)
	A = sum(group.nuclA)
	DensityShuf = float(NbpG4rShuf)/Tot *1000
	DensityWt = float(NbpG4rWt)/Tot *1000
	row = {'nuclG' : [G],
			'nuclC' : [C],
			'nuclT' : [T],
			'nuclA' : [A],
			'nbTr' : [nbTr],
			'NbpG4rWt' : [NbpG4rWt],
			'NbpG4rShuf' : [NbpG4rShuf],
			'NbTrpG4Wt' : [NbTrpG4Wt],
			'NbTrpG4Shuf' : [NbTrpG4Shuf],
			'Tot' : [Tot]}
	return row

def computePercentage(df, type):
	"""Computes densities of pG4.

	:param df: contains number of pG4.
	:type df: DataFrame
	:param type: location type, either point or segment.
	:type type: string

	:returns: tmp, updated with pG4 densities.
	:rtype: DataFrame
	"""
	tmp = pd.DataFrame()
	tmp = tmp.append(df)
	locType = 'NbpG4rLoc'+type
	namePercent = 'Percent'+type
	tmp[namePercent] = 0.0
	for index, row in tmp.iterrows():
		locId = row.LocID
		if row[locType] != 0:
			d = float(row[locType])/row.NbLocation *100
		else:
			d = 0
		tmp.loc[ tmp['LocID'] == locId, [namePercent] ] = d
		# tmp[namePercent].iloc[index] = d
	return tmp

def computeDensity(df, type):
	"""Computes densities of pG4.

	:param df: contains number of pG4.
	:type df: DataFrame
	:param type: location type, either point or segment.
	:type type: string

	:returns: tmp, updated with pG4 densities.
	:rtype: DataFrame
	"""
	tmp = pd.DataFrame()
	tmp = tmp.append(df)
	tmp['GC'] = (tmp['nuclG'] + tmp['nuclC']) / tmp['Tot'] *100
	if type == 'Point':
		tmp['DensityWt'] = tmp['NbpG4rWt'] / tmp['NbLocation']
		tmp['DensityShuf'] = tmp['NbpG4rShuf'] / tmp['NbLocation']
	else:
		tmp['DensityWt'] = tmp['NbpG4rWt'] / tmp['Tot'] * 1000
		tmp['DensityShuf'] = tmp['NbpG4rShuf'] / tmp['Tot'] * 1000
	return tmp

def computePercent(dico):
	"""Computes percentage of transcript with pG4.

	:param dico: number of transcript with at least one pG4 and tot tr.
	:type dico: dictionnary

	:returns: dico, update with the percentage of transcript with pG4.
	:rtype: dictionnary
	"""
	if dico['nbTr'] == 0:
		dico['PercentWt'] = 0
		dico['PercentShuf'] = 0
	else:
		dico['PercentWt'] = float(dico['NbTrpG4Wt']) / float(dico['nbTr']) * 100
		dico['PercentShuf'] = float(dico['NbTrpG4Shuf']) / float(dico['nbTr']) * 100
	return dico

def getFigBysubClass(df, dfCsv, dfWT, dfShuffle, nameClass):
	"""Same as for figure 3 but for subclass level.

	This function aims to filter the input dataFrame to get information about
	pG4r at class and subclass level. The final dataFrame contains
	information to make lolipop plot for GC content, densities and number
	of transcript.
	The glaobal statistic are cuomputed for each location and each informations
	(GC, densities, percent).

	:param df: contains all informations for all biotype locations.
	:type df: dataFrame
	:param dfCsv: contains all locations.
	:type dfCsv: dataFrame
	:param dfWT: contains all WT pG4.
	:type dfWT: dataFrame
	:param dfShuffle: contains all shuffled pG4.
	:type dfShuffle: dataFrame
	:param nameClass: name of the class we want to get the data (Coding, LongNC,
		Pseudogene).
	:type nameClass: string

	:returns: classDf, contains data for a location figure.
	:rtype: dataFrame
	"""
	tmp = pd.DataFrame()
	tmp = tmp.append(df)
	dicoNbTrClass = getNbTrByFeature(dfCsv, dfWT, dfShuffle, 'Class')
	dicoNbTrBt = getNbTrByFeature(dfCsv, dfWT, dfShuffle, 'Biotype')
	del tmp['Type']
	classDf = pd.DataFrame()
	classDftmp = tmp[ tmp.Class == nameClass]
	groups = classDftmp.groupby('Biotype')
	for name, group in groups:
		groupFilter = group[ group.Location == 'intron' ]
		groupFilter = groupFilter.append( group[ group.Location == 'exon' ])
		row = sumSubTable(groupFilter, name)
		row['Biotype'] = name
		row['Class'] = nameClass
		if name not in dicoNbTrBt['Tot']:
			dicoNbTrBt['Tot'][name] = 0
		if name not in dicoNbTrBt['WT']:
			dicoNbTrBt['WT'][name] = 0
		if name not in dicoNbTrBt['Shuff']:
			dicoNbTrBt['Shuff'][name] = 0
		row['nbTr'] = dicoNbTrBt['Tot'][name]
		row['NbTrpG4Wt'] = dicoNbTrBt['WT'][name]
		row['NbTrpG4Shuf'] = dicoNbTrBt['Shuff'][name]
		row.update(computePercent(row))
		row = pd.DataFrame(row, index=[len(classDftmp)+1])
		classDf = classDf.append(row)
	row = {'Class' : nameClass,
			'Biotype' : nameClass,
			'nuclG' : sum(classDftmp.nuclG),
			'nuclC' : sum(classDftmp.nuclC),
			'nuclT' : sum(classDftmp.nuclT),
			'nuclA' : sum(classDftmp.nuclA),
			'nbTr' : dicoNbTrClass['Tot'][nameClass],
			'NbpG4rWt' : sum(classDftmp.NbpG4rWt),
			'NbpG4rShuf' : sum(classDftmp.NbpG4rShuf),
			'NbTrpG4Wt' : dicoNbTrClass['WT'][nameClass],
			'NbTrpG4Shuf' : dicoNbTrClass['Shuff'][nameClass],
			'Tot' : sum(classDftmp.Tot)}
	row.update(computePercent(row))
	row = pd.DataFrame(row, index=[len(classDf)+1])
	classDf = classDf.append(row)
	classDf = computeDensity(classDf, 'Segment')
	return classDf

def getAllLocNum(df, suf, type, dftype):
	"""Computes the number of locations.

	:param df: contains pG4 informations depending on the dataset.
	:type df: dataFrame
	:param suf: either pG4 or tr
	:type suf: string
	:param type: either class or biotype
	:type type: string
	:param dftype: either 'Shuff' or 'WT'
	:type dftype: string
	:param nameClass: name of the class we want to get the data.
	:type nameClass: string

	:returns: dicoNbLoc, contains data for a location figure.
	:rtype: dictionnary
	"""
	dicoNbLoc = {}
	df['LocUniqId'] = df['Transcript'].astype(str) + '~' + \
		        df[suf+'Start'].astype(str)+'~'+df[suf+'End'].astype(str) +'~'+\
		        df['Strand'].astype(str) +'~'+ df['Chromosome'].astype(str)
	tmp = pd.DataFrame()
	if dftype == 'csv':
		tmp = tmp.append(df[['LocUniqId', 'Biotype', 'Location', suf+'Start', suf+'End', 'Class']])
	else:
		tmp = tmp.append(df[['LocUniqId', 'Biotype', 'Location', 'locStart', 'locEnd', 'Class']])
	tmp = tmp.drop_duplicates(subset=None, keep='first', inplace=False)
	if type =='ClassBiotype':
		groups = tmp.groupby(['Biotype','Location', 'Class'])
	else:
		groups = tmp.groupby([type,'Location'])
	for name, group in groups:
		if type =='ClassBiotype':
			name = (name[0]+'-'+name[2], name[1])
		if name[0] not in dicoNbLoc:
			dicoNbLoc[ name[0] ] = {name[1] : []}
		if name[1] not in dicoNbLoc[ name[0] ]:
			dicoNbLoc[ name[0] ][ name[1] ] = []
		dicoNbLoc[ name[0] ][ name[1] ] += list(group.LocUniqId)
	for cl in dicoNbLoc:
		for loc in dicoNbLoc[cl]:
			dicoNbLoc[cl][loc] = len(set(dicoNbLoc[cl][loc]))
	return dicoNbLoc

def fillEmptyClassLoc(dico1, dico2):
	"""Add 0 when no pG4 exist in a location/transcript class.

	:param dico1: csv informations.
	:type dico1: dictionnary
	:param dico2: pG4 informations.
	:type di2co: dictionnary

	:returns: dico2, contains data for a location figure.
	:rtype: dictionnary
	"""
	for cl in dico1:
		if cl not in dico2:
			dico2[cl] = {}
		for loc in dico1[cl]:
			if loc not in dico2[cl]:
				dico2[cl][loc] = 0
	return dico2

def getFigDensityByLocation(df, dfCsv, dfWT, dfShuffle, nameClass):
	"""Gets the data frame for figure with pG4r densities by location.

	This function aims to filter the input dataFrame to get information about
	pG4r at location level. The final dataFrame contains information to make
	lolipop plot for GC content, densities and number of location figure for
	segmental and point locations
	The glaobal statistic are cuomputed for each location and each informations
	(GC, densities, percent).

	:param df: contains all informations for all biotype locations.
	:type df: dataFrame
	:param dfCsv: contains all locations.
	:type dfCsv: dataFrame
	:param dfWT: contains all WT pG4.
	:type dfWT: dataFrame
	:param dfShuffle: contains all shuffled pG4.
	:type dfShuffle: dataFrame
	:param nameClass: name of the class we want to get the data (Coding, LongNC,
		Pseudogene).
	:type nameClass: string

	:returns: classDf, contains data for a location figure.
	:rtype: dataFrame
	"""
	tmp =  pd.DataFrame()
	dicoNbLoc = {}
	dicoNbLoc['Tot'] = getAllLocNum(dfCsv, '', 'Class', 'csv')
	dicoNbLoc['WT'] = getAllLocNum(dfWT, 'pG4', 'Class', 'wt')
	dicoNbLoc['Shuff'] = getAllLocNum(dfShuffle, 'pG4', 'Class', 'shuff')
	dicoNbLoc['WT'] = fillEmptyClassLoc(dicoNbLoc['Tot'], dicoNbLoc['WT'])
	dicoNbLoc['Shuff'] = fillEmptyClassLoc(dicoNbLoc['Tot'], dicoNbLoc['Shuff'])
	classDf = df[ df.Class == nameClass]
	dicoNbLoc['Tot'] = {k.split('-')[1] if '-' in k else k:v for k,v in dicoNbLoc['Tot'].items()}
	dicoNbLoc['WT'] = {k.split('-')[1] if '-' in k else k:v for k,v in dicoNbLoc['WT'].items()}
	dicoNbLoc['Shuff'] = {k.split('-')[1] if '-' in k else k:v for k,v in dicoNbLoc['Shuff'].items()}
	if nameClass in ['Pseudogene', 'LongNC', 'ShortNC']:
		#Those class should not get those location
		classDf = classDf[ classDf.Location != 'StartCodon']
		classDf = classDf[ classDf.Location != 'StopCodon']
	classDf = classDf.drop(['Class'], axis = 1)
	groups = classDf.groupby('Location')
	for name, group in groups:
		row = sumSubTable(group, name)
		row['Biotype'] = nameClass
		row['LocID'] = name+'-'+nameClass
		row['Location'] = name
		row['NbLocation'] = dicoNbLoc['Tot'][nameClass][name]
		row['NbpG4rLocWt'] = dicoNbLoc['WT'][nameClass][name]
		row['NbpG4rLocShuf'] = dicoNbLoc['Shuff'][nameClass][name]
		row['PercentWt'] = float(dicoNbLoc['WT'][nameClass][name]) / float(dicoNbLoc['Tot'][nameClass][name]) * 100
		row['PercentShuf'] = float(dicoNbLoc['Shuff'][nameClass][name]) / float(dicoNbLoc['Tot'][nameClass][name]) * 100
		if name in ['donor', 'acceptor', 'junction', 'StartCodon', 'StopCodon']:
			row['Type'] = 'Point'
		else:
			row['Type'] = 'Segment'
		row = pd.DataFrame(row, index=[len(tmp)+1])
		tmp = tmp.append(row)
	tmp1 = computeDensity(tmp[ tmp.Type == 'Point'], 'Point')
	tmp2 = computeDensity(tmp[ tmp.Type == 'Segment'], 'Segment')
	classDf = classDf.append(tmp1)
	classDf = classDf.append(tmp2)
	classDf = computePercentage(classDf, 'Wt')
	classDf = computePercentage(classDf, 'Shuf')
	return classDf

def getGroupTrNb(df, feature, dftype):
	""" Computes densities for all class/biotype of a species.

	When looking at a class, only exon and intron are concider to avoid
	redondance.

	:param df: contains all pG4 of a feature and dataset type.
	:type df: DataFrame
	:param feature: either class or biotype.
	:type feature: string
	:param dftype: either 'Shuff' or 'WT'
	:type dftype: string

	:param dicoNbTrClass: contains number of transcript by class/biotype.
	:type dicoNbTrClass: dictionnary
	"""
	dicoNbTr = {'Global' : []}
	dfClass = df.groupby(feature)
	for name, group in dfClass:
		dicoNbTr[name] = len(set(group['Transcript']))
		dicoNbTr['Global'] += list(group['Transcript'])
	dicoNbTr['Global'] = len(set(dicoNbTr['Global']))
	return dicoNbTr

def getNbTrByFeature(dfCsv, dfWT, dfShuffle, feature):
	""" Computes densities for all class/biotype of a species.

	When looking at a class, only exon and intron are concider to avoid
	redondance.

	:param dfCsv: contains all locations information of a species.
	:type dfCsv: DataFrame
	:param dfWT: contains all pG4 found in WT dataset.
	:type dfWT: DataFrame
	:param dfShuffle: contains all pG4 found in Shuffled dataset.
	:type dfShuffle: DataFrame
	:param feature: either class or biotype.
	:type feature: string

	:param dicoNbTrClass: contains number of transcript by class/biotype
	:type dicoNbTrClass: dictionnary
	"""
	dicoNbTrClass = {}
	if feature == 'Class':
		tmpWt = pd.DataFrame()
		tmpWt = tmpWt.append(dfWT[dfWT.Location == 'exon'])
		tmpWt = tmpWt.append(dfWT[dfWT.Location == 'intron'])
		tmpShuf = pd.DataFrame()
		tmpShuf = tmpShuf.append(dfShuffle[dfShuffle.Location == 'exon'])
		tmpShuf = tmpShuf.append(dfShuffle[dfShuffle.Location == 'intron'])
		dicoNbTrClass['Tot'] = getGroupTrNb(dfCsv, feature, 'csv')
		dicoNbTrClass['WT'] = getGroupTrNb(tmpWt, feature, 'wt')
		dicoNbTrClass['Shuff'] = getGroupTrNb(tmpShuf, feature, 'shuf')
		for name in dicoNbTrClass['Tot']:
			if name not in dicoNbTrClass['WT']:
				dicoNbTrClass['WT'][name] = 0
			if name not in dicoNbTrClass['Shuff']:
				dicoNbTrClass['Shuff'][name] = 0
	else:
		dicoNbTrClass['Tot'] = getGroupTrNb(dfCsv, feature, 'csv')
		dicoNbTrClass['WT'] = getGroupTrNb(dfWT, feature, 'wt')
		dicoNbTrClass['Shuff'] = getGroupTrNb(dfShuffle, feature, 'shuf')
		for name in dicoNbTrClass['Tot']:
			if name not in dicoNbTrClass['WT']:
				dicoNbTrClass['WT'][name] = 0
			if name not in dicoNbTrClass['Shuff']:
				dicoNbTrClass['Shuff'][name] = 0
	return dicoNbTrClass

def getDensitiesByClass(df, dfCsv, dfWT, dfShuffle):
	""" Computes densities for all transcripts class of a species.

	When looking at a class, only exon and intron are concider to avoid
	redondance.

	:param df: contains densities of all loc-bt of a species.
	:type df: DataFrame
	:param dfCsv: contains all locations information of a species.
	:type dfCsv: DataFrame
	:param dfWT: contains all pG4 found in WT dataset.
	:type dfWT: DataFrame
	:param dfShuffle: contains all pG4 found in Shuffled dataset.
	:type dfShuffle: DataFrame

	:param Global: contains all transcripts classes of a species.
	:type Global: DataFrame
	"""
	tmp = pd.DataFrame()
	tmp = tmp.append(df[df.Location == 'exon'])
	tmp = tmp.append(df[df.Location == 'intron'])
	Global = pd.DataFrame()
	groups = tmp.groupby('Class')
	dicoNbTrClass = getNbTrByFeature(dfCsv, dfWT, dfShuffle, 'Class')
	# print(tmp)
	for name, group in groups:
		row = sumSubTable(group, name)
		row['Class'] = [name]
		# pprint(dicoNbTrClass['Tot'])
		row['nbTr'] = [dicoNbTrClass['Tot'][name]]
		row['NbTrpG4Wt'] = [dicoNbTrClass['WT'][name]]
		row['NbTrpG4Shuf'] = [dicoNbTrClass['Shuff'][name]]
		row = pd.DataFrame.from_dict(row)
		Global = Global.append(row)
	row = {'Class' : 'Global',
			'nuclG' : sum(Global.nuclG),
			'nuclC' : sum(Global.nuclC),
			'nuclT' : sum(Global.nuclT),
			'nuclA' : sum(Global.nuclA),
			'nbTr': [dicoNbTrClass['Tot']['Global']],
			'NbTrpG4Wt' : [dicoNbTrClass['WT']['Global']],
			'NbTrpG4Shuf': [dicoNbTrClass['Shuff']['Global']],
			'NbpG4rWt' : sum(Global.NbpG4rWt),
			'NbpG4rShuf' : sum(Global.NbpG4rShuf),
			'Tot' : sum(Global.Tot)}
	row = pd.DataFrame(row, index=[len(Global)+1])
	Global = Global.append(row)
	Global['PercentWt'] = Global['NbTrpG4Wt'] / Global['nbTr'] * 100
	Global['PercentShuf'] = Global['NbTrpG4Shuf'] / Global['nbTr'] * 100
	Global = computeDensity(Global, 'Segment')
	return Global

def main(path, sp):
	results = path + sp + '/allDensitiesByLoc.csv'
	print(sp)
	try:
		df = pd.read_csv(results, sep='\t', index_col=0)
		dfCsv = pd.read_csv(path + sp + '/' + sp + '.csv', sep='\t')
		dfWT = pd.read_csv(path + sp + '/' + 'pG4WT.csv', sep='\t')
		dfShuffle = pd.read_csv(path + sp + '/' + 'pG4Shuffled.csv', sep='\t')
	except:
		print("This file couldn't be converted in data frame : " + results)
	else:
		df = removeBiotype(df)
		df = df.reset_index()
		dfCsv = removeBiotype(dfCsv)
		dfCsv = dfCsv.reset_index()
		dfWT = removeBiotype(dfWT)
		dfWT = dfWT.reset_index()
		dfShuffle = removeBiotype(dfShuffle)
		dfShuffle = dfShuffle.reset_index()
		fig3 = getDensitiesByClass(df, dfCsv, dfWT, dfShuffle)
		fig3.to_csv(path_or_buf=path+sp+'/Figures/DataFig3.csv', header=True, index=None, sep='\t')
		fig4 = getFigBysubClass(df, dfCsv, dfWT, dfShuffle, 'Coding')
		fig4.to_csv(path_or_buf=path+sp+'/Figures/DataFig4.csv', header=True, index=None, sep='\t')
		fig5 = getFigDensityByLocation(df, dfCsv, dfWT, dfShuffle, 'Coding')
		fig5.to_csv(path_or_buf=path+sp+'/Figures/DataFig5.csv', header=True, index=None, sep='\t')
		print(set(df.Class))
		if 'LongNC' in list(df['Class']) and 'ShortNC' in list(df['Class']):
			fig6 = getFigBysubClass(df, dfCsv, dfWT, dfShuffle, 'LongNC')
			fig6 = fig6.append(getFigBysubClass(df, dfCsv, dfWT, dfShuffle, 'ShortNC'))
			fig6.to_csv(path_or_buf=path+sp+'/Figures/DataFig6.csv', header=True, index=None, sep='\t')
		elif 'LongNC' in list(df['Class']) and 'ShortNC' not in list(df['Class']):
			fig6 = getFigBysubClass(df, dfCsv, dfWT, dfShuffle, 'LongNC')
			fig6.to_csv(path_or_buf=path+sp+'/Figures/DataFig6.csv', header=True, index=None, sep='\t')
		elif 'LongNC' not in list(df['Class']) and 'ShortNC' in list(df['Class']):
			fig6 = getFigBysubClass(df, dfCsv, dfWT, dfShuffle, 'ShortNC')
			fig6.to_csv(path_or_buf=path+sp+'/Figures/DataFig6.csv', header=True, index=None, sep='\t')
		if 'LongNC' in list(df['Class']):
			fig7 = getFigDensityByLocation(df, dfCsv, dfWT, dfShuffle, 'LongNC')
			fig7.to_csv(path_or_buf=path+sp+'/Figures/DataFig7.csv', header=True, index=None, sep='\t')
		if 'Pseudogene' in list(df['Class']):
			fig8 = getFigDensityByLocation(df, dfCsv, dfWT, dfShuffle, 'Pseudogene')
			fig8.to_csv(path_or_buf=path+sp+'/Figures/DataFig8.csv', header=True, index=None, sep='\t')
			supfig1 = getFigBysubClass(df, dfCsv, dfWT, dfShuffle, 'Pseudogene')
			supfig1.to_csv(path_or_buf=path+sp+'/Figures/DataSupFig1.csv', header=True, index=None, sep='\t')

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'getDataFig')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	parser.add_argument ('-sp', '--specie', default = 'yersinia_pestis_biovar_microtus_str_91001')
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path + 'data/'
	# path = arg.path + 'data/'
	sp = arg.specie
	main(path, sp)
