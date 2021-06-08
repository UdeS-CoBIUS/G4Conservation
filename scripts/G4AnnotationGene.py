#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import os
import argparse
import Parser_gtf as pgtf
import pandas as pd
from pprint import pprint

"""

Copyright:
    Copyright Universite of Sherbrooke, departement of biochemistry and
    departement of computation.

Date:
    Jully 2020

Description:
    This script read all ouput files from G4RNA Screener under the name
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
    seq = str(dfTmp.seqG4.iloc[0])
    dfTmp = dfTmp.sort_values(by=['wStart'])
    for w in range(1,len(dfTmp)):
        step = int(dfTmp.wStart.iloc[w] - dfTmp.wStart.iloc[w-1])
        # convert to int elsewise it's a float
        wSeq = dfTmp.seqG4.iloc[w]
        seq += wSeq[-step:]
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
    if len(geneDescSplit) == 4:
        start = geneDescSplit[2].split('~')[0]
        end = geneDescSplit[2].split('~')[1]
        dico = {'Strand' : [geneDescSplit[3]],
                'Chromosome' : [geneDescSplit[1]],
                'locStart' : [int(float(start))], #for junction it's the
                'locEnd' : [int(float(end))], #intron start and end.
                'Gene' : [geneDescSplit[0]],
                'meancGcC' : [df.cGcC.mean()],
                'meanG4H' : [df.G4H.mean()],
                'meanG4NN' : [df.G4NN.mean()],
                'pG4Start' : [min(df.wStart)],
                'pG4End' : [max(df.wEnd)]}
    elif len(geneDescSplit) == 5:
        start = geneDescSplit[3].split('~')[0]
        end = geneDescSplit[3].split('~')[1]
        dico = {'Strand' : [geneDescSplit[4]],
                'Chromosome' : [geneDescSplit[2]],
                'locStart' : [int(float(start))], #for junction it's the
                'locEnd' : [int(float(end))], #intron start and end.
                'Gene' : [geneDescSplit[0]],
                'meancGcC' : [df.cGcC.mean()],
                'meanG4H' : [df.G4H.mean()],
                'meanG4NN' : [df.G4NN.mean()],
                'pG4Start' : [min(df.wStart)],
                'pG4End' : [max(df.wEnd)]}
    else:
        start = geneDescSplit[4].split('~')[0]
        end = geneDescSplit[4].split('~')[1]
        dico = {'Strand' : [geneDescSplit[5]],
                'Chromosome' : [geneDescSplit[3]],
                'locStart' : [int(float(start))], #for junction it's the
                'locEnd' : [int(float(end))], #intron start and end.
                'Gene' : [geneDescSplit[0]],
                'meancGcC' : [df.cGcC.mean()],
                'meanG4H' : [df.G4H.mean()],
                'meanG4NN' : [df.G4NN.mean()],
                'pG4Start' : [min(df.wStart)],
                'pG4End' : [max(df.wEnd)]}
    return dico

def mergeWindows(df, extension):
    """Merge overlaping windows.

    :param df: contain overlaping windows.
    :type df: dataFrame
    :param extension: length of unction, by default it's 100 nt.
    :type extension: integer

    :returns: pG4, contains the pG4 which is the merge of overlaping windows.
    :rtype: dictionary
    """
    pG4rSeq = mergeOverlappingSequences(df)
    dicoInfo = getInfo(df)
    if dicoInfo['Gene'][0] == 'CHLRE_01g049450v5':
        pprint(dicoInfo)
    if dicoInfo['Strand'][0] in ['1', '+']:
        pG4Start = dicoInfo['locStart'][0] + dicoInfo['pG4Start'][0] -1
        pG4End = dicoInfo['locStart'][0] + dicoInfo['pG4End'][0] -1
        pG4 = True
    else:
        pG4Start = dicoInfo['locEnd'][0] - dicoInfo['pG4Start'][0] +1
        pG4End = dicoInfo['locEnd'][0] - dicoInfo['pG4End'][0] +1
        pG4 = True
    if pG4:
        pG4 = {}
        pG4 = dicoInfo
        pG4['pG4Start'] = [min(pG4Start, pG4End)]
        pG4['pG4End'] = [max(pG4Start, pG4End)]
        pG4['Sequence'] = [pG4rSeq]
    if dicoInfo['Gene'][0] == 'CHLRE_01g049450v5':
        pprint(dicoInfo)
        pprint(pG4)
        print('-------')
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
        pG4 = mergeWindows(dfTmp, dicoParam["extension"])
        dfTmp = pd.DataFrame.from_dict(pG4)
        dfpG4 = dfpG4.append(dfTmp)
    else:
        for w in range(1,len(df)): # w for window
            # browses all windows over thresholds, exept the first one
            if (df.geneDesc.iloc[w] == df.geneDesc.iloc[w-1] and
              (df.wStart.iloc[w] >= df.wStart.iloc[w-1] and \
              df.wStart.iloc[w] <= df.wEnd.iloc[w-1])):
                # if window overlap, add window at the current pG4
                dfTmp = dfTmp.append(df[w:w+1])
                if w == len(df)-1 :
                    pG4 = mergeWindows(dfTmp, dicoParam["extension"])
                    dfTmp = pd.DataFrame.from_dict(pG4)
                    dfpG4 = dfpG4.append(dfTmp)
            else: # new pG4
                pG4 = mergeWindows(dfTmp, dicoParam["extension"])
                dfTmp = pd.DataFrame.from_dict(pG4)
                dfpG4 = dfpG4.append(dfTmp)
                dfTmp = df.iloc[w:w+1]
                if w == len(df)-1 :
                    pG4 = mergeWindows(dfTmp, dicoParam["extension"])
                    dfTmp = pd.DataFrame.from_dict(pG4)
                    dfpG4 = dfpG4.append(dfTmp)
    return dfpG4

def merge(filename, dicoParam):
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
        return dfpG42

def main(dicoParam, mainPath):
    output1 = mainPath + '/Gene_pG4WT.csv'
    output2 = mainPath + '/Gene_pG4Shuffled.csv'
    directory = mainPath + 'CSVFile'
    dfpG4WT = pd.DataFrame()
    dfpG4Shuffle = pd.DataFrame()
    for path, dirs, files in os.walk(directory):
        for filename in files:
            inputfile = directory + '/' + filename
            if ('Sequences_Gene_WT' in filename):
                dfpG4WT = dfpG4WT.append(merge(inputfile, dicoParam))
                dfpG4WT = dfpG4WT.reset_index(drop=True)
            elif ('Sequences_Gene_Shuffled' in filename):
                dfpG4Shuffle = dfpG4Shuffle.append(merge(inputfile, dicoParam))
                dfpG4Shuffle = dfpG4Shuffle.reset_index(drop=True)
    if len(dfpG4WT) > 0:
        dfpG4WT = dfpG4WT.drop_duplicates(subset=None, keep='first', inplace=False)
        dfpG4WT = dfpG4WT.reset_index(drop=True)
    if len(dfpG4Shuffle) > 0:
        dfpG4Shuffle = dfpG4Shuffle.drop_duplicates(subset=None, keep='first', inplace=False)
        dfpG4Shuffle = dfpG4Shuffle.reset_index(drop=True)
    if dfpG4WT.shape[0] > 1:
        dfpG4WT.to_csv(path_or_buf=output1, header=True, index=None, sep='\t')
    if dfpG4Shuffle.shape[0] > 1:
        dfpG4Shuffle.to_csv(path_or_buf=output2, header=True, index=None, sep='\t')

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
                "extension" : int(arg.extension),
                "windowLength" : int(arg.WINDOW),
                "step" : int(arg.STEP)}
    return dicoParam

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'G4Annotation')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-sp', '--specie', default = \
        'yersinia_pestis_biovar_microtus_str_91001')
    parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
    parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
    parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
    parser.add_argument ('-E', '--extension', default = 40)
    parser.add_argument ('-W', '--WINDOW', default = 60)
    parser.add_argument ('-S', '--STEP', default = 10)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    sp = arg.specie
    path = arg.path + 'data/' + sp + '/'
    print("Specie : " + sp)
    dicoParam = createDicoParam(arg)
    main(dicoParam, path)
    print("\tDone")
