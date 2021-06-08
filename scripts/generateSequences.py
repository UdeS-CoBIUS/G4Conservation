#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import re
import os
import math
import random
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
    This script generates a fasta with all possible locations from a speices csv
    file. Two output can be generated: a 'wt' fasta file corresponding to
    real sequences, and a shuffled fasta file with the exact same sequences as wt
    ones but with the nucleotides order shuffled.

"""

def writeFasta(fasta, outputDir, opt):
    """From a dictionary {id : seq}, write a fasta file.

    :param fasta: {id : seq}
    :type fasta: dictionary
    :param outputDir: name of the output directory.
    :type outputDir: string
    :param opt: Both/Shuff/Wt.
    :type opt: string
    """
    if opt == 'Both':
        for type in fasta:
            output = open(outputDir + 'Sequences_' + type + '.fa', "w")
            for id in fasta[type]:
                output.write(id + "\n")
                nbLine = math.ceil( float( len(fasta[type][id]) ) / 60 )
                cpt1 = 0
                cpt2 = 60
                for i in range(0,int(nbLine)) :
                    output.write(fasta[type][id][cpt1:cpt2] + "\n")
                    # to have a new line after 60 characters
                    cpt1 += 60
                    cpt2 += 60
            output.close()
    else:
        output = open(outputDir + 'Sequences_' + opt + '.fa', "w")
        for id in fasta[type]:
            output.write(id + "\n")
            nbLine = math.ceil( float( len(fasta[type][id]) ) / 60 )
            cpt1 = 0
            cpt2 = 60
            for i in range(0,int(nbLine)) :
                output.write(fasta[type][id][cpt1:cpt2] + "\n")
                # to have a new line after 60 characters
                cpt1 += 60
                cpt2 += 60
        output.close()

def shuffleSeq(seq):
    """Shuffle a sequence.

    This function aims to shuffle a fasta sequence to randomize its sequence.
    The fasta sequence is imported from a fasta file, then converted into a
    list (one element corrresponds to one nucleotide), the list is shuffled and
    tehn joined with nothing to recreate the sequence.

    :param seq: sequence to shuffle.
    :type seq: string

    :returns: seq, sequence shuffled
    :rtype: string
    """
    seq = list(seq)
    random.shuffle(seq)
    seq = ''.join(seq)
    return seq

def reverseSequence(Sequence):
    """Reverse complement a DNA sequence.

    :param Sequence: DNA sequence that will be reversed.
    :type Sequence: string

    :returns: Sequence, the initial DNA sequence but reverse complemented.
    :rtype: string
    """
    reverse = ""
    nucleotides = {'A' : 'T',
                    'T' : 'A',
                    'C' : 'G',
                    'G' : 'C'}
    for n in Sequence:
        if n in nucleotides:
            tmp = nucleotides[n]
        else :
            tmp = n # in some sequences there is many N or other letter
        reverse += tmp
    reverse = reverse[::-1]
    return reverse

def importFastaChromosome(filename, chr):
    """Import the fasta file containing an entire chromosome sequence.

    :param filename: file name containing of a chromosome fasta.
    :type filename: string
    :param chr: chromosome wanted.
    :type chr: string

    :returns: sequence, entire chromosome sequence.
    :rtype: string
    """
    with open(filename) as f: # file opening
        content = f.read()
        l = content.split('\n')
        if l[0].startswith('>'):
            header = l[0]
            sequence = "".join(l[1:])
    return sequence

def getDicoLocTrBtCl(dfChr):
    """Gets a dictionary of locations with their bt-cl and tr.

    :param dfChr: contains information of only one chromosome.
    :type dfChr: dataFrame

    :returns: dicoLocTrBtCl, {locaction Id : transcriptId - biotype - class}
    :rtype: dictionary
    """
    dicoLocTrBtCl = {}
    listTest = []
    listTest = dfChr['Chromosome'].astype(str) + '~' + \
        dfChr['Start'].astype(str) + '~' + \
        dfChr['End'].astype(str) + '~' + \
        dfChr['Strand'].astype(str) + '~' + \
        dfChr['Transcript'].astype(str) + '~' + \
        dfChr['Biotype'].astype(str) + '~' + \
        dfChr['Location'].astype(str) + '~' + \
        dfChr['Class'].astype(str)
    for loc in listTest:
        w = loc.split('~')
        locId = w[6] + ':' + w[0] + ':' + w[1] + '~' + w[2] + ':' + w[3]
        if locId not in dicoLocTrBtCl:
            dicoLocTrBtCl[locId] = []
        # add Tr + biotype + class
        dicoLocTrBtCl[locId].append(w[4] + '~' + w[5] + '~' + w[-1])
    return dicoLocTrBtCl

def getJunctionSeq(chrSeq, dicoLoc):
    """Gets the sequence of junction.

    :param chrSeq: sequence of one chromosome.
    :type chrSeq: string
    :param dicoLoc: contains all locations.
    :type dicoLoc: dictionary

    :returns: seq, sequence of a junction location.
    :rtype: string
    """
    start =  dicoLoc['Start'].split('|')
    end = dicoLoc['End'].split('|')
    start = list(map(int, start))
    end = list(map(int, end))
    if  start[0] == start[1] and end[0] == end[1]:
        seqUpstream = chrSeq[ start[0]-1 - 40 : start[0] ]
        seqDownstream = chrSeq[ end[0]+1 : end[0] + 40]
    elif start[0] != start[1] and end[0] == end[1]:
        seqUpstream = chrSeq[ start[0]-1 : start[1] ]
        seqDownstream = chrSeq[ end[0]+1 : end[0] + 40]
    elif start[0] == start[1] and end[0] != end[1]:
        seqUpstream = chrSeq[ start[0]-1 - 40 : start[0] ]
        seqDownstream = chrSeq[ end[0]+1 : end[1] ]
    else:
        seqUpstream = chrSeq[ start[0]-1 : start[1] ]
        seqDownstream = chrSeq[ end[0]+1 : end[1] ]
    seq = seqUpstream + seqDownstream
    return seq

def getOriginSeq(chrSeq, dicoLoc):
    """Gets the sequences for Origin overlaping location.

    :param chrSeq: sequence of one chromosome.
    :type chrSeq: string
    :param dicoLoc: contains all locations.
    :type dicoLoc: dictionary

    :returns: seq, sequence of a location.
    :rtype: string
    """
    negSeq = chrSeq[- -dicoLoc['Start'] - 1 :]
    # sequence before the origin of replication
    posSeq = chrSeq[ 0 : dicoLoc['End'] ]
    # sequence after the origin of replication
    seq = negSeq + posSeq
    return seq

def getLocationSeq(w, chrSeq, dicoLoc):
    """Gets the sequences of a location.

    :param w: row of a dataFrame that contains information about a location.
    :type w: list
    :param chrSeq: sequence of one chromosome.
    :type chrSeq: string
    :param dicoLoc: contains all locations.
    :type dicoLoc: dictionary

    :returns: seq, sequence of site location.
    :rtype: string
    """
    if w[4] == 'junction':
        seq = getJunctionSeq(chrSeq, dicoLoc)
    else:
        seq = chrSeq[ dicoLoc['Start'] - 1 : \
            dicoLoc['End'] ]
    return seq

def createUniqID(dfChr):
    """Gets the sequences of a location.

    :param dfChr: contains information of only one chromosome.
    :type dfChr: dataFrame

    :returns: colUniqID, combination of many column to get a uniq ID.
    :rtype: column
    """
    colUniqID = dfChr['Gene'].astype(str) + '~' + \
        dfChr['Start'].astype(str) + '~' + \
        dfChr['End'].astype(str) + '~' + \
        dfChr['Strand'].astype(str) + '~' + \
        dfChr['Location'].astype(str)
    return colUniqID

def getDicoLoc(w):
    """Creates a dictionary with information for one location.

    :param w: row of a dataFrame that contains information about a location.
    :type w: list

    :returns: d, contains coords, strand and gene id of a location.
    :rtype: dictionary
    """
    d = {'gene' : w[0],
        'Start' : w[1],
        'End' : w[2],
        'Strand' : w[3]}
    return d

def createFasta(dfGTF, pathFasta, outputDir, opt):
    """Creates a dictionary with information for one location.

    :param dfGTF: row of a dataFrame that contains information about a location.
    :type dfGTF: list
    :param pathFasta: path of the directory containing chrmosomes fasta files.
    :type pathFasta: string
    :param outputDir: name of the output directory.
    :type outputDir: string
    :param opt: Both/Shuff/Wt.
    :type opt: string

    :returns: d, contains coords, strand and gene id of a location.
    :rtype: dictionary
    """
    dicoFasta = { 'WT' : {}, 'Shuffled' : {} }
    fastaFiles = [f for f in os.listdir(pathFasta) if \
        os.path.isfile( os.path.join(pathFasta, f) )]
    dfGTF.Chromosome = dfGTF.Chromosome.apply(str)
    chrList = list( set(dfGTF.Chromosome) )
    for chr in chrList:
        chr = str(chr)
        print(chr)
        dfChr = pd.DataFrame()
        for file in fastaFiles:
            if (chr + '.fa' == file and chr != 'MT'):
                fastaFile = pathFasta + file
        chrSeq = importFastaChromosome(fastaFile, chr)
        dfChr = dfChr.append(dfGTF[ dfGTF.Chromosome == chr])
        dfChr['UniqID'] =  createUniqID(dfChr)
        dicoLocTrBt = getDicoLocTrBtCl(dfChr)
        Loc = list(set(dfChr['UniqID']))
        for loc in Loc:
            w = loc.split('~')
            dicoLoc = getDicoLoc(w)
            locId = w[4] + ':' + chr  + ':' + w[1] + '~' + w[2] + ':' + w[3]
            if dicoLoc['End'] != dicoLoc['Start']:
                if w[4] != 'junction':
                    dicoLoc['Start'] = int(float(dicoLoc['Start']))
                    dicoLoc['End'] = int(float(dicoLoc['End']))
                    start = int(float(dicoLoc['Start']))
                else:
                    start = int(float(dicoLoc['Start'].split('|')[0]))
                if start > 0:
                    seq = getLocationSeq(w, chrSeq, dicoLoc)
                else:
                    seq = getOriginSeq(chrSeq, dicoLoc)
                if dicoLoc['Strand'] == '-':
                    seq = reverseSequence(seq)
                if seq != '':
                    randomSeq = shuffleSeq(seq)
                    listTr = list( set( dicoLocTrBt[locId] ) )
                    if len(w[0].split(':')) > 1:
                        w[0] = w[0].split(':')[0]
                    id = '>' + w[0] + ':'+ locId +':'+ '|'.join(listTr)
                    dicoFasta['WT'][id] = seq
                    dicoFasta['Shuffled'][id] = randomSeq
    writeFasta(dicoFasta, outputDir, opt)

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'generateRandom')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-o', '--option', default = 'Both')
    parser.add_argument ('-sp', '--specie', default = 'yersinia_pestis_biovar_microtus_str_91001')
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    sp = arg.specie
    path = arg.path
    opt = arg.option
    fastaFile = path + 'data/' + sp + '/Fasta/'
    outputDir = path + 'data/' + sp + '/'
    gtfParsed = path + 'data/' + sp + '/' + sp + '.csv'
    try:
        df = pd.read_csv(gtfParsed, sep='\t', header=0)
    except:
        print("This file couldn't be converted in data frame : " + gtfParsed)
    else:
        df = df[ df.Biotype != 'Biotype' ]
        createFasta(df, fastaFile, outputDir, opt)
