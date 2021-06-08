#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

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
    Same as generateSequences.py but for genes. The main differences are about
    the header.

"""

def addGene(attributes):
    """Apply function to retrieve the transcript id of a feature.

    :param attributes: all attributes from the gtf file, corresponds to the
        ninth column of the file.
    :type attributes: string

    :returns: transcript id.
    :rtype: string
    """
    attributes = attributes.split(';')
    geneID = attributes[0].split('"')[1]
    return geneID

def parseDF(df):
    """Parses a dataFrame to retrieve important data and compute all locations.

    We only kept some locations : exon, CDS, 5UTR, 3UTR and codons. From
    attributes we kept : gene and transcript id, Biotype of the transcript and
    the type of the transcript (coding or not). Then we compute the coords of
    other locations (introns, junctions) that does not exist. Codons coords are
    changed to fit with the length of point location (40 nucleotides upstream
    and downstream of the site). The case were this window would go over the
    transcript is taken in account and coords are made to not go over those of
    a transcript.

    :param df: contains all transcripts feature for a species but non parsed.
    :type df: dataFrame

    :returns: dfTmp, parsed df.
    :rtype: dataFrame
    """
    dfTmp = pd.DataFrame()
    dfTmp = dfTmp.append(df[ df.Location.str.contains('gene') ].dropna())
    dfTmp['Gene'] = dfTmp.Attributes.apply(addGene)
    return dfTmp

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
            output = open(outputDir + 'Sequences_Gene_' + type + '.fa', "w")
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
        output = open(outputDir + 'Sequences_Gene_' + opt + '.fa', "w")
        for id in fasta[opt]:
            output.write(id + "\n")
            nbLine = math.ceil( float( len(fasta[opt][id]) ) / 60 )
            cpt1 = 0
            cpt2 = 60
            for i in range(0,int(nbLine)) :
                output.write(fasta[opt][id][cpt1:cpt2] + "\n")
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
    """ Reverse complement a DNA sequence.

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

    :param directory: directory name containing all chromosome fasta.
    :type directory: string
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

def getOriginSeq(row, chrSeq):
    """Gets the sequences for Origin overlaping location.

    :param chrSeq: sequence of one chromosome.
    :type chrSeq: string
    :param dicoLoc: contains all locations.
    :type dicoLoc: dictionary

    :returns: seq, sequence of a location.
    :rtype: string
    """
    negSeq = chrSeq[- -row.Start - 1 :]
    # sequence before the origin of replication
    posSeq = chrSeq[ 0 : row.End ]
    # sequence after the origin of replication
    seq = negSeq + posSeq
    return seq

def getLocationSeq(row, chrSeq):
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
    seq = chrSeq[ row.Start -1 : row.End]
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
        if (chr != 'MT' and chr != 'Mt' and chr != 'MtDNA' and chr != 'pES100' and chr != 'Pt'):
            print(chr)
            dfChr = pd.DataFrame()
            for file in fastaFiles:
                if (chr + '.fa' == file):
                    fastaFile = pathFasta + file
            chrSeq = importFastaChromosome(fastaFile, chr)
            dfChr = dfChr.append(dfGTF[ dfGTF.Chromosome == chr])
            for index, row in dfChr.iterrows():
                locId = chr  + ':' + str(row.Start) + '~' + str(row.End) + ':' + row.Strand
                if row.Start > 0:
                    seq = getLocationSeq(row, chrSeq)
                else:
                    seq = getOriginSeq(row, chrSeq)
                if row.Strand == '-':
                    seq = reverseSequence(seq)
                if seq != '':
                    randomSeq = shuffleSeq(seq)
                    id = '>' + row.Gene +':' + locId
                    dicoFasta['WT'][id] = seq
                    dicoFasta['Shuffled'][id] = randomSeq
    writeFasta(dicoFasta, outputDir, opt)

def importGTFdf(filename, sp, outputDir, opt):
    """Imports a gtf file into a dataframe.

    Read the gtf file in a dataFrame, then change the columns names. NI goes for
    'Not Important' (for me), those are columms that will be deleted. I also
    retrieve only what I need from attributes and then delate it. Once the
    dataFrame is parsed (all location computed mainly), the dataFrame is
    returned to be saved as a csv file.

    :param filename: name of the gtf file.
    :type filename: string

    :returns: df, contains all transcripts feature for a species.
    :rtype: dataFrame
    """
    try:
        if sp in ['leishmania_major', 'gasterosteus_aculeatus', 'pongo_abelii', 'vitis_vinifera']:
            df = pd.read_csv(filename, sep='\t', index_col=0, skiprows=4)
        elif sp == 'drosophila_melanogaster':
            df = pd.read_csv(filename, sep='\t', index_col=0, skiprows=3)
        else:
            df = pd.read_csv(filename, sep='\t', index_col=0, skiprows=5)
    except:
        print("This file couldn't be converted in data frame : " + filename)
    else:
        # dataFrame with all windows from G4RNA Screener
        df.columns = ['Source', 'Location','Start','End',
                    'NI1', 'Strand', 'NI2', 'Attributes']
        df.drop(['NI1', 'NI2', 'Source'], axis=1, inplace=True)
        df['Chromosome'] = df.index
        df = parseDF(df)
        del df['Attributes']
        createFasta(df, fastaFile, outputDir, opt)

def build_arg_parser():
    GITDIR = os.getcwd()+'/'
    parser = argparse.ArgumentParser(description = 'Parser_gtf')
    # parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-p', '--path', default ='/home/vana2406/scratch/G4Conservation/')
    parser.add_argument ('-o', '--option', default = 'Both')
    parser.add_argument ('-sp', '--specie', default = 'yersinia_pestis_biovar_microtus_str_91001')
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    sp = arg.specie    # specie to parse
    path = arg.path
    opt = arg.option
    gtfFile = path + "data/" + sp + "/" + sp + ".gtf"
    fastaFile = path + 'data/' + sp + '/Fasta/'
    outputDir = path + 'data/' + sp + '/'
    importGTFdf(gtfFile, sp, outputDir, opt)
