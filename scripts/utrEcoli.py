#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import os
import math
import argparse
import pandas as pd
import Parser_gtf as pgtf
import generateSequences as gS
import getMainDensities as gMD
from pprint import pprint

def applyClass(Bt):
    dicoClass = pgtf.createDicoType()
    if Bt in dicoClass:
        classe = dicoClass[Bt]
    else:
        classe = 'None' #TEC Bt
    return classe

def applyByotipe(desc):
    w = desc.split(':')
    w = w[5].split('~')
    return w[1]

def applyEnd(desc):
    w = desc.split(':')
    w = w[3].split('~')
    return w[1]

def applyStart(desc):
    w = desc.split(':')
    w = w[3].split('~')
    return w[0]

def applyStrand(desc):
    w = desc.split(':')
    return w[4]

def applyStrand(desc):
    w = desc.split(':')
    return w[4]

def applyLoc(desc):
    w = desc.split(':')
    return w[1]

def applyChromosome(desc):
    w = desc.split(':')
    return w[2]

def applyTranscript(desc):
    w = desc.split(':')
    w = w[5].split('~')
    return w[0]

def applyGene(desc):
    w = desc.split(':')
    return w[0]

def writeFasta(fasta, outFasta):
    """
    """
    output = open(outFasta, "w")
    for id in fasta:
        output.write(id + "\n")
        nbLine = math.ceil( float( len(fasta[id]) ) / 60 )
        cpt1 = 0
        cpt2 = 60
        for i in range(0,int(nbLine)) :
            output.write(fasta[id][cpt1:cpt2] + "\n")
            # to have a new line after 60 characters
            cpt1 += 60
            cpt2 += 60
    output.close()

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

def getLocationSeq(chrSeq, coord):
    """
    """
    seq = chrSeq[ coord[0] - 1 : coord[1] -1]
    return seq

def createFasta(dfUTR, pathFasta, w5, w3):
    """
    """
    dicoFastaWt = {}
    dicoFastaShuf = {}
    fastaFiles = [f for f in os.listdir(pathFasta) if \
        os.path.isfile( os.path.join(pathFasta, f) )]
    dfUTR.Chromosome = dfUTR.Chromosome.apply(str)
    chrList = list( set(dfUTR.Chromosome) )
    for chr in chrList:
        chr = str(chr)
        print(chr)
        dfChr = pd.DataFrame()
        for file in fastaFiles:
            if (chr + '.fa' == file and chr != 'MT'):
                fastaFile = pathFasta + file
        chrSeq = importFastaChromosome(fastaFile, chr)
        dfChr = dfChr.append(dfUTR[ dfUTR.Chromosome == chr])

        groups = dfUTR.groupby('Transcript')
        for name, group in groups:
            if list(group.Chromosome)[0] == chr:
                if len(group) == 1:
                    lowerCDS = group
                    higherCDS = group
                else:
                    group = group.sort_values(by=['Start'])
                    lowerCDS = group.iloc[0:]
                    higherCDS = group.iloc[-1:]
                if list(lowerCDS.Strand)[0] == '+':
                    utr5 = [ int(float(list(lowerCDS.Start)[0])) -w5, int(float(list(lowerCDS.Start)[0])) ]
                    utr3 = [ int(float(list(higherCDS.End)[0])), int(float(list(higherCDS.End)[0])) +w3 ]
                elif list(lowerCDS.Strand)[0] == '-':
                    utr3 = [ int(float(list(lowerCDS.Start)[0])) -w3, int(float(list(lowerCDS.Start)[0])) ]
                    utr5 = [ int(float(list(higherCDS.End)[0])), int(float(list(higherCDS.End)[0])) +w5 ]
                if utr5[0] < 0:
                    utr5[0] = 1
                if utr3[0] < 0:
                    utr3[0] = 1
                utr5Seq = getLocationSeq(chrSeq, utr5)
                utr3Seq = getLocationSeq(chrSeq, utr3)
                utr5SeqShuf = gS.shuffleSeq(utr5Seq)
                utr3SeqShuf = gS.shuffleSeq(utr3Seq)
                if list(lowerCDS.Strand)[0] == '-':
                    utr5Seq = reverseSequence(utr5Seq)
                    utr3Seq = reverseSequence(utr3Seq)
                if utr5[0] != utr5[1] and (utr5[1] - utr5[0]) > 20 and utr5[1] <= len(chrSeq)  and utr5[0] <= len(chrSeq):
                    id5utr = '>' + list(lowerCDS.Gene)[0] + ':five_prime_utr:'+ list(lowerCDS.Chromosome)[0] +':'+ \
                        str(utr5[0]) +'~'+ str(utr5[1]) +':'+list(lowerCDS.Strand)[0]+':'+ list(lowerCDS.Transcript)[0] +'~'+ \
                        list(lowerCDS.Biotype)[0]
                    dicoFastaWt[id5utr] = utr5Seq
                    dicoFastaShuf[id5utr] = utr5SeqShuf
                if utr3[0] != utr3[1] and (utr3[1] - utr3[0]) > 20 and utr3[1] <= len(chrSeq) and utr3[0] <= len(chrSeq):
                    id3utr = '>' + list(lowerCDS.Gene)[0] + ':three_prime_utr:'+ list(lowerCDS.Chromosome)[0] +':'+ \
                        str(utr3[0]) +'~'+ str(utr3[1])  +':'+list(lowerCDS.Strand)[0]+':'+ list(lowerCDS.Transcript)[0] +'~'+ \
                        list(lowerCDS.Biotype)[0]
                    dicoFastaWt[id3utr] = utr3Seq
                    dicoFastaShuf[id3utr] = utr3SeqShuf
    return dicoFastaWt, dicoFastaShuf

def parseCsv(filename):
    df = pd.read_csv(filename, sep='\t', header=0)
    df['Gene'] = df.description.apply(applyGene)
    df['Transcript'] = df.description.apply(applyTranscript)
    df['Strand'] = df.description.apply(applyStrand)
    df['Chromosome'] = df.description.apply(applyChromosome)
    df['Strand'] = df.description.apply(applyStrand)
    df['locStart'] = df.description.apply(applyStart)
    df['locEnd'] = df.description.apply(applyEnd)
    df['Location'] = df.description.apply(applyLoc)
    df['Biotype'] = df.description.apply(applyByotipe)
    df['Class'] = df.Biotype.apply(applyClass)

    listCoord = []
    for index, row in df.iterrows():
        if row['Strand'] == '+':
            pG4Start = int(row['locStart']) + int(row['start']) - 1
            pG4End = int(row['locStart']) + int(row['end']) -1
            pG4 = True
        else:
            pG4End = int(row['locEnd']) - int(row['start']) - 1
            pG4Start = int(row['locEnd']) - int(row['end']) -1
            pG4 = True
        listCoord.append( (pG4Start, pG4End, row.sequence) )
    dftmp =  pd.DataFrame(listCoord, columns=['pG4Start','pG4End', 'sequence'])

    del df['description']
    del df['start']
    del df['end']
    del df['Unnamed: 0']
    df = pd.merge(df, dftmp, on='sequence')
    col = ['meancGcC', 'meanG4H', 'Sequence', 'meanG4NN', 'Gene', 'Transcript',
           'Strand', 'Chromosome', 'locStart', 'locEnd', 'Location',
           'Biotype', 'Class', 'pG4Start', 'pG4End']
    df.columns = col

    # print(df)
    col = ['Biotype', 'Chromosome', 'Class', 'Gene', 'Location',
           'Sequence', 'Strand', 'Transcript', 'locEnd', 'locStart', 'meanG4H',
           'meanG4NN', 'meancGcC', 'pG4End', 'pG4Start']
    df = df.reindex(columns=col)
    return df

def stat(df, dfShuff, fastaFile):
    dfData = gMD.importData(fastaFile)
    dfData = gMD.sumNt(dfData)
    pG4Df = []
    groups = dfData.groupby('Location')
    for name, group in groups:
        print(name)
        group = df[df.Location == name]
        groupShuf = dfShuff[dfShuff.Location == name]
        pG4Df.append( ( float(len(group)), name, len(set(group.Transcript)),  float(len(groupShuf)), len(set(groupShuf.Transcript))  ) )
    dftmp =  pd.DataFrame(pG4Df, columns=['nbpG4','Location', 'nbTrpG4', 'nbpG4Shuf', 'nbTrShuf'])
    dfData = pd.merge(dfData, dftmp, on='Location')
    dfData['Tot'] = dfData['nuclC'] + dfData['nuclT'] + dfData['nuclG'] + dfData['nuclA']
    dfData['DensityWT'] = dfData['nbpG4'] / dfData['Tot'] *1000
    dfData['DensityShuf'] = dfData['nbpG4Shuf'] / dfData['Tot'] *1000
    return dfData

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'generateRandom')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-w5', '--w5UTR', default = 100)
    parser.add_argument ('-w3', '--w3UTR', default = 300)
    parser.add_argument ('-sp', '--specie', default = 'escherichia_coli_str_k_12_substr_mg1655')
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    sp = arg.specie
    path = arg.path
    w5 = int(arg.w5UTR)
    w3 = int(arg.w3UTR)
    print(sp)
    fastaFile = path + 'data/' + sp + '/Fasta/'
    outFasta = path + 'data/' + sp + '/Sequences_UTR.fa'
    gtfParsed = path + 'data/' + sp + '/' + sp + '.csv'
    try:
        df = pd.read_csv(gtfParsed, sep='\t', header=0)
    except:
        print("This file couldn't be converted in data frame : " + gtfParsed)
    else:
        df = df[ df.Biotype != 'Biotype' ]
        df = df[ df.Location == 'CDS' ]

        utrFastaWT,  utrFastaShuf= createFasta(df, fastaFile, w5, w3)
        writeFasta(utrFastaWT, outFasta)
        writeFasta(utrFastaShuf, path + 'data/' + sp + '/Sequences_UTR_Shuf.fa')

        cmd = '~/software/g4rna_screener/screen.py '+outFasta+ \
            ' -a ~/software/g4rna_screener/G4RNA_2016-11-07.pkl -w 60 -s 10 -c description cGcC G4H G4NN sequence start end -e >'+\
            path + 'data/' + sp + '/Sequences_UTR.csv'
        os.system(cmd)
        cmd = '~/software/g4rna_screener/screen.py '+ path + 'data/' + sp + '/Sequences_UTR_Shuf.fa -a ~/software/g4rna_screener/G4RNA_2016-11-07.pkl -w 60 -s 10 -c description cGcC G4H G4NN sequence start end -e >'+\
            path + 'data/' + sp + '/Sequences_UTR_Shuf.csv'
        os.system(cmd)

        cmd = '~/software/g4rna_screener/merge.py '+path + 'data/' + sp + '/'+'Sequences_UTR.csv'+ \
            ' --cGcC 4.5 --G4H 0.9 --G4NN 0.5 -a mean >'+\
            path + 'data/' + sp + '/pG4_UTR.csv'
        os.system(cmd)
        cmd = '~/software/g4rna_screener/merge.py '+ path + 'data/' + sp + '/Sequences_UTR_Shuf.csv'+ \
            ' --cGcC 4.5 --G4H 0.9 --G4NN 0.5 -a mean >'+\
            path + 'data/' + sp + '/pG4_UTR_Shuf.csv'
        os.system(cmd)

        dfShuff = parseCsv(path + 'data/' + sp + '/pG4_UTR_Shuf.csv')
        dfShuff.to_csv(path_or_buf=path+'data/'+sp+'/pG4_UTR_Shuf.csv', header=True, index=None, sep='\t')
        df = parseCsv(path + 'data/' + sp + '/pG4_UTR.csv')
        df.to_csv(path_or_buf=path+'data/'+sp+'/pG4_UTR.csv', header=True, index=None, sep='\t')

        df = stat(df, dfShuff, outFasta)
        df.to_csv(path_or_buf=path+'data/'+sp+'/pG4_UTRStat.csv', header=True, index=None, sep='\t')
