#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import re
import os
import sys
import math
from Bio import SeqIO
from pprint import pprint
from ushuffle import shuffle, Shuffler

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

def writeFasta(outFilename, seq, name):
    output = open(outFilename, "a")
    output.write(">"+name + "\n")
    nbLine = math.ceil( float( len(seq) ) / 60 )
    cpt1 = 0
    cpt2 = 60
    for i in range(0,int(nbLine)) :
        output.write(seq[cpt1:cpt2] + "\n")
        # to have a new line after 60 characters
        cpt1 += 60
        cpt2 += 60
    output.close()

def createFasta(path, sp, fileName):
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
    centroFile = path+sys.argv[1]+'/'+sys.argv[2]
    dicoInfoCentro = {}

    with open(centroFile) as f:
        content = f.read()
        lines = content.split('\n')
        header = lines[0]
        lines = lines[1:]
        for l in lines:
            w = l.split('\t')
            # print(w)
            if w[0] != '':
                if sp == 'saccharomyces_cerevisiae':
                    if 'centromere' in w[1]:
                        chr = w[7]
                        start = w[8]
                        end = w[9]
                        strand = w[10]
                        dicoInfoCentro[chr+':'+start+'-'+end+':'+strand] = {'Chromosome': chr,
                            'Start': start, 'End': end}
                elif sp == 'homo_sapiens':
                    chr = w[1]
                    start = w[2]
                    end = w[3]
                    dicoInfoCentro[chr+':'+start+'-'+end] = {'Chromosome': chr,
                        'Start': start, 'End': end}
                else:
                    # print(w)
                    if len(w) > 5:
                        if 'centromere' == w[7]:
                            chr = w[1]
                            start = w[2]
                            end = w[3]
                            dicoInfoCentro[chr+':'+start+'-'+end] = {'Chromosome': chr,
                                'Start': start, 'End': end}

    for centro in dicoInfoCentro:
        if 'chr' in dicoInfoCentro[centro]['Chromosome']:
            chr = dicoInfoCentro[centro]['Chromosome'].split('r')[1]
            chrSeq = importFastaChromosome(path+sp+'/Fasta/'+chr+'.fa', chr)
            seqForw = chrSeq[ int(dicoInfoCentro[centro]['Start']) - 1 : \
                int(dicoInfoCentro[centro]['End']) ]
            seqRev = reverseSequence(seqForw)

            writeFasta(path+sp+'/Sequences_centromere.fa', seqForw, centro+':+')
            writeFasta(path+sp+'/Sequences_centromere.fa', seqRev, centro+':-')

            seqForw = bytes(seqForw, "utf8")
            shuffler = Shuffler(seqForw, 1)
            for repro in range(1, 51):
                seqres = shuffler.shuffle()
                seqres = seqres.decode("utf-8")
                writeFasta(path+sp+'/Repro'+str(repro)+'/Sequences_centromere_Shuffled_Mono.fa',
                    seqres, centro+':+')

            shuffler = Shuffler(seqForw, 3)
            for repro in range(1, 51):
                seqres = shuffler.shuffle()
                seqres = seqres.decode("utf-8")
                writeFasta(path+sp+'/Repro'+str(repro)+'/Sequences_centromere_Shuffled_Tri.fa',
                    seqres, centro+':+')

            seqRev = bytes(seqRev, "utf8")
            shuffler = Shuffler(seqRev, 1)
            for repro in range(1, 51):
                seqres = shuffler.shuffle()
                seqres = seqres.decode("utf-8")
                writeFasta(path+sp+'/Repro'+str(repro)+'/Sequences_centromere_Shuffled_Mono.fa',
                    seqres, centro+':-')

            shuffler = Shuffler(seqRev, 3)
            for repro in range(1, 51):
                seqres = shuffler.shuffle()
                seqres = seqres.decode("utf-8")
                writeFasta(path+sp+'/Repro'+str(repro)+'/Sequences_centromere_Shuffled_Tri.fa',
                    seqres, centro+':-')

if __name__ == '__main__':
    path = '/home/anais/Documents/Projet/G4Conservation/reviewTRCentro/'
    createFasta(path, sys.argv[1], sys.argv[2])
