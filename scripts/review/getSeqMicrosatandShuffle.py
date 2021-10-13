#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import os
import math
import argparse
from Bio import SeqIO
from pprint import pprint
# from ushuffle import shuffle, Shuffler

def shuff(dicoFasta, path, sp):
    for fasta in dicoFasta:
        seq = bytes(dicoFasta[fasta], "utf8")

        shuffler = Shuffler(seq, 1)
        for repro in range(1, 11):
            seqres = shuffler.shuffle()
            seqres = seqres.decode("utf-8")

            outFilename = path+'reviewTRCentro/'+sp+'/Repro'+str(repro)+'/Sequences_Shuffled_Mono_Micro.fa'
            output = open(outFilename, "a")
            output.write(">"+fasta + "\n")
            nbLine = math.ceil( float( len(seqres) ) / 60 )
            cpt1 = 0
            cpt2 = 60
            for i in range(0,int(nbLine)) :
                output.write(seqres[cpt1:cpt2] + "\n")
                # to have a new line after 60 characters
                cpt1 += 60
                cpt2 += 60
            output.close()

        shuffler = Shuffler(seq, 3)
        for repro in range(1, 11):
            seqres = shuffler.shuffle()
            seqres = seqres.decode("utf-8")

            outFilename = path+'reviewTRCentro/'+sp+'/Repro'+str(repro)+'/Sequences_Shuffled_Tri_Micro.fa'
            output = open(outFilename, "a")
            output.write(">"+fasta + "\n")
            nbLine = math.ceil( float( len(seqres) ) / 60 )
            cpt1 = 0
            cpt2 = 60
            for i in range(0,int(nbLine)) :
                output.write(seqres[cpt1:cpt2] + "\n")
                # to have a new line after 60 characters
                cpt1 += 60
                cpt2 += 60
            output.close()

def updateStat(tag, stat, seq):
    stat[tag]['A'] += seq.count('A')
    stat[tag]['T'] += seq.count('T')
    stat[tag]['G'] += seq.count('G')
    stat[tag]['C'] += seq.count('C')
    if 'nbSeq' in stat[tag]:
        stat[tag]['nbSeq'] += 1
    return stat

def importFile(path, sp, filename):
    microsatFile = path +'reviewTRCentro/'+sp+'/'+filename
    dicoF = {}
    stat = {'micro sat G': {'A': 0, 'T': 0, 'G': 0, 'C': 0},
        'micro sat no G': {'A': 0, 'T': 0, 'G': 0, 'C': 0},
        'tot': {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'nbSeq': 0}}
    with open(microsatFile) as f:
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            w = l.split('\t')
            if w[0] != '':
                chr = w[0]
                startMicrosat = w[1]
                endMicrosat = w[2]
                totLength = int(w[4])
                strand = w[5]
                nbRep = int(w[6])
                seqRep = w[7]
                annotation = w[12]

                if 'G' in seqRep and totLength >= 20:
                    seq = seqRep * nbRep
                    if len(seq) != totLength:
                        diff = totLength - len(seq)
                        seq = seq+seq[0:diff]
                        id = chr+':'+startMicrosat+'-'+endMicrosat+':'+strand+\
                            ':'+annotation
                        dicoF[id] = seq
                        stat = updateStat('micro sat G', stat, seq)
                        stat = updateStat('tot', stat, seq)
                    else:
                        stat = updateStat('micro sat G', stat, seq)
                        stat = updateStat('tot', stat, seq)
                elif totLength >= 20:
                    seq = seqRep * nbRep
                    if len(seq) != totLength:
                        diff = totLength - len(seq)
                        seq = seq+seq[0:diff]
                        stat = updateStat('tot', stat, seq)
                        stat = updateStat('micro sat no G', stat, seq)
                    else:
                        stat = updateStat('tot', stat, seq)
                        stat = updateStat('micro sat no G', stat, seq)
    pprint(stat)
    return dicoF

def main(path):
    spDico = {'anolis_carolinensis': 'MSDBG000014_perf.tsv',
        'arabidopsis_thaliana': 'MSDBG000832_perf.tsv',
        'caenorhabditis_elegans': 'MSDBG001685_perf.tsv',
        'danio_rerio': 'MSDBG000039_perf.tsv',
        'drosophila_melanogaster': 'MSDBG000009_perf.tsv',
        'gallus_gallus': 'MSDBG000006_perf.tsv',
        'gasterosteus_aculeatus': 'MSDBG019944_perf.tsv',
        'homo_sapiens':'MSDBG000001_perf.tsv',
        'monodelphis_domestica': 'MSDBG000070_perf.tsv',
        'mus_musculus': 'MSDBG000071_perf.tsv',
        'ornithorhynchus_anatinus':'MSDBG000079_perf.tsv',
        'pan_troglodytes': 'MSDBG000002_perf.tsv',
        'pongo_abelii': 'MSDBG000088_perf.tsv',
        'saccharomyces_cerevisiae': 'MSDBG000011_perf.tsv'}
    for sp in spDico:
        print(sp)
        dicoFasta = {}
        dicoFasta = importFile(path, sp, spDico[sp])

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'G4Annotation')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    main(path)
