#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import os
import math
import argparse
from Bio import SeqIO
from pprint import pprint
from ushuffle import shuffle, Shuffler

def shuff(dicoFasta, path):
    for fasta in dicoFasta:
        seq = bytes(dicoFasta[fasta], "utf8")

        shuffler = Shuffler(seq, 1)
        for repro in range(1, 11):
            seqres = shuffler.shuffle()
            seqres = seqres.decode("utf-8")

            outFilename = path+'reviewTRCentro/ManySpecies/Repro'+str(repro)+'/Sequences_Shuffled_Mono_TE.fa'
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

            outFilename = path+'reviewTRCentro/ManySpecies/Repro'+str(repro)+'/Sequences_Shuffled_Tri_TE.fa'
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

def deallWithBac(path):
    bacteriaFasta = path +'reviewTRCentro/ManySpecies/BacTE.fasta'
    bacteriaCSV = path +'reviewTRCentro/ManySpecies/BacTE.csv'

    dicoBac = {}
    dicoF = {}

    with open(bacteriaCSV) as f:
        content = f.read()
        lines = content.split('\n')
        lines = lines[1:]
        for l in lines:
            w = l.split('\t')
            if w[0] != '':
                sp = w[-2][:-1]
                id = w[0].split(' ')[0]
                type = w[3]
                dicoBac[id] = [sp, type]

    fastaOrigin = SeqIO.parse(open(bacteriaFasta),'fasta')
    for fasta in fastaOrigin:
        name, seq = fasta.id, str(fasta.seq)
        if name in dicoBac:
            # name = name+':'+dicoBac[name][0]+':'+dicoBac[name][1]
            print(name+'::'+dicoBac[name][1]+'\t'+dicoBac[name][0])
            name = name+'::'+dicoBac[name][1]
            dicoF[name] = seq
    return dicoF

def importGacuFasta(path):
    GacuFasta = path +'reviewTRCentro/ManySpecies/Gasterosteus aculeatus.txt'
    dicoF = {}

    fastaOrigin = SeqIO.parse(open(GacuFasta),'fasta')
    for fasta in fastaOrigin:
        name, seq = fasta.id, str(fasta.seq)
        dicoF[name] = seq
        print(name+'\tGacu')
    return dicoF

def importRestEuk(path):
    eukFasta = path +'reviewTRCentro/ManySpecies/trep-db_complete_Rel-19.fasta'
    dicoF = {}
    listSp = ['Atal', 'Cele', 'Dmel', 'Ecol', 'Hsap', 'Osat', 'Ppat', 'Scer']

    fastaOrigin = SeqIO.parse(open(eukFasta),'fasta')
    for fasta in fastaOrigin:
        name, seq = fasta.id, str(fasta.seq)
        if name.split('_')[1] in listSp:
            print(name+'\t'+name.split('_')[1])
            dicoF[name] = seq
            # print(name+'::'+dicoBac[name][1]+'\tGacu')
    return dicoF

def main(path):
    dicoFasta = {}
    dicoFasta = deallWithBac(path)
    dicoFasta.update(importGacuFasta(path))
    dicoFasta.update(importRestEuk(path))

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
