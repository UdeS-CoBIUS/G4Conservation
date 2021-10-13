#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import os
import math
import argparse
from Bio import SeqIO
from pprint import pprint
from ushuffle import shuffle, Shuffler

def shuff(fileseq, path, type, sp):
    fastaOrigin = SeqIO.parse(open(fileseq),'fasta')
    for fasta in fastaOrigin:
        name, seq = fasta.id, str(fasta.seq)
        seq = bytes(seq, "utf8")

        shuffler = Shuffler(seq, 1)
        for repro in range(1, 51):
            seqres = shuffler.shuffle()
            seqres = seqres.decode("utf-8")

            outFilename = path+'reviewShuffle/'+sp+'/Repro'+str(repro)+'/Sequences_Shuffled_Mono_'+type+'.fa'
            output = open(outFilename, "a")
            output.write(">"+name + "\n")
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
        for repro in range(1, 51):
            seqres = shuffler.shuffle()
            seqres = seqres.decode("utf-8")

            outFilename = path+'reviewShuffle/'+sp+'/Repro'+str(repro)+'/Sequences_Shuffled_Tri_'+type+'.fa'
            output = open(outFilename, "a")
            output.write(">"+name + "\n")
            nbLine = math.ceil( float( len(seqres) ) / 60 )
            cpt1 = 0
            cpt2 = 60
            for i in range(0,int(nbLine)) :
                output.write(seqres[cpt1:cpt2] + "\n")
                # to have a new line after 60 characters
                cpt1 += 60
                cpt2 += 60
            output.close()

def main(path):
    spList = ['leishmania_major', 'mycobacterium_tuberculosis_h37rv',
    'myxococcus_xanthus_dk_1622', 'candidatus_korarchaeum_cryptofilum_opf8',
    'pyrococcus_horikoshii_ot3', 'chlamydomonas_reinhardtii',
    'saccharomyces_cerevisiae', 'escherichia_coli_str_k_12_substr_mg1655',
    'halobacterium_salinarum_r1', 'thermus_thermophilus_hb8']
    for sp in spList:
    # sp = 'chlamydomonas_reinhardtii'
        print(sp)
        fileGeneSeq = path + 'data/' +sp+ '/Sequences_Gene_WT.fa'
        # fileLocSeq = path + 'data/' +sp+ '/Sequences_WT.fa'

        shuffDico = {}
        shuffDico['Gene'] = shuff(fileGeneSeq, path, 'Gene', sp)

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
