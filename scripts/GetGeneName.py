#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import os
import argparse

def main(path):
    pG4 = path + 'data/homo_sapiens/TranscriptDensity.csv'
    geneNameF = '/home/anais/Downloads/mart_export.txt'
    dicoName = {}
    with open(geneNameF) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l:
                w = l.split('\t')
                dicoName[w[2]] = w[4]
    out = []
    with open(pG4) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l:
                w = l.split('\t')
                if w[0] in dicoName:
                    geneName = dicoName[w[0]]
                else:
                    geneName = ''
                out.append('\t'.join(w)+'\t'+geneName)
    output = open(path + 'data/homo_sapiens/TranscriptDensity_GeneName.csv', "w")
    output.write('\n'.join(out))
    output.close()



def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'analyseGC')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    main(path)