#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main(path):
    dfTranscripts = pd.DataFrame()
    dfLoc= pd.DataFrame()
    df = pd.read_csv(path + 'data/homo_sapiens/pG4WT.csv', sep='\t')
    dfLoc = dfLoc.append(df[['Transcript', 'locStart', 'locEnd', 'Location', 'Biotype']])
    dfLoc = dfLoc.drop_duplicates(subset=None, keep='first', inplace=False)
    del df
    dfLoc = dfLoc[dfLoc.Location != 'acceptor']
    dfLoc = dfLoc[dfLoc.Location != 'CDS']
    dfLoc = dfLoc[dfLoc.Location != 'donor']
    dfLoc = dfLoc[dfLoc.Location != 'five_prime_utr']
    dfLoc = dfLoc[dfLoc.Location != 'junction']
    dfLoc = dfLoc[dfLoc.Location != 'start_codon']
    dfLoc = dfLoc[dfLoc.Location != 'stop_codon']
    dfLoc = dfLoc[dfLoc.Location != 'three_prime_utr']
    dfLoc = dfLoc[dfLoc.Biotype == 'lncRNA']
    dfLoc = dfLoc.drop_duplicates(subset=None, keep='first', inplace=False)
    print(len(dfLoc))
    groups = dfLoc.groupby('Transcript')
    for name, group in groups:
        trLength = group['locEnd'].max() - group['locStart'].min()
        row = {'Transcript': [name], 'Length': [trLength], 'nbpG4': [len(group)]}
        dfTranscripts = dfTranscripts.append(pd.DataFrame(data=row))
    dfTranscripts['Density'] = dfTranscripts.nbpG4 / dfTranscripts.Length *1000
    dfTranscripts = dfTranscripts.sort_values(['Density'], ascending=[False])
    print(dfTranscripts)
    dfTranscripts.to_csv(path_or_buf=path + 'data/homo_sapiens/TranscriptDensity.csv', header=True, index=None, sep='\t')
    sns.kdeplot(dfTranscripts['Density'])
    plt.show()
    sns.boxplot( y=dfTranscripts["Density"] )
    plt.show()
    sns.violinplot(y=dfTranscripts["Density"])
    plt.show()




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