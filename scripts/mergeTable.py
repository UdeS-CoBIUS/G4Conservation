#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import os
import argparse
import pandas as pd

def appendDf(filename, sp):
    dftmp = pd.read_csv(filename, sep='\t')
    dftmp['Sp'] = sp
    return dftmp

def main(path):
    dfGene = pd.DataFrame()
    dfTranscript = pd.DataFrame()
    for path2, dirs, files in os.walk(path+'toDl/'):
        for sp in dirs:
            if sp not in ['Figures', 'Fasta', 'CSVFile', 'SplitFile']:
                try :
                    dfGene = dfGene.append(appendDf(path+'toDl/'+sp+'/Gene_pG4WT.csv', sp))
                except:
                    print(sp + ' doesnt have dfGene.')
                try :
                    dfTranscript = dfTranscript.append(appendDf(path + 'toDl/'+sp+'/pG4WT.csv', sp))
                except:
                    print(sp + ' doesnt have dfTranscript.')

    outfile = path+'GenepG4.csv'
    dfGene.to_csv(path_or_buf=outfile, header=True, index=None, sep='\t')
    outfile = path+'TranscriptpG4.csv'
    dfTranscript.to_csv(path_or_buf=outfile, header=True, index=None, sep='\t')


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
