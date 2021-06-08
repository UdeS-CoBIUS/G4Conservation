#!/usr/bin/python
from Bio import SeqIO
import argparse
import re
import os

def splitter(fasta_file, output, limit, large_handling=False):
    """
    Splits a large fasta_file in sub-files created in a given directory.
    """
    file_ = open(fasta_file, 'r')
    file_count = 1
    outfile = open(output.rstrip("/")+"/%s_%05d.fa"%(
        fasta_file.split('/')[-1].split('.')[0],file_count),'w')
    nt_count = 0
    for seq in SeqIO.parse(fasta_file, 'fasta'):
        if large_handling == True and len(str(seq.seq)) >= int(limit):
            file_count += 1
            largefile = open(output.rstrip("/")+"/%s_%05d_XL.fa"%(
                fasta_file.split('/')[-1].split('.')[0],file_count),'w')
            largefile.write(">"+str(seq.description)+"\n"+"\n".join(
                str(seq.seq)[i:i+50]for i in range(0,len(seq.seq),50))+"\n")
            largefile.close()
        else:
            nt_count += len(str(seq.seq))
            outfile.write(">"+str(seq.description)+"\n"+"\n".join(
                str(seq.seq)[i:i+50]for i in range(0,len(seq.seq),50))+"\n")    
            if nt_count >= int(limit):
                outfile.close()
                file_count += 1
                nt_count = 0
                outfile = open(output.rstrip("/")+"/%s_%05d.fa"%(
                    fasta_file.split('/')[-1].split('.')[0],file_count),'w')
    outfile.close()

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'generateRandom')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-f', '--fasta')
    parser.add_argument ('-n', '--ntLimit')
    parser.add_argument ('-o', '--output')
    parser.add_argument ('-l', '--large')
    return parser

def main():
    """
    Handles arguments.
    """
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    fasta = path + arg.fasta
    ntLimit = arg.ntLimit
    output = path + arg.output
    large = arg.large
    splitter(fasta, output, ntLimit, large)

if __name__ == '__main__':
    main()
