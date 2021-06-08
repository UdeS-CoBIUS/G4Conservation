#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

"""

Copyright:
    Copyright Universite of Sherbrooke, departement of biochemistry and
    departement of computation.

Date:
    Jully 2020

Description:
    This script aim to parse a gtf file into a tsv file. It imports a gtf file,
    and computes all locations used in this study (intron, codon start and stop,
    junction). Once those locations are computed, all informations are exported
    as a csv to get a faster access for further scripts. Also, to be more
    efficient, this script is run for each chromosome of a specie and then
    agglamerate together.

"""

import re
import os
import argparse
import pandas as pd
from pprint import pprint

def createDicoType():
    """Creates a dictionnary with all sublclasses and classes of transcripts.

    :returns: dicoFam, contains all classes and subclasses.
    :rtype: dictionary
    """
    dicoFam = {'IG_C_gene' : 'Coding',
            'IG_D_gene' : 'Coding',
            'IG_J_gene' : 'Coding',
            'IG_LV_gene' : 'Coding',
            'IG_M_gene' : 'Coding',
            'IG_V_gene' : 'Coding',
            'IG_Z_gene' : 'Coding',
            'nonsense_mediated_decay' : 'Coding',
            'nontranslating_CDS' : 'Coding',
            'non_stop_decay' : 'Coding',
            'protein_coding' : 'Coding',
            'TR_C_gene' : 'Coding',
            'TR_D_gene' : 'Coding',
            'TR_gene' : 'Coding',
            'TR_J_gene' : 'Coding',
            'TR_V_gene' : 'Coding',
            'transcribed_unitary_pseudogene' : 'Pseudogene',
            'disrupted_domain' : 'Pseudogene',
            'IG_C_pseudogene' : 'Pseudogene',
            'IG_J_pseudogene' : 'Pseudogene',
            'IG_pseudogene' : 'Pseudogene',
            'IG_V_pseudogene' : 'Pseudogene',
            'processed_pseudogene' : 'Pseudogene',
            'pseudogene' : 'Pseudogene',
            'transcribed_processed_pseudogene' : 'Pseudogene',
            'transcribed_unprocessed_pseudogene' : 'Pseudogene',
            'translated_processed_pseudogene' : 'Pseudogene',
            'translated_unprocessed_pseudogene' : 'Pseudogene',
            'TR_J_pseudogene' : 'Pseudogene',
            'TR_V_pseudogene' : 'Pseudogene',
            'unitary_pseudogene' : 'Pseudogene',
            'unprocessed_pseudogene' : 'Pseudogene',
            'polymorphic_pseudogene' : 'Pseudogene',
            'macro_lncRNA' : 'LongNC',
            'bidirectional_promoter_lncRNA' : 'LongNC',
            'sense_intronic' : 'LongNC',
            '3prime_overlapping_ncRNA' : 'LongNC',
            'ambiguous_orf' : 'LongNC',
            'antisense' : 'LongNC',
            'antisense_RNA' : 'LongNC',
            'lincRNA' : 'LongNC',
            'ncrna_host''non_coding' : 'LongNC',
            'processed_transcript' : 'LongNC',
            'lncRNA' : 'LongNC',
            'retained_intron' : 'LongNC',
            'sense_overlapping' : 'LongNC',
            'vaultRNA' : 'ShortNC',
            'scaRNA' : 'ShortNC',
            'sRNA' : 'ShortNC',
            'miRNA' : 'ShortNC',
            'miRNA_pseudogene' : 'ShortNC',
            'misc_RNA' : 'ShortNC',
            'misc_RNA_pseudogene' : 'ShortNC',
            'Mt_rRNA' : 'ShortNC',
            'Mt_tRNA' : 'ShortNC',
            'Mt_tRNA_pseudogene' : 'ShortNC',
            'ncRNA' : 'ShortNC',
            'pre_miRNA' : 'ShortNC',
            'RNase_MRP_RNA' : 'ShortNC',
            'piRNA' : 'ShortNC',
            'RNase_P_RNA' : 'ShortNC',
            'rRNA' : 'ShortNC',
            'rRNA_pseudogene' : 'ShortNC',
            'scRNA' : 'ShortNC',
            'scRNA_pseudogene' : 'ShortNC',
            'snlRNA' : 'ShortNC',
            'snoRNA' : 'ShortNC',
            'snoRNA_pseudogene' : 'ShortNC',
            'snRNA' : 'ShortNC',
            'snRNA_pseudogene' : 'ShortNC',
            'SRP_RNA' : 'ShortNC',
            'tmRNA' : 'ShortNC',
            'tRNA' : 'ShortNC',
            'tRNA_pseudogene' : 'Pseudogene',
            'ribozyme' : 'ShortNC'}
    return dicoFam

def retrieveBiotypeFronAttributes(attributes, feature):
    """Gets the biotype from attributes.

    :param feature: transcript, gene, exon, or utr.
    :type feature: string
    :param attributes: last colons of gtf file, contains a lot of informations
        but is different depending on the feature.
    :type attributes: list

    :returns: biotype of the feature. The biotype of a gene can be different
        compared to the transcript biotype.
    :rtype: string
    """
    biotype = ''
    for attribute in attributes:
        if re.search(feature+"_biotype", attribute):
            biotype = attribute.split('"')[1]
    return biotype

def retrieveIdTrFronAttributes(attributes):
    """Gets the id transcript from attributes.

    :param attributes: last colons of gtf file, contains a lot of informations
        but is different depending on the feature.
    :type attributes: list

    :returns: idTr, transctipt identifier.
    :rtype: string
    """
    idTr = ''
    for attribute in attributes:
        if re.search('transcript_id', attribute):
            idTr = attribute.split('"')[1]
    return idTr

def addTranscript(attributes):
    """Apply function to retrieve the transcript id of a feature.

    :param attributes: all attributes from the gtf file, corresponds to the
        ninth column of the file.
    :type attributes: string

    :returns: transcript id.
    :rtype: string
    """
    attributes = attributes.split(';')
    return retrieveIdTrFronAttributes(attributes)

def addBiotype(attributes):
    """Apply function to retrieve the biotype of a transcript.

    :param attributes: all attributes from the gtf file, corresponds to the
        ninth column of the file.
    :type attributes: string
    :returns: biotype of the transcript.
    :rtype: string
    """
    attributes = attributes.split(';')
    return retrieveBiotypeFronAttributes(attributes, 'transcript')

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

def addTypeTr(tr, dicoClass):
    """Apply function to retrieve the type of a transcript depending on biotype.

    :param tr: transcript id.
    :type tr: string
    :param dicoClass: {transcript id : class}.
    :type dicoClass: dictionary

    :returns: type, coding or non coding.
    :rtype: string
    """
    return dicoClass[tr]

def addAcctepor(df,  dicoTr, start, end):
    """Creates a row for a new acceptor junction.

    Creates a dictionary with all information about the new acceptor. This dico
    will be parsed into a dataFrame row.

    :param df: df with gtf informations.
    :type df: dataFrame
    :param dicoTr: contains transcript informations to find if the acceptor
        coords are over the transcript.
    :type dicoTr: dictionary
    :param start: intron strat.
    :type start: int
    :param end: intron end.
    :type end: int

    :returns:row, contains all information about the new acceptor.
    :rtype: dictionary
    """
    if (end - 40 < dicoTr[ df.Transcript[0] ]['Start']):
        #if start Acceptor would be before the start of the tr
        startN = dicoTr[ df.Transcript[0] ]['Start']
    else:
        startN = end -40
    if (end + 40 > dicoTr[ df.Transcript[0] ]['End']):
        #if end Acceptor would be after the start of the tr
        endN = dicoTr[ df.Transcript[0] ]['End']
    else:
        endN = end + 40
    if df.Strand[0] == '+':
        row = {'Location' : ['acceptor'], 'Start' : startN,
        'End' : endN, 'Strand' : df.Strand[0],
        'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
        'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
        'Biotype' : df.Biotype[0], 'Class' : df.Class[0]}
    else:
        row = {'Location' : ['donor'], 'Start' : startN,
        'End' : endN, 'Strand' : df.Strand[0],
        'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
        'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
        'Biotype' : df.Biotype[0], 'Class' : df.Class[0]}
    return row

def addDonor(df,  dicoTr, start, end):
    """Creates a row for a new donor junction.

    Creates a dictionary with all information about the new donor. This dico
    will be parsed into a dataFrame row.

    :param df: df with gtf informations.
    :type df: dataFrame
    :param dicoTr: contains transcript informations to find if the donor coords
        are over the transcript.
    :type dicoTr: dictionary
    :param start: intron strat.
    :type start: int
    :param end: intron end.
    :type end: int

    :returns:row, contains all information about the new donor.
    :rtype: dictionary
    """
    if (start - 40 < dicoTr[ df.Transcript[0] ]['Start']):
        startN = dicoTr[ df.Transcript[0] ]['Start']
    else:
        startN = start - 40
    if (start + 40 > dicoTr[ df.Transcript[0] ]['End']):
        endN = dicoTr[ df.Transcript[0] ]['End']
    else:
        endN = start + 40
    if df.Strand[0] == '+':
        row = {'Location' : ['donor'], 'Start' : startN,
        'End' : endN, 'Strand' : df.Strand[0],
        'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
        'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
        'Biotype' : df.Biotype[0], 'Class' : df.Class[0]}
    else:
        row = {'Location' : ['acceptor'], 'Start' : startN,
        'End' : endN, 'Strand' : df.Strand[0],
        'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
        'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
        'Biotype' : df.Biotype[0], 'Class' : df.Class[0]}
    return row

def addJunction(df, dicoTr, start, end):
    """Creates a row for a new junction splice site.

    Creates a dictionary with all information about the new junction. This dico
    will be parsed into a dataFrame row.

    :param df: df with gtf informations.
    :type df: dataFrame
    :param dicoTr: contains transcript informations to find if the junction
        coords are over the transcript.
    :type dicoTr: dictionary
    :param start: intron strat.
    :type start: int
    :param end: intron end.
    :type end: int

    :returns:row, contains all information about the new junction splice site.
    :rtype: dictionary
    """
    if (start - 40 < dicoTr[ df.Transcript[0] ]['Start']):
        startN = dicoTr[ df.Transcript[0] ]['Start']
    else:
        startN = start
    if (end + 40 > dicoTr[ df.Transcript[0] ]['End']):
        endN = dicoTr[ df.Transcript[0] ]['End']
    else:
        endN = end
    row = {'Location' : ['junction'], 'Start' : str(startN)+'|'+str(start),
    'End' : str(end)+'|'+str(endN), 'Strand' : df.Strand[0],
    'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
    'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
    'Biotype' : df.Biotype[0], 'Class' : df.Class[0]}
    return row

def addIntron(df, start, end):
    """Creates a row for a new intron.

    Creates a dictionary with all information about the new intron. This dico
    will be parsed into a dataFrame row.

    :param df: df with gtf informations.
    :type df: dataFrame
    :param start: intron strat.
    :type start: int
    :param end: intron end.
    :type end: int

    :returns:row, contains all information about the new intron.
    :rtype: dictionary
    """
    row = {'Location' : ['intron'], 'Start' : start,
    'End' : end, 'Strand' : df.Strand[0],
    'Attributes' : df.Attributes[0], 'Chromosome' : df.Chromosome[0],
    'Transcript' : df.Transcript[0], 'Gene' : df.Gene[0],
    'Biotype' : df.Biotype[0], 'Class' : df.Class[0]}
    return row

def getPointLocation(df, dicoTr):
    """Computes all point location that were not in the native gtf.

    Introns are firstly computed, then from them other locations are found. One
    row is added for each new location added.

    :param df: df with gtf informations.
    :type df: dataFrame
    :param dicoTr: contains transcript informations to find if the new coords
        of codons are over the transcript.
    :type dicoTr: dictionary

    :returns:dfTmpCodon, contains all point locations (introns, junction,
        acceptor and donor).
    :rtype: dataFrame
    """
    dfIntron = pd.DataFrame()
    dfExon = pd.DataFrame()
    dfExon = dfExon.append(df[ df.Location == 'exon'])
    groups = dfExon.groupby('Transcript')
    for name, group in groups:
        if len(group) > 1:
            # if more then 1 exon, then there is an intron
            group = group.sort_values(by=['Start'])
            group = group.reset_index(drop=True)
            starts = group.Start
            ends = group.End
            cptIntron = 0
            while cptIntron < len(starts) - 1:
                start = ends[cptIntron] + 1
                end = starts[cptIntron + 1] - 1
                if start != end and start < end:
                    row = addIntron(group, start, end)
                    dfRow = pd.DataFrame.from_dict(row)
                    dfIntron = dfIntron.append(dfRow)
                    row = addJunction(group,  dicoTr, start, end)
                    dfRow = pd.DataFrame.from_dict(row)
                    dfIntron = dfIntron.append(dfRow)
                    row = addDonor(group,  dicoTr, start, end)
                    dfRow = pd.DataFrame.from_dict(row)
                    dfIntron = dfIntron.append(dfRow)
                    row = addAcctepor(group,  dicoTr, start, end)
                    dfRow = pd.DataFrame.from_dict(row)
                    dfIntron = dfIntron.append(dfRow)
                cptIntron += 1
    dfIntron = dfIntron.reset_index(drop=True)
    dfIntron['Coords'] = [ [dfIntron.Start[x], dfIntron.End[x] ] for x in range(0,len(dfIntron))]
    return(dfIntron)

def getExceptionCoords(row):
    """ Get coords of unusual start/stop informations.

    Some chromosome are circular and genes can overlap the origin of
    replication which give strange coords. Sometimes coords are negative, other
    times they are positiv but over the max length. To correct that, we deal
    with each case indepently because their jsut a few case and it would have
    been a waste of time to try to check all possibility and compute the rigth
    coords.

    :param row: row with the problem (location with a strange coords).
    :type row: dataFrame

    :returns: start, end, right coords.
    :rtype: int
    """
    if row.Transcript == 'AAC68473' and row.Start == 1041920:
        start = 1041920
        end = 1041962
    elif row.Transcript == 'AAC68473':
        start = 1134
        end = 1176
    elif row.Transcript == 'OE_7223F':
        start = -784
        end = 283
    elif row.Transcript == 'BAA31944' and row.Start == 1738396:
        start = 1738396
        end = 1738438
    elif row.Transcript == 'BAA31944':
        start = 185
        end = 217
    elif row.Transcript == 'BAA31943' and row.Start == 1738133:
        start = 1738133
        end = 1738175
    elif row.Transcript == 'BAA31943':
        start = 31
        end = 71
    elif row.Transcript == 'BAA31942':
        start = 1737907
        end = 1737949
    elif row.Transcript == 'AAK43339' and row.Start == 2991448:
        start = 2991448
        end = 2991490
    elif row.Transcript == 'AAK43339':
        start = 212
        end = 252
    elif row.Transcript == 'AAR38856':
        start = row.End
        end = row.Start
    else:
        print(row)
    return start, end


def addCodon(df, dicoTr):
    """Adds codons location and change their coords.

    Codons are added appart from other location basiquely present in the gtf
    to easly modified their coords. Codons coords are changed to fit with the
    length of point location (40 nucleotides upstream and downstream of the
    site). The case were this window would go over the transcript is taken in
    account and coords are made to not go over those of a transcript.

    :param df: df with gtf informations.
    :type df: dataFrame
    :param dicoTr: contains transcript informations to find if the new coords
        of codons are over the transcript.
    :type dicoTr: dictionary

    :returns:dfTmpCodon, contains all codons locations.
    :rtype: dataFrame
    """
    dfTmpCodon = pd.DataFrame()
    dfTmpCodon = dfTmpCodon.append(df[ df.Location.str.contains('start_codon') ].dropna())
    dfTmpCodon = dfTmpCodon.append(df[ df.Location.str.contains('stop_codon') ].dropna())
    dfTmpCodon['Transcript'] = dfTmpCodon.Attributes.apply(addTranscript)
    dfTmpCodon = dfTmpCodon.reset_index(drop=True)
    for index, row in dfTmpCodon.iterrows():
        if dicoTr[row.Transcript]['Start'] < 0:
            start, end = getExceptionCoords(row)
            dfTmpCodon['Start'].iloc[index] = start
            dfTmpCodon['End'].iloc[index] = end
        else:
            # we want to extend +40 from each side of the codon but we need to
            # check also if the extend don't go over the transcript length
            if (row.Start - 40 < dicoTr[row.Transcript]['Start']):
                dfTmpCodon['Start'].iloc[index] = dicoTr[row.Transcript]['Start']
            else:
                dfTmpCodon['Start'].iloc[index] = row.Start - 40
            if (row.End + 40 > dicoTr[row.Transcript]['End']):
                dfTmpCodon['End'].iloc[index] = dicoTr[row.Transcript]['End']
            else:
                dfTmpCodon['End'].iloc[index] = row.End + 40
    return dfTmpCodon

def parseDF(df, dicoClass):
    """Parses a dataFrame to retrieve important data and compute all locations.

    We only kept some locations : exon, CDS, 5UTR, 3UTR and codons. From
    attributes we kept : gene and transcript id, transcript's biotype and
    the type of the transcript (coding or not). Then we compute the coords of
    other locations (introns, junctions) that does not exist. Codons coords are
    changed to fit with the length of point location (40 nucleotides upstream
    and downstream of the site). The case were this window would go over the
    transcript is taken in account and coords are made to not go over those of
    a transcript.

    :param df: contains all transcripts feature for a specie but non parsed.
    :type df: dataFrame
    :param dicoClass: {transcript id : class}.
    :type dicoClass: dictionary

    :returns: dfTmp, parsed df.
    :rtype: dataFrame
    """
    # create a special df to make a dictionary of it. This dictionary contains
    # informations for transcripts.
    dfTmpTr = pd.DataFrame()
    dfTmpTr = dfTmpTr.append(df[ df.Location.str.contains('transcript') ].dropna())
    dfTmpTr['Transcript'] = dfTmpTr.Attributes.apply(addTranscript)
    dicoTr = dfTmpTr.set_index('Transcript').to_dict('index')
    del dfTmpTr

    # make a new df to get only location we are interested
    dfTmp = pd.DataFrame()
    dfTmp = dfTmp.append(df[ df.Location.str.contains('exon') ].dropna())
    dfTmp = dfTmp.append(df[ df.Location.str.contains('CDS') ].dropna())
    dfTmp = dfTmp.append(df[ df.Location.str.contains('five_prime_utr') ].dropna())
    dfTmp = dfTmp.append(df[ df.Location.str.contains('three_prime_utr') ].dropna())
    dfTmpCodon = addCodon(df, dicoTr)
    dfTmp = dfTmp.append(dfTmpCodon.dropna())
    dfTmp = dfTmp.reset_index(drop=True)

    # add some columns to get information more easly
    dfTmp['Coords'] = [ [dfTmp.Start[x], dfTmp.End[x]] for x in range(0,len(dfTmp))]
    dfTmp['Transcript'] = dfTmp.Attributes.apply(addTranscript)
    dfTmp['Gene'] = dfTmp.Attributes.apply(addGene)
    dfTmp['index1'] = dfTmp.index
    dfTmp['Biotype'] = dfTmp.Attributes.apply(addBiotype)
    dfTmp['Class'] = dfTmp['Transcript'].apply(addTypeTr, args=(dicoClass, ))
    # dfTmp['Class'] = dfTmp.Tr.apply(addTypeTr)

    #compute missing location (intron and junctions)
    dfTmp = dfTmp.append( getPointLocation(dfTmp, dicoTr) )
    return dfTmp

def importGTFdf(filename, sp, chr, dicoClass):
    """Imports a gtf file into a dataframe.

    Read the gtf file in a dataFrame, then change the columns names. NI goes for
    'Not Important' (for me), those are columms that will be deleted. I also
    retrieve only what I need from attributes and then delate it. Once the
    dataFrame is parsed (all location computed mainly), the dataFrame is
    returned to be saved as a csv file.

    :param filename: name of the gtf file.
    :type filename: string
    :param sp: name of the species for which the script run.
    :type sp: string
    :param chr: name of the chr for which the script run.
    :type chr: string
    :param dicoClass: {transcript id : class}.
    :type dicoClass: dictionary

    :returns: df, contains all transcripts feature for a species.
    :rtype: dataFrame
    """
    try:
        #Firts lines of gtf file are not the same for each species so depending
        # on the species we can skip a certain number
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
        df['Chromosome'] = df['Chromosome'].astype(str)
        df['Transcript'] = df.Attributes.apply(addTranscript)

        #filter the df to keep only one chromosome
        df = df[df['Chromosome'] == chr]
        results = pd.DataFrame()
        results = parseDF(df, dicoClass)
        del results['Attributes']
        return results

def getTrClass(filename):
    """Read gtf file and get the class of all transcripts.

    Read the gtf file line by line to get all transcripts class. A lot of
    classes can be determined using the biotype, but for some short non coding
    transcripts, the length of the transcript is way higher than 200nt, so
    their class is changed to LongNC.

    :param filename: name of the gtf file.
    :type filename: string

    :returns: dicoCl, {transcript id : class}
    :rtype: dictionary
    """
    dicobt = createDicoType()
    dicoCl = {}
    with open(filename) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            if l:
                w = l.split('\t')
                if len(w)>3:
                    if w[2] == 'transcript':
                        bt = retrieveBiotypeFronAttributes(w[8].split(';'), 'transcript')
                        tr = retrieveIdTrFronAttributes(w[8].split(';'))
                        if bt in dicobt:
                            cl = dicobt[bt]
                            if cl == 'ShortNC' and (int(w[4]) - int(w[3])) > 200:
                                cl = 'LongNC'
                        else:
                            cl = 'None'
                        dicoCl[tr] = cl
    return dicoCl

def build_arg_parser():
    GITDIR = os.getcwd()+'/'
    parser = argparse.ArgumentParser(description = 'Parser_gtf')
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-chr', '--chromosome', default = '1')
    parser.add_argument ('-sp', '--specie', default = 'leishmania_major')
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    sp = arg.specie    # specie to parse
    path = arg.path
    chr = arg.chromosome
    print(path)
    gtfFile = path + "data/" + sp + "/" + sp + ".gtf"
    output = path + "data/" + sp + "/chrGTF/" + sp + "_chr" + str(chr) + ".csv"
    dicoClass = getTrClass(gtfFile)
    df = importGTFdf(gtfFile, sp, chr, dicoClass)
    print(df.columns)
    del df['index1']
    df.to_csv(path_or_buf=output, header=True, index=None, sep='\t')
