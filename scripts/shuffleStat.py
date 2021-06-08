#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import os
import argparse
import numpy as np
import pandas as pd

def getShortName():
    dicoSp = {'pan_troglodytes' : 'Ptro',
    'homo_sapiens' : 'Hsap',
    'pongo_abelii' : 'Pabe',
    'mus_musculus' : 'Mmus',
    'monodelphis_domestica' : 'Mdom',
    'ornithorhynchus_anatinus' : 'Oana',
    'anolis_carolinensis' : 'Acar',
    'gallus_gallus' : 'Ggal',
    'danio_rerio' : 'Drer',
    'gasterosteus_aculeatus' : 'Gacu',
    'drosophila_melanogaster' : 'Dmel',
    'apis_mellifera' : 'Amel',
    'caenorhabditis_elegans' : 'Cele',
    'neurospora_crassa' : 'Ncra',
    'aspergillus_nidulans' : 'Anid',
    'saccharomyces_cerevisiae' : 'Scer',
    'schizosaccharomyces_pombe' : 'Spom',
    'dictyostelium_discoideum' : 'Ddis',
    'arabidopsis_thaliana' : 'Atha',
    'vitis_vinifera' : 'Vvin',
    'solanum_lycopersicum' : 'Slyc',
    'oryza_sativa' : 'Osat',
    'physcomitrella_patens' : 'Ppat',
    'chlamydomonas_reinhardtii' : 'Crei',
    'leishmania_major' : 'Lmaj',
    'methanosarcina_acetivorans_c2a' : 'Mace',
    'halobacterium_salinarum_r1' : 'Hsal',
    'hyperthermus_butylicus_dsm_5456' : 'Hbut',
    'archaeoglobus_fulgidus_dsm_4304' : 'Aful',
    'methanobrevibacter_smithii_atcc_35061' : 'Msmi',
    'pyrococcus_horikoshii_ot3' : 'Phor',
    'thermoplasma_acidophilum_dsm_1728' : 'Taci',
    'sulfolobus_solfataricus_p2' : 'Ssol',
    'pyrobaculum_aerophilum_str_im2' : 'Paer',
    'nanoarchaeum_equitans_kin4_m' : 'Nequ',
    'candidatus_korarchaeum_cryptofilum_opf8' : 'Ckor',
    'cenarchaeum_symbiosum_a' : 'Csym',
    'aquifex_aeolicus_vf5' : 'Aaeo',
    'mycoplasma_pneumoniae_m129' : 'Mpne',
    'staphylococcus_aureus_subsp_aureus_n315' : 'Saur',
    'bacillus_subtilis_subsp_subtilis_str_168' : 'Bsub',
    'enterococcus_faecalis_v583' : 'Efae',
    'streptococcus_pneumoniae_tigr4' : 'Spne',
    'chloroflexus_aurantiacus_j_10_fl' : 'Caur',
    'mycobacterium_tuberculosis_h37rv' : 'Mtub',
    'thermus_thermophilus_hb8' : 'Tthe',
    'chlamydia_trachomatis_d_uw_3_cx' : 'Ctra',
    'borrelia_burgdorferi_b31' : 'Bbur',
    'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819' : 'Cjej',
    'myxococcus_xanthus_dk_1622' : 'Mxan',
    'geobacter_sulfurreducens_pca' : 'Gsul',
    'wolbachia_endosymbiont_of_drosophila_melanogaster' : 'Wend',
    'anaplasma_phagocytophilum_str_hz' : 'Apha',
    'brucella_abortus_bv_1_str_9_941' : 'Babo',
    'neisseria_meningitidis_z2491' : 'Nmen',
    'legionella_pneumophila_str_paris' : 'Lpne',
    'francisella_tularensis_subsp_tularensis_schu_s4' : 'Ftul',
    'vibrio_cholerae_o1_biovar_el_tor_str_n16961' : 'Vcho',
    'haemophilus_influenzae_rd_kw20' : 'Hinf',
    'yersinia_pestis_biovar_microtus_str_91001' : 'Ypes',
    'escherichia_coli_str_k_12_substr_mg1655' : 'Ecol'}
    return dicoSp

def applyShortName(sp):
    dico = getShortName()
    return dico[sp]

def appendDf(filename, sp, nbRepro):
    dftmp = pd.read_csv(filename, sep='\t')
    dftmp['Sp'] = sp
    dftmp['NbRepro'] = nbRepro
    return dftmp

def writeDf(df, outfile):
    dfTmp = pd.DataFrame()
    if 'Biotype' in df.columns and 'Location' in df.columns:
        grouped_multiple = df.groupby(['Location', 'Biotype', 'Sp'])
    elif 'Biotype' in df.columns:
        grouped_multiple = df.groupby(['Biotype', 'Sp'])
    else:
        grouped_multiple = df.groupby(['Class', 'Sp'])
    for name, group in grouped_multiple:
        avrg = group['DensityShuf'].mean()
        std = group['DensityShuf'].std()
        row = group[group.NbRepro == '3']
        row['meanShuffle'] = avrg
        row['stdShuffle'] = std
        dfTmp = dfTmp.append(row)
        # print(row)
        # row.to_csv(path_or_buf='Test', header=True, index=None, sep='\t')

    dfTmp['ShortName'] = dfTmp.Sp.apply(applyShortName)
    dfTmp.to_csv(path_or_buf=outfile, header=True, index=None, sep='\t')

def main(path):
    dfAllDens = pd.DataFrame()
    dfGlobal = pd.DataFrame() #figure 3
    dfCoding = pd.DataFrame() #figure 4
    dfCodingLoc = pd.DataFrame() #figure 5
    dfNC = pd.DataFrame() #figure 6
    dfNCLoc = pd.DataFrame() #figure 6
    dfPseudoLoc = pd.DataFrame() #figure 6
    dfPseudo = pd.DataFrame() #figure 6
    for path2, dirs, files in os.walk(path+'ReproShuffleG4Conservation/'):
        for repro in dirs:
            if 'Repro' in repro:
                print('-----------------------')
                print('|       '+repro+'        |')
                print('-----------------------')
                nbRepro = repro.split('o')[1]
                for path3, dirs2, files2 in os.walk(path+'ReproShuffleG4Conservation/'+repro):
                    for sp in dirs2:
                        if sp not in ['Figures', 'Fasta', 'CSVFile', 'SplitFile']:
                            # print(path + 'ReproShuffleG4Conservation/'+repro+'/'+sp + '/')
                            try :
                                dfAllDens = dfAllDens.append(appendDf(path + 'ReproShuffleG4Conservation/' + \
                                    repro+'/'+sp + '/allDensitiesByLoc.csv', sp, nbRepro))
                            except:
                                print(sp + ' doesnt have dfAllDens.')
                            try :
                                dfGlobal = dfGlobal.append(appendDf(path + 'ReproShuffleG4Conservation/'+  \
                                    repro+'/'+sp +'/Figures/DataFig3.csv', sp, nbRepro))
                            except:
                                print(sp + ' doesnt have dfGlobal.')
                            try :
                                dfCoding = dfCoding.append(appendDf(path +'ReproShuffleG4Conservation/' +  \
                                    repro+'/'+sp +'/Figures/DataFig4.csv', sp, nbRepro))
                            except:
                                print(sp + ' doesnt have dfCoding.')
                            try :
                                dfCodingLoc = dfCodingLoc.append(appendDf(path +'ReproShuffleG4Conservation/' + \
                                    repro+'/'+sp +'/Figures/DataFig5.csv', sp, nbRepro))
                            except:
                                print(sp + ' doesnt have dfCodingLoc.')
                            try :
                                dfNC = dfNC.append(appendDf(path +'ReproShuffleG4Conservation/' +  \
                                    repro+'/'+sp +'/Figures/DataFig6.csv', sp, nbRepro))
                            except:
                                print(sp + ' doesnt have non coding transcripts.')
                            try :
                                dfNCLoc = dfNCLoc.append(appendDf(path +'ReproShuffleG4Conservation/' +  \
                                    repro+'/'+sp +'/Figures/DataFig7.csv', sp, nbRepro))
                            except:
                                print(sp + ' doesnt have non coding transcripts.')
                            try :
                                dfPseudo = dfPseudo.append(appendDf(path +'ReproShuffleG4Conservation/' +  \
                                    repro+'/'+sp +'/Figures/DataSupFig1.csv', sp, nbRepro))
                            except:
                                print(sp + ' doesnt have pseudogene transcripts.')
                            try :
                                dfPseudoLoc = dfPseudoLoc.append(appendDf(path +'ReproShuffleG4Conservation/' +  \
                                    repro+'/'+sp +'/Figures/DataFig8.csv', sp, nbRepro))
                            except:
                                print(sp + ' doesnt have loc pseudogene transcripts.')
    dfCodingLoc['Class'] = 'Coding'
    dfGlobal['ShortName'] = dfGlobal.Sp.apply(applyShortName)
    dfNC['ShortName'] = dfNC.Sp.apply(applyShortName)
    writeDf(dfAllDens, 'Results/AllDensities.csv')
    writeDf(dfGlobal, 'Results/Global.csv')
    writeDf(dfCoding, 'Results/Coding.csv')
    writeDf(dfCodingLoc, 'Results/CodingLocation.csv')
    writeDf(dfNC, 'Results/NonCoding.csv')
    writeDf(dfNCLoc, 'Results/NonCodingLocation.csv')
    writeDf(dfPseudo, 'Results/Pseudo.csv')
    writeDf(dfPseudoLoc, 'Results/PseudoLoc.csv')


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
