#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import os
import argparse
import pandas as pd
from pprint import pprint
import matplotlib.pyplot as plt

def main(path):
    spList = ['anaplasma_phagocytophilum_str_hz', 'leishmania_major', 'anolis_carolinensis', 'methanobrevibacter_smithii_atcc_35061', 'apis_mellifera', 'methanosarcina_acetivorans_c2a', 'aquifex_aeolicus_vf5', 'monodelphis_domestica', 'arabidopsis_thaliana', 'mus_musculus', 'archaeoglobus_fulgidus_dsm_4304', 'mycobacterium_tuberculosis_h37rv', 'aspergillus_nidulans', 'mycoplasma_pneumoniae_m129', 'bacillus_subtilis_subsp_subtilis_str_168', 'myxococcus_xanthus_dk_1622', 'borrelia_burgdorferi_b31', 'nanoarchaeum_equitans_kin4_m', 'brucella_abortus_bv_1_str_9_941', 'neisseria_meningitidis_z2491', 'caenorhabditis_elegans', 'neurospora_crassa', 'campylobacter_jejuni_subsp_jejuni_nctc_11168_atcc_700819', 'ornithorhynchus_anatinus', 'candidatus_korarchaeum_cryptofilum_opf8', 'oryza_sativa', 'cenarchaeum_symbiosum_a', 'pan_troglodytes', 'chlamydia_trachomatis_d_uw_3_cx', 'physcomitrella_patens', 'chlamydomonas_reinhardtii', 'pongo_abelii', 'chloroflexus_aurantiacus_j_10_fl', 'pyrobaculum_aerophilum_str_im2', 'danio_rerio', 'pyrococcus_horikoshii_ot3', 'dictyostelium_discoideum', 'saccharomyces_cerevisiae', 'drosophila_melanogaster', 'schizosaccharomyces_pombe', 'enterococcus_faecalis_v583', 'solanum_lycopersicum', 'escherichia_coli_str_k_12_substr_mg1655', 'staphylococcus_aureus_subsp_aureus_n315', 'francisella_tularensis_subsp_tularensis_schu_s4', 'streptococcus_pneumoniae_tigr4', 'gallus_gallus', 'sulfolobus_solfataricus_p2', 'gasterosteus_aculeatus', 'thermoplasma_acidophilum_dsm_1728', 'geobacter_sulfurreducens_pca', 'thermus_thermophilus_hb8', 'haemophilus_influenzae_rd_kw20', 'vibrio_cholerae_o1_biovar_el_tor_str_n16961', 'halobacterium_salinarum_r1', 'vitis_vinifera', 'homo_sapiens', 'wolbachia_endosymbiont_of_drosophila_melanogaster', 'hyperthermus_butylicus_dsm_5456', 'yersinia_pestis_biovar_microtus_str_91001', 'legionella_pneumophila_str_paris']
    dicoLengthGene = {}
    for sp in spList:
        file = path + 'data/' + sp + '/' + sp + '.gtf'
        # try:
        with open(file) as f:
            content = f.read()
            lines = content.split('\n')
            for l in lines:
                if not l.startswith('#'):
                    w = l.split('\t')
                    if w[0] != '':
                        # print(w)
                        if w[2] == 'gene':
                            if sp not in dicoLengthGene:
                                dicoLengthGene[sp] = []
                            if int(w[4]) > 0 and int(w[3])>0:
                                dicoLengthGene[sp].append( int(w[4]) - int(w[3]) +1 )
        # except:
        #     print('this file is not on my computer : '+sp)
    data = []
    label = []
    for sp in dicoLengthGene:
        data.append(dicoLengthGene[sp])
        label.append(sp)
    plt.boxplot(data, labels=label)
    plt.xticks(rotation=90)
    plt.savefig(path+'temp.png')
    pprint(dicoLengthGene)

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
