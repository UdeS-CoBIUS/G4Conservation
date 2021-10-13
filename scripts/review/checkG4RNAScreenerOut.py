import os

def getStat(directory, filetype):
    GPos = []
    GNeg = []
    CPos = []
    CNeg = []
    GCPos = []
    GCNeg = []
    for paths, dirs, files in os.walk(directory):
        for file in files:
            if filetype in file and '.csv' in file:
                filename = directory+file
                with open(filename) as f:
                    content = f.read()
                    lines = content.split('\n')
                    for l in lines:
                        w = l.split('\t')
                        if w[0] != '':
                            if float(w[2]) >= 4.5 and float(w[3]) >= 0.9 and float(w[-1]) >= 0.5:
                                GPos.append(w[4].count('G')/float(60))
                                CPos.append(w[4].count('C')/float(60))
                                GCPos.append((w[4].count('G') + w[4].count('C')) /float(60))
                            else:
                                GNeg.append(w[4].count('G')/float(60))
                                CNeg.append(w[4].count('C')/float(60))
                                GCNeg.append((w[4].count('G') + w[4].count('C')) /float(60))
    print('----------Positive window----------')
    print(len(GPos))
    print('G   ->    ', sum(GPos)/float(len(GPos)))
    print('GC  ->    ', sum(GCPos)/float(len(GCPos)))
    print('C   ->    ', sum(CPos)/float(len(CPos)))
    print('----------Negative window----------')
    print(len(GNeg))
    print('G   ->    ', sum(GNeg)/float(len(GNeg)))
    print('GC  ->    ', sum(GCNeg)/float(len(GCNeg)))
    print('C   ->    ', sum(CNeg)/float(len(CNeg)))

path = '/home/anais/Documents/Projet/G4Conservation/'

print('chlamydomonas_reinhardtii')
getStat(path+'data/chlamydomonas_reinhardtii/CSVFile/', 'Sequences_Gene_WT_00')
getStat(path+'reviewShuffle/chlamydomonas_reinhardtii/Repro9/CSVFile/', 'Sequences_Shuffled_Mono_Gene_00')
getStat(path+'reviewShuffle/chlamydomonas_reinhardtii/Repro9/CSVFile/', 'Sequences_Shuffled_Tri_Gene_00')

print('leishmania_major')
getStat(path+'data/leishmania_major/CSVFile/', 'Sequences_Gene_WT_00')
getStat(path+'reviewShuffle/leishmania_major/Repro9/CSVFile/', 'Sequences_Shuffled_Mono_Gene_00')
getStat(path+'reviewShuffle/leishmania_major/Repro9/CSVFile/', 'Sequences_Shuffled_Tri_Gene_00')

print('candidatus_korarchaeum_cryptofilum_opf8')
getStat(path+'data/candidatus_korarchaeum_cryptofilum_opf8/CSVFile/', 'Sequences_Gene_WT_00')
getStat(path+'reviewShuffle/candidatus_korarchaeum_cryptofilum_opf8/Repro9/CSVFile/', 'Sequences_Shuffled_Mono_Gene_00')
getStat(path+'reviewShuffle/candidatus_korarchaeum_cryptofilum_opf8/Repro9/CSVFile/', 'Sequences_Shuffled_Tri_Gene_00')

print('thermus_thermophilus_hb8')
getStat(path+'data/thermus_thermophilus_hb8/CSVFile/', 'Sequences_Gene_WT_00')
getStat(path+'reviewShuffle/thermus_thermophilus_hb8/Repro9/CSVFile/', 'Sequences_Shuffled_Mono_Gene_00')
getStat(path+'reviewShuffle/thermus_thermophilus_hb8/Repro9/CSVFile/', 'Sequences_Shuffled_Tri_Gene_00')

print('escherichia_coli_str_k_12_substr_mg1655')
getStat(path+'data/escherichia_coli_str_k_12_substr_mg1655/CSVFile/', 'Sequences_Gene_WT_00')
getStat(path+'reviewShuffle/escherichia_coli_str_k_12_substr_mg1655/Repro9/CSVFile/', 'Sequences_Shuffled_Mono_Gene_00')
getStat(path+'reviewShuffle/escherichia_coli_str_k_12_substr_mg1655/Repro9/CSVFile/', 'Sequences_Shuffled_Tri_Gene_00')
