configfile: "config.json"

rule all:
    input:
        expand("data/{specie}/Figures", specie = config["Species"])

rule download_gtf:
    output:
        "data/{specie}/{specie}.gtf"
    params:
        lambda wildcards: config["links"]["gtf"][wildcards.specie]
    shell:
        "wget -O {output}.gz {params} && gunzip {output}.gz"

rule download_fasta:
    output:
        directory("data/{specie}/Fasta/")
    params:
        link = lambda wildcards: config["links"]["fasta"][wildcards.specie],
        chrz = lambda wildcards: config["Chromosomes"][wildcards.specie]
    shell:
        "for chr in {params.chrz}; do wget -O {output}$chr.fa.gz \
        {params.link}$chr.fa.gz && gunzip {output}$chr.fa.gz; done"

rule parser_gtf:
    input:
        "data/{specie}/{specie}.gtf"
    output:
        directory("data/{specie}/chrGTF/")
    params:
        chrz = lambda wildcards: config["Chromosomes"][wildcards.specie]
    conda:
        "./envs/pypackages.yaml"
    shell:
        "mkdir --parents {output} &&\
        for chr in {params.chrz}; do sbatch scripts/Parser_gtf.sh \
        {wildcards.specie} $chr; done"

rule mergeGtf:
    input:
        "data/{specie}/chrGTF/"
    output:
        "data/{specie}/{specie}.csv"
    params:
        chrz = lambda wildcards: config["Chromosomes"][wildcards.specie]
    shell:
        "for chr in {params.chrz}; do cat {input}{wildcards.specie}_chr$chr.csv; done >> {output}"

rule generate_sequences:
    input:
        "data/{specie}/{specie}.csv", "data/{specie}/Fasta/"
    output:
        "data/{specie}/Sequences_WT.fa", "data/{specie}/Sequences_Shuffled.fa"
    shell:
        "python3 scripts/generateSequences.py -sp {wildcards.specie}"

rule geneSequences:
    input:
        "data/{specie}/{specie}.gtf", "data/{specie}/Fasta/"
    output:
        "data/{specie}/Sequences_Gene_WT.fa", "data/{specie}/Sequences_Gene_Shuffled.fa"
    shell:
        "python3 scripts/GeneSequences.py -sp {wildcards.specie}"

rule split_fasta:
    input:
        wt = "data/{specie}/Sequences_WT.fa",
        shuffle = "data/{specie}/Sequences_Shuffled.fa",
        genesWt = "data/{specie}/Sequences_Gene_WT.fa",
        genesShuffled = "data/{specie}/Sequences_Gene_Shuffled.fa"
    output:
        directory("data/{specie}/SplitFile")
    conda:
        "./envs/pypackages.yaml"
    shell:
        "mkdir --parents {output} && \
        python3 scripts/fasta_enumerator.py --fasta {input.wt} --ntLimit 5000000 \
        --large False --output {output} && \
        python3 scripts/fasta_enumerator.py --fasta {input.shuffle} \
        --ntLimit 5000000 --large False --output {output} && \
        python3 scripts/fasta_enumerator.py --fasta {input.genesWt} --ntLimit 5000000 \
        --large False --output {output} && \
        python3 scripts/fasta_enumerator.py --fasta {input.genesShuffled} --ntLimit 5000000 \
        --large False --output {output}"

rule G4RNA_screener:
    input:
        "data/{specie}/SplitFile"
    output:
        directory("data/{specie}/CSVFile")
    shell:
        """mkdir --parents {output} && \
        startTime=`date +%s` && \
        list_file=$(ls {input}) && \
        for file in $(ls {input} | cut -d '.' -f 1);
        do
            jobs_nb=$(squeue -u vana2406 | wc -l) && \
            if [ $jobs_nb -le 750 ]; then
                echo 'sbatch  ~/scratch/G4Conservation/scripts/G4RNAscreener.sh {input} {output} $file' && \
                sbatch ~/scratch/G4Conservation/scripts/G4RNAscreener.sh {input} {output} $file && \
                sleep 0.1
            else
                sleep 42
            fi;
        done"""

rule getLocations:
    input:
        "data/{specie}/Sequences_WT.fa"
    output:
        "data/{specie}/Locations.txt"
    shell:
        "cat {input} | grep '>' > {output}"

rule G4Annotation:
    input:
        "data/{specie}/CSVFile"
    output:
        "data/{specie}/pG4Shuffled.csv", "data/{specie}/pG4WT.csv"
    conda:
        "./envs/pypackages.yaml"
    shell:
        "python3 ./scripts/G4Annotation.py -sp {wildcards.specie}"

rule G4AnnotationGene:
    input:
        "data/{specie}/CSVFile"
    output:
        "data/{specie}/Gene_pG4Shuffled.csv", "data/{specie}/Gene_pG4WT.csv"
    conda:
        "./envs/pypackages.yaml"
    shell:
        "python3 ./scripts/G4AnnotationGene.py -sp {wildcards.specie}"

rule getMAinDensities:
    input:
        "data/{specie}/pG4Shuffled.csv", "data/{specie}/pG4WT.csv", "data/{specie}/Locations.txt"
    output:
        "data/{specie}/allDensitiesByLoc.csv"
    conda:
        "./envs/pypackages.yaml"
    shell:
        "python3 ./scripts/getMainDensities.py -sp {wildcards.specie}"

rule getDataByClass:
    input:
        "data/{specie}/pG4Shuffled.csv",
        "data/{specie}/pG4WT.csv",
        "data/{specie}/allDensitiesByLoc.csv",
        "data/{specie}/{specie}.csv"
    output:
        directory("data/{specie}/Figures")
    conda:
        "./envs/pypackages.yaml"
    shell:
        "mkdir --parents {output} && \
        python3 ./scripts/getDataFig.py -sp {wildcards.specie}"
