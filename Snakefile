#/project/gbru_wheat2/fhb/conda/exomecluster_env

wildcards = glob_wildcards("data/raw/{sample}_interleaved.fastq.gz")
SAMPLES = wildcards.sample

#SAMPLES = [
#    "AGS2000-RALEIGH-SRWW_S1_L001", 
#    "AGS2000-RALEIGH-SRWW_S9_L002", 
#    "AGS2000_S1_L001", 
#    "AGS2000_S9_L002", 
#    "HILLIARD_S27_L001", 
#    "HILLIARD_S27_L002"
#]

REF_ACCESSION = config['REF_ACCESSION']
REF_FASTA = config['REF_FASTA']

rule all:
    input:
        expand("fastqc_out/raw/{sample}_interleaved_fastqc.html", sample = SAMPLES), 
        expand("data/trimmed/{sample}_interleaved.trimmed.fastq", sample = SAMPLES), 
        expand("fastqc_out/trimmed/{sample}_interleaved.trimmed_fastqc.html", sample = SAMPLES), 
        expand("{accession}.{int}.ht21", accession = REF_ACCESSION)

rule fastqc_raw:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        "fastqc_out/raw/{sample}_fastqc.html"
    params:
        outdir = "fastqc_out/raw"
    log:
        "logs/fastqc_raw/{sample}.log"
    shell:
        "fastqc {input} -o {params.outdir} 2> {log}"

rule fastp_trim:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        "data/trimmed/{sample}.trimmed.fastq", 
    params:
        json = "reports/fastqc_trimmed/{sample}.json", 
        html = "reports/fastqc_trimmed/{sample}.html"
    log:
        "logs/fastp/{sample}.log"
    shell:
        "fastp --interleaved_in -j {params.json} -h {params.html} "
        "-i {input} --stdout > {output} 2> {log}"

rule fastqc_trimmed:
    input:
        "data/trimmed/{sample}.fastq"
    output:
        "fastqc_out/trimmed/{sample}_fastqc.html"
    params:
        outdir = "fastqc_out/trimmed"
    log:
        "logs/fastqc_trimmed/{sample}.log"
    shell:
        "fastqc {input} -o {params.outdir} 2> {log}"

rule hisat2_index:
    input:
        fasta = REF_FASTA
    output:
        "{accession}.{int}.ht21"
    params:
        accession = REF_ACCESSION
    threads: 12
    shell:
        "hisat2-build -p {threads} {input} {params.accession}"
