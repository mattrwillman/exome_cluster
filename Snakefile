#/project/gbru_wheat2/fhb/conda/exomecluster_env

SAMPLES = [
    "AGS2000-RALEIGH-SRWW_S1_L001", 
    "AGS2000-RALEIGH-SRWW_S9_L002", 
    "AGS2000_S1_L001", 
    "AGS2000_S9_L002", 
    "HILLIARD_S27_L001", 
    "HILLIARD_S27_L002"
]

rule all:
    input:
        expand("fastqc_out/{sample}_interleaved_fastqc.html", sample = SAMPLES), 
        expand("data/trimmed/{sample}_interleaved.trimmed.fastq.gz", sample = SAMPLES)

rule fastqc_raw:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        "fastqc_out/{sample}_fastqc.html"
    params:
        outdir = "fastqc_out"
    shell:
        "fastqc {input} -o {params.outdir}"

rule fastp_trim:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        "data/trimmed/{sample}.trimmed.fastq.gz"
    shell:
        "fastp --interleaved_in -i {input} -o {output}"

rule fastqc_trimmed:
    input:
        "data/trimmed/{sample}.fastq.gz"
    output:
        "fastqc_trimmed/{sample}_fastqc.html"
    params:
        outdir = "fastqc_trimmed"
    shell:
        "fastqc {input} -o {params.outdir}"