#/project/gbru_wheat2/fhb/conda/exomecluster_env

import os

#wildcards = glob_wildcards("data/raw/{sample}_interleaved.fastq.gz")
#SAMPLES = wildcards.sample

SAMPLES = [
    "AGS2000-RALEIGH-SRWW_S1_L001", 
    "AGS2000-RALEIGH-SRWW_S9_L002", 
    "AGS2000_S1_L001", 
    "AGS2000_S9_L002", 
    "HILLIARD_S27_L001", 
    "HILLIARD_S27_L002"
]

#SAMPLES = "AGS2000-RALEIGH-SRWW_S1_L001"

REF_ACCESSION = config['REF_ACCESSION']
REF_FASTA = config['REF_FASTA']
REF_DICT = os.path.splitext(REF_FASTA)[0] + '.dict'
PICARD_JAR = config['PICARD_JAR']

rule all:
    input:
        expand("fastqc_out/raw/{sample}_interleaved_fastqc.html", sample = SAMPLES), 
        expand("data/trimmed/{sample}.F.trimmed.fastq", sample = SAMPLES), 
        expand("data/trimmed/{sample}.R.trimmed.fastq", sample = SAMPLES), 
        expand("fastqc_out/trimmed/{sample}_interleaved.trimmed_fastqc.html", sample = SAMPLES), 
        expand("{accession}.{int}.ht21", accession = REF_ACCESSION, int = "1"), 
        expand("alignment/{sample}.sam", sample = SAMPLES), 
        expand("samtools_stats/{sample}.txt", sample = SAMPLES), 
        expand("calls/{sample}.g.vcf.gz", sample = SAMPLES), 
        REF_DICT, 
        "db",
        "srww_exome.vcf"

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
        "data/raw/{sample}_interleaved.fastq.gz"
    output:
        out1 = "data/trimmed/{sample}.F.trimmed.fastq",
        out2 = "data/trimmed/{sample}.R.trimmed.fastq" 
    params:
        json = "reports/fastqc_trimmed/{sample}.json", 
        html = "reports/fastqc_trimmed/{sample}.html"
    log:
        "logs/fastp/{sample}.log"
    shell:
        "fastp --interleaved_in -j {params.json} -h {params.html} "
        "-i {input} --out1 {output.out1} --out2 {output.out2} 2> {log}"

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

rule hisat2_align:
    input:
        in1 = "data/trimmed/{sample}.F.trimmed.fastq",
        in2 = "data/trimmed/{sample}.R.trimmed.fastq" 
    output:
        sam = "alignment/{sample}.sam", 
        bam = "alignment/{sample}.bam", 
        sorted_bam = "alignment/{sample}.sorted.bam"
    params:
        accession = REF_ACCESSION
    log:
        "logs/hisat2/{sample}.log"
    threads: 8
    shell:
        "hisat2 -x {params.accession} -p {threads} --rg-id {wildcards.sample} --rg SM:{wildcards.sample} "
        "-1 {input.in1} -2 {input.in2} > {output.sam} 2> {log} ; "
        "samtools view -Sb  -F 256 {output.sam} > {output.bam} 2>> {log} ; " #-F 256 removes non-primary alignments
        "samtools sort {output.bam} -o {output.sorted_bam} 2>> {log} ; "
        "samtools index -c {output.sorted_bam} 2>> {log}"

rule mark_dups:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        bam = "alignment/{sample}.dedup.bam", 
        metrics = "logs/markdups/{sample}_metrics.txt"
    log:
        "logs/markdups/{sample}.log"
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} 2> {log} ; "
        "samtools index -c {output.bam} 2>> {log}"

rule samtools_stats:
    input:
        "alignment/{sample}.dedup.bam"
    output:
        "samtools_stats/{sample}.txt"
    log:
        "logs/samtools_stats/{sample}.log"
    threads: 8
    shell:
        "samtools stats {input} --threads {threads} > {output} 2> {log}"

rule create_sequence_dictionary:
    input:
        ref = REF_FASTA
    output:
        dict = REF_DICT
    params:
        picard = PICARD_JAR
    shell:
        "java -Xms5g -jar {params.picard} CreateSequenceDictionary -R {input} -O {output}"

rule haplotypecaller:
    input:
        ref = REF_FASTA, 
        dict = REF_DICT, 
        bam = "alignment/{sample}.dedup.bam"
    output:
        "calls/{sample}.g.vcf"
    log:
        "logs/haplotypecaller/{sample}.log"
    threads: 4
    shell:
        """
            gatk --java-options "-Xmx32g" HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
            -ERC GVCF \
            2> {log}
        """

rule makedatabase:
    input:
        gvcf = expand("calls/{sample}.g.vcf", sample = SAMPLES), 
        intervals = "files/intervals.list", 
#        gvcf = list(map("--variant {}".format, expand("calls/{sample}.g.vcf", sample = SAMPLES)))
    output:
        db = directory("db"), 
        store = "db/something"
    log:
        "logs/gatk/genomicsdbimport.log"
    threads: 48
    run:
        flagged_gvcf = ["--variant " + f for f in input.gvcf]
        shell("""
            gatk --java-options "-Xmx300g -Xms100g" GenomicsDBImport {flagged_gvcf} \
            --genomicsdb-workspace-path {output.db} -L {input.intervals} 2> {log}
        """)

rule genotype:
    input:
        ref = REF_FASTA, 
        db = "db"
    output:
        "srww_exome.vcf"
    threads: 48
    shell:
        """
            gatk --java-option "-Xmx300g -Xms100g" GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.db} \
            -O {output}
        """
