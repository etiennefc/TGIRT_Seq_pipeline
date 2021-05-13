import os
from pathlib import Path

configfile: "config.json"

original_id = list(config['dataset'].values())
simple_id = list(config['dataset'].keys())


rule all:
    input:
        qc_before_trim = expand("logs/fastqc/before_trim/{id}.log",
            id=simple_id),
        qc_after_trim = expand("logs/fastqc/after_trim/{id}.log", id=simple_id),
        coco_cc = expand('results/coco/{id}.tsv', id=simple_id)


rule all_downloads:
    input:
        samples = expand('data/references/fastq/{id}_{pair}.fastq.gz',
            id=original_id, pair=[1, 2]),
        reference_genome = config['path']['reference_genome'],
        standard_gtf = config['path']['standard_annotation'],
        reference_gtf = config['path']['reference_annotation'],
        reference_intron_gtf = config['path']['reference_intron_annotation'],
        coco_git = 'git_repos/coco'




rule rename_samples:
    """Rename samples with a nicer understandable name"""
    input:
        fastq = expand(
            "data/references/fastq/{id}_{pair}.fastq.gz",
            id=original_id, pair=[1, 2])
    output:
        renamed_fastq = expand(
            "data/references/fastq/{id}_R{pair}.fastq.gz",
            id=simple_id, pair=[1, 2])
    run:
        for new_name, old_name in config['dataset'].items():
            for num in [1, 2]:
                old = "data/references/fastq/{}_{}.fastq.gz".format(old_name, num)
                new = "data/references/fastq/{}_R{}.fastq.gz".format(new_name, num)
                print(old, new)
                os.rename(old, new)


rule trimming:
    """Trims the input FASTQ files using Trimmomatic"""
    input:
        fastq1 = "data/references/fastq/{id}_R1.fastq.gz",
        fastq2 = "data/references/fastq/{id}_R2.fastq.gz"
    output:
        fastq1 = "data/Trimmomatic/trimmed_reads/{id}_R1.fastq.gz",
        fastq2 = "data/Trimmomatic/trimmed_reads/{id}_R2.fastq.gz",
        unpaired_fastq1 = "data/Trimmomatic/trimmed_reads/{id}_R1.unpaired.fastq.gz",
        unpaired_fastq2 = "data/Trimmomatic/trimmed_reads/{id}_R2.unpaired.fastq.gz"
    threads:
        32
    params:
        options = [
            "ILLUMINACLIP:data/Trimmomatic/Adapters-PE_NextSeq.fa:2:12:10:8:true",
            "TRAILING:30", "LEADING:30", "MINLEN:20"
        ]
    log:
        "logs/trimmomatic/{id}.log"
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "-phred33 "
        "{input.fastq1} {input.fastq2} "
        "{output.fastq1} {output.unpaired_fastq1} "
        "{output.fastq2} {output.unpaired_fastq2} "
        "{params.options} "
        "&> {log}"


rule qc_before_trim:
    "Assess fastq quality before trimming reads"
    input:
        fastq1 = "data/references/fastq/{id}_R1.fastq.gz",
        fastq2 = "data/references/fastq/{id}_R2.fastq.gz"
    output:
        qc_report1 = "data/FastQC/Before_trim/{id}_R1_fastqc.html",
        qc_report2 = "data/FastQC/Before_trim/{id}_R2_fastqc.html"
    threads:
        32
    params:
        out_dir = "data/FastQC/Before_trim"
    log:
        "logs/fastqc/before_trim/{id}.log"
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc "
        "-f fastq "
        "-t {threads} "
        "-o {params.out_dir} "
        "{input.fastq1} "
        "{input.fastq2} "
        "&> {log}"


rule qc_after_trim:
    "Assess fastq quality after trimming reads"
    input:
        fastq1 = rules.trimming.output.fastq1,
        fastq2 = rules.trimming.output.fastq2
    output:
        qc_report1 = "data/FastQC/After_trim/{id}_R1_fastqc.html",
        qc_report2 = "data/FastQC/After_trim/{id}_R2_fastqc.html"
    threads:
        32
    params:
        out_dir = "data/FastQC/After_trim"
    log:
        "logs/fastqc/after_trim/{id}.log"
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc "
        "-f fastq "
        "-t {threads} "
        "-o {params.out_dir} "
        "{input.fastq1} "
        "{input.fastq2} "
        "&> {log}"


rule star_index:
    """Generate the genome index needed for STAR alignment"""
    input:
        fasta = config['path']['reference_genome'],
        standard_gtf = config['path']['standard_annotation']
    output:
        chrNameLength = "data/star_index/chrNameLength.txt"
    threads:
        32
    params:
        index_dir = "data/star_index/"
    log:
        "logs/star/star_index.log"
    conda:
        "envs/star.yaml"
    shell:
        "STAR --runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.standard_gtf} "
        "--sjdbOverhang 74 "
        "&> {log}"


rule star_align:
    """Align reads to reference genome using STAR"""
    input:
        fastq1 = rules.trimming.output.fastq1,
        fastq2 = rules.trimming.output.fastq2,
        idx = rules.star_index.output
    output:
        bam = "results/star/{id}/Aligned.sortedByCoord.out.bam"
    threads:
        32
    params:
        outdir = "results/star/{id}/",
        index_dir = "data/star_index/"
    log:
        "logs/star/star_align_{id}.log"
    conda:
        "envs/star.yaml"
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {params.index_dir} "
        "--readFilesIn {input.fastq1} {input.fastq2} "
        "--runThreadN {threads} "
        "--readFilesCommand zcat "
        "--outReadsUnmapped Fastx "
        "--outFilterType BySJout "
        "--outStd Log "
        "--outSAMunmapped None "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.outdir} "
        "--outFilterScoreMinOverLread 0.3 "
        "--outFilterMatchNminOverLread 0.3 "
        "--outFilterMultimapNmax 100 "
        "--winAnchorMultimapNmax 100 "
        "--alignEndsProtrude 5 ConcordantPair"
        "&> {log}"


rule coco_cc:
    """Quantify the number of counts, counts per million (CPM) and transcript
        per million (TPM) for each gene using CoCo correct_count (cc). If you
        want to use your own annotation file, you must implement in this
        workflow the correct_annotation function of CoCo in order to use CoCo's
        correct_count function. Otherwise, by default, you will use the given 
        annotation file used in this analysis"""
    input:
        gtf = config['path']['reference_annotation'],
        intron_gtf = config['path']['reference_intron_annotation'],
        bam = rules.star_align.output.bam
    output:
        counts = Path("results/coco/", "{id}.tsv")
    threads:
        32
    params:
        coco_path = "git_repos/coco/bin"
    log:
        "logs/coco/coco_{id}.log"
    conda:
        "envs/coco.yaml"
    shell:
        "python {params.coco_path}/coco.py cc "
        "--countType both "
        "--thread {threads} "
        "--strand 1 "
        "--paired "
        "{input.gtf} "
        "{input.bam} "
        "{output.counts} "
        "&> {log}"


#Include download_all rules to download necessary datasets and references
include: "rules/download_all.smk"
