rule rename_samples_human:
    """Rename human samples with a nicer understandable name"""
    input:
        fastq = expand(
            "data/references/fastq/human/{id}_{pair}.fastq.gz",
            id=original_id, pair=[1, 2])
    output:
        renamed_fastq = expand(
            "data/references/fastq/human/{id}_R{pair}.fastq.gz",
            id=simple_id, pair=[1, 2])
    run:
        for new_name, old_name in config['dataset_human'].items():
            for num in [1, 2]:
                old = "data/references/fastq/human/{}_{}.fastq.gz".format(old_name, num)
                new = "data/references/fastq/human/{}_R{}.fastq.gz".format(new_name, num)
                print(old, new)
                os.rename(old, new)


rule qc_before_trim:
    "Assess fastq quality before trimming reads"
    input:
        fastq1 = "data/references/fastq/{species}/{id}_R1.fastq.gz",
        fastq2 = "data/references/fastq/{species}/{id}_R2.fastq.gz"
    output:
        qc_report1 = "data/FastQC/Before_trim/{species}/{id}_R1_fastqc.html",
        qc_report2 = "data/FastQC/Before_trim/{species}/{id}_R2_fastqc.html"
    threads:
        32
    params:
        out_dir = "data/FastQC/Before_trim/{species}"
    log:
        "logs/fastqc/before_trim/{species}/{id}.log"
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "-f fastq "
        "-t {threads} "
        "-o {params.out_dir} "
        "{input.fastq1} "
        "{input.fastq2} "
        "&> {log}"


rule trimming:
    """Trims the input FASTQ files using Trimmomatic"""
    input:
        fastq1 = "data/references/fastq/{species}/{id}_R1.fastq.gz",
        fastq2 = "data/references/fastq/{species}/{id}_R2.fastq.gz"
    output:
        fastq1 = "data/Trimmomatic/trimmed_reads/{species}/{id}_R1.fastq.gz",
        fastq2 = "data/Trimmomatic/trimmed_reads/{species}/{id}_R2.fastq.gz",
        unpaired_fastq1 = "data/Trimmomatic/trimmed_reads/{species}/{id}_R1.unpaired.fastq.gz",
        unpaired_fastq2 = "data/Trimmomatic/trimmed_reads/{species}/{id}_R2.unpaired.fastq.gz"
    threads:
        32
    params:
        options = [
            "ILLUMINACLIP:data/Trimmomatic/Adapters-PE_NextSeq.fa:2:12:10:8:true",
            "TRAILING:30", "LEADING:30", "MINLEN:20"
        ]
    log:
        "logs/trimmomatic/{species}/{id}.log"
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "-phred33 "
        "{input.fastq1} {input.fastq2} "
        "{output.fastq1} {output.unpaired_fastq1} "
        "{output.fastq2} {output.unpaired_fastq2} "
        "{params.options} "
        "&> {log}"


rule qc_after_trim:
    "Assess fastq quality after trimming reads"
    input:
        fastq1 = rules.trimming.output.fastq1,
        fastq2 = rules.trimming.output.fastq2
    output:
        qc_report1 = "data/FastQC/After_trim/{species}/{id}_R1_fastqc.html",
        qc_report2 = "data/FastQC/After_trim/{species}/{id}_R2_fastqc.html"
    threads:
        32
    params:
        out_dir = "data/FastQC/After_trim/{species}"
    log:
        "logs/fastqc/after_trim/{species}/{id}.log"
    conda:
        "../envs/fastqc.yaml"
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
        fasta = get_species_genome,
        standard_gtf = get_species_gtf
    output:
        chrNameLength = "data/star_index/{species}/chrNameLength.txt"
    threads:
        32
    params:
        index_dir = "data/star_index/{species}/",
        temp_dir = "./_STARtmp_{species}/"
    log:
        "logs/star_index/{species}.log"
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.standard_gtf} "
        "--sjdbOverhang 74 "
        "--outTmpDir {params.temp_dir} "
        "&> {log}"


rule star_align:
    """Align reads to reference genome using STAR"""
    input:
        fastq1 = rules.trimming.output.fastq1,
        fastq2 = rules.trimming.output.fastq2,
        idx = rules.star_index.output
    output:
        bam = "results/star/{species}/{id}/Aligned.sortedByCoord.out.bam"
    threads:
        32
    params:
        outdir = "results/star/{species}/{id}/",
        index_dir = "data/star_index/{species}/"
    log:
        "logs/star_align/{species}/{id}.log"
    conda:
        "../envs/star.yaml"
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
        "--alignEndsProtrude 5 ConcordantPair "
        "--limitBAMsortRAM 6802950316 "
        "&> {log}"


rule coco_ca:
    """ Generate corrected annotation from the gtf."""
    input:
        gtf = get_species_gtf,
        paired_bam_to_bed12_dependency = rules.install_pairedBamToBed12.output,
        coco_dir = rules.download_coco_git.output.git_coco_folder
    output:
        gtf_corrected = "data/references/gtf/{species}.correct_annotation.gtf"
    params:
        rules.install_pairedBamToBed12.params.pairedBamToBed12_bin
    conda:
        "../envs/coco.yaml"
    shell:
        "export PATH=$PWD/{params}:$PATH && "
        "python3 git_repos/coco/bin/coco.py ca {input.gtf} -o {output.gtf_corrected}"


rule coco_cc:
    """Quantify the number of counts, counts per million (CPM) and transcript
        per million (TPM) for each gene using CoCo correct_count (cc)."""
    input:
        gtf = rules.coco_ca.output.gtf_corrected,
        bam = rules.star_align.output.bam
    output:
        counts = "results/coco_cc/{species}/{id}.tsv"
    threads:
        32
    params:
        coco_path = "git_repos/coco/bin"
    log:
        "logs/coco/{species}/coco_cc_{id}.log"
    conda:
        "../envs/coco.yaml"
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


rule merge_coco_cc_output:
    """ Merge CoCo correct count outputs into one count, cpm or tpm file (all
        samples merged inside one dataframe). This rule takes in input the
        output result directory of CoCo within the TGIRT-Seq pipeline."""
    input:
        counts = get_coco_cc_output_per_species 
    output:
        merged_counts = "results/coco_cc/{species}_merged_counts.tsv",
        merged_cpm = "results/coco_cc/{species}_merged_cpm.tsv",
        merged_tpm = "results/coco_cc/{species}_merged_tpm.tsv"
    conda:
        "../envs/coco.yaml"
    script:
        "../scripts/merge_coco_cc_output.py"


rule coco_cb:
    """ Create a bedgraph from the bam files """
    input:
        bam = rules.star_align.output.bam,
        chrNameLength = rules.star_index.output.chrNameLength
    output:
        unsorted_bedgraph = "results/coco_cb/{species}/{id}_unsorted.bedgraph"
    params:
        pb = rules.install_pairedBamToBed12.params.pairedBamToBed12_bin,
        coco_path = "git_repos/coco/bin"
    conda:
        "../envs/coco.yaml"
    threads:
        28
    shell:
        "export PATH=$PWD/{params.pb}:$PATH && "
        "python {params.coco_path}/coco.py cb "
        "-u " # UCSC compatible (adds a chr and track info)
        "-t {threads} "
        "-c 2500000 " # Chunk size, default value
        "{input.bam} "
        "{output.unsorted_bedgraph} "
        "{input.chrNameLength}"


rule sort_bg:
    """ Sort the new bedgraphs and change the chrM to chrMT to work with bedGraphToBigWig """
    input:
        unsorted_bedgraph = rules.coco_cb.output.unsorted_bedgraph
    output:
        sorted_bedgraph = "results/coco_cb/{species}/{id}.bedgraph"
    shell:
        "sort -k1,1 -k2,2n {input.unsorted_bedgraph} "
        "| sed 's/chrM/chrMT/g' > {output.sorted_bedgraph}"
