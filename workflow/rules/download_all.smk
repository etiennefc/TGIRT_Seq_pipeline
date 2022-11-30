import os

rule download_human_sample_fastq:
    """Download fastq of human samples from GEO. This step
        might take long, depending on your downloading speed. Other samples were 
	manually added on the cluster. To automate the download of ALL samples, 
	simply add their SRR id in the config file (use SRA explorer to find the URL)"""
    output:
        fake_output = "data/references/geo_download.txt"
    params:
        sample_list = "data/references/sra_id_wget.txt"
    shell:
        "mkdir -p data/references/fastq/human/ && "
        "cd data/references/fastq/human/ && "
        "wget -i ../../references/sra_id_wget.txt && "
        "touch ../../../../{params.sample_list}"


rule download_sample_fastq:
    """Download expression datasets of all tissue samples from GEO. This step
        might take long, depending on your downloading speed."""
    output:
        samples_fastq_1_gz = "data/references/fastq/{id}_1.fastq.gz",
        samples_fastq_2_gz = "data/references/fastq/{id}_2.fastq.gz"
    params:
        samples_fastq_1 = "data/references/fastq/{id}_1.fastq",
        samples_fastq_2 = "data/references/fastq/{id}_2.fastq"
    conda:
        "../envs/geo_download.yaml"
    shell:
        "mkdir -p data/references/fastq/ && "
        "fasterq-dump --split-files {wildcards.id} "
        "-O data/references/fastq/ && "
        "gzip -f {params.samples_fastq_1} > {output.samples_fastq_1_gz} && "
        "gzip -f {params.samples_fastq_2} > {output.samples_fastq_2_gz}"


rule download_mammal_genome:
    """Download the reference genome (fasta file) of human and mouse
        from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["mus_musculus",    
              "homo_sapiens"])
    params:
        link = "ftp://ftp.ensembl.org/pub/release-108/fasta/{species}/dna/*dna.primary_assembly.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"


rule download_yeast_genome:
    """Download the reference genome (fasta file) of S. pombe and 
        S. cerevisiae from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["saccharomyces_cerevisiae",
              "schizosaccharomyces_pombe"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-55/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"


rule download_tetrahymena_genome:
    """ Download fasta of Tetrahymena thermophila genome per species from Zenodo/TGD"""
    output:
        genome = 'data/references/genome_fa/tetrahymena_thermophila_genome.fa'
    params:
        link = config['download']['t_thermophila_genome']
    shell:
        "wget -O {output.genome} {params.link}"


rule download_mouse_gtf:
    """Download the annotation of the mouse genome (gtf)
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/mus_musculus.gtf'
    params:
        link = config['download']['mouse_gtf'] 
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"


rule download_yeast_gtf:
    """Download the reference genome (fasta file) of S. pombe and
        S. cerevisiae from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ["saccharomyces_cerevisiae",
              "schizosaccharomyces_pombe"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-55/gtf/{species}/*5.gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"


rule download_tetrahymena_gtf:
    """ Download gtf of Tetrahymena thermophila genome from Zenodo/TGD"""
    output:
        gtf = 'data/references/gtf/tetrahymena_thermophila.gtf'
    params:
        link = config['download']['t_thermophila_gtf']
    shell:
        "wget -O {output.gtf} {params.link}"


rule download_human_gtf:
    """ Download gtf of human genome from Zenodo"""
    output:
        gtf = 'data/references/gtf/homo_sapiens.gtf'
    params:
        link = config['download']['human_gtf']
    shell:
        "wget -O {output.gtf} {params.link}"


rule download_coco_git:
    """Download git repository of CoCo."""
    output:
        git_coco_folder = directory('git_repos/coco')
    params:
        git_coco_link = config['path']['coco_git_link']
    conda:
        '../envs/git.yaml'
    shell:
        'mkdir -p {output.git_coco_folder} '
        '&& git clone {params.git_coco_link} {output.git_coco_folder}'

rule install_pairedBamToBed12:
    output:
        directory("scripts/pairedBamToBed12/")
    params:
        pairedBamToBed12_bin = 'scripts/pairedBamToBed12/bin'
    conda:
        "../envs/coco.yaml"
    shell:
        'mkdir -p scripts && cd scripts && pwd && '
        'git clone https://github.com/Population-Transcriptomics/pairedBamToBed12 && '
        'cd pairedBamToBed12 && '
        'make '
