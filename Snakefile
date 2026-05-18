configfile: "config.yaml"

# function to list the outputs of GO term enrichment per maSigPro module
def get_go_outputs(wildcards):
    checkpoint_output = checkpoints.run_maSigPro.get().output.modules_dir
    return expand("05_masigpro_analysis/01_genes_per_module/module_{mSigPro_modules}_genes_GOterms.tsv",
                  mSigPro_modules = glob_wildcards(checkpoint_output + "/module_{mSigPro_modules}_genes.ls").mSigPro_modules)

rule all:
    input:
        get_go_outputs


##########################
#     DOWNLOAD READS     #
##########################

# download reads through the NCBI utilities and run QC
rule download_reads:
    output:
        fastq_1 = "01_raw_reads/{library}/{library}_1.fastq.gz",
        fastq_2 = "01_raw_reads/{library}/{library}_2.fastq.gz",
        fastqc_1 = "01_raw_reads/01_fastqc/{library}_1_paired_fastqc.html",
        fastqc_2 = "01_raw_reads/01_fastqc/{library}_2_paired_fastqc.html"
    conda:
        "conda_envs/ncbi_env.yaml"
    log:
        stdout = "01_raw_reads/{library}/{library}_download.stdout",
        stderr = "01_raw_reads/{library}/{library}_download.stderr"
    shell:
        """
        python3 scripts/01_download_reads.py \
            -i {wildcards.library} \
            -o 01_raw_reads
        """

# aggregate QC reports
rule multiqc_raw_reads:
    input:
        expand("01_raw_reads/01_fastqc/{library}_{read}_paired_fastqc.html",
            library = config["libraries"],
            read = config["read"])
    output:
        "01_raw_reads/01_fastqc/multiqc_report.html"
    log:
        stdout = "01_raw_reads/01_fastqc/multiqc.stdout",
        stderr = "01_raw_reads/01_fastqc/multiqc.stderr"
    conda:
        "conda_envs/multiqc_env.yaml"
    shell:
        "multiqc -o 01_raw_reads/01_fastqc/ 01_raw_reads/01_fastqc/" 


######################
#     TRIM READS     #
######################

# trim reads from adapters and low-quality bases, and run QC
rule trim_reads:
    input:
        fastq_1 = "01_raw_reads/{library}/{library}_1.fastq.gz",
        fastq_2 = "01_raw_reads/{library}/{library}_2.fastq.gz"
    output:
        fastq_1 = "02_trimmed_reads/{library}_trimmed/{library}_1_paired.fastq.gz",
        fastq_2 = "02_trimmed_reads/{library}_trimmed/{library}_2_paired.fastq.gz",
        fastqc_1 = "02_trimmed_reads/01_fastqc/{library}_1_paired_fastqc.html",
        fastqc_2 = "02_trimmed_reads/01_fastqc/{library}_2_paired_fastqc.html"
    conda:
        "conda_envs/readtrimming_env.yaml"
    log:
        stdout = "02_trimmed_reads/{library}_trimmed/{library}_trim.stdout",
        stderr = "02_trimmed_reads/{library}_trimmed/{library}_trim.stderr"
    resources:
        runtime = 300,
        mem_mb = 30000,
        cpus_per_task = 20
    shell:
        """
        python3 scripts/02_trim_reads.py \
            -d 01_raw_reads/{wildcards.library} \
            -adapt 00_input_files/contaminants2trimm.fa
        """

# aggregate QC reports
rule multiqc_trimmed_reads:
    input:
        expand("02_trimmed_reads/01_fastqc/{library}_{read}_paired_fastqc.html",
            library = config["libraries"],
            read = config["read"])
    output:
        "02_trimmed_reads/01_fastqc/multiqc_report.html"
    log:
        stdout = "02_trimmed_reads/01_fastqc/multiqc.stdout",
        stderr = "02_trimmed_reads/01_fastqc/multiqc.stderr"
    conda:
        "conda_envs/multiqc_env.yaml"
    shell:
        "multiqc -o 02_trimmed_reads/01_fastqc/ 02_trimmed_reads/01_fastqc/" 


###############################
#     MAP READS WITH STAR     #
###############################

# extract the genome files
rule tarunzip_genome:
    input:
        targzipped_genome = "00_input/mgal_genome.tar.gz"
    output:
        genome_fna = "00_input/mgal_genome/GCA.900618805.1_mgal_genomic.fna",
        genome_gff = "00_input/mgal_genome/GCA.900618805.1_mgal_genomic_noIso.gff",
        genome_cds = "00_input/mgal_genome/GCA.900618805.1_mgal_cds_noIso.fna"
    shell:
        "tar -xvzf {input.targzipped_genome} -C 00_input"

# map reads with STAR against the annotated genome
rule map_reads_star:
    input:
        fastq_1 = "02_trimmed_reads/{library}_trimmed/{library}_1_paired.fastq.gz",
        fastq_2 = "02_trimmed_reads/{library}_trimmed/{library}_2_paired.fastq.gz",
        genome_fna = "00_input/mgal_genome/GCA.900618805.1_mgal_genomic.fna",
        genome_gff = "00_input/mgal_genome/GCA.900618805.1_mgal_genomic_noIso.gff"
    output:
        mapped_bam = "03a_mapped_reads_STAR/{library}_trimmed_Aligned.sortedByCoord.out.bam"
    log:
        stdout = "03a_mapped_reads_STAR/{library}_star.stdout",
        stderr = "03a_mapped_reads_STAR/{library}_star.stderr"
    params:
        fastq_directory = "02_trimmed_reads/{library}_trimmed"
    conda:
        "conda_envs/mapreads_star_env.yml"
    resources:
        runtime = 1380,
        mem_mb = 200000,
        cpus_per_task = 20
    shell:
        """
        python scripts/04a_map_reads_star.py \
            -d {params.fastq_directory} \
            -i 03a_mapped_reads_STAR/01_genome_index \
            -r {input.genome_fna} \
            -a {input.genome_gff} \
            -o 03a_mapped_reads_STAR
        """

# merge STAR raw counts into one file
rule merge_star_rawcounts_stringtie:
    input:
        mapped_bam = expand("03a_mapped_reads_STAR/{library}_trimmed_Aligned.sortedByCoord.out.bam", library = config["libraries"]),
        genome_gff = "00_input/mgal_genome/GCA.900618805.1_mgal_genomic_noIso.gff"
    output:
        mapped_reads = "03a_mapped_reads_STAR/transcript_count_matrix.csv"
    log:
        stdout = "03a_mapped_reads_STAR/stringtie.stdout",
        stderr = "03a_mapped_reads_STAR/stringtie.stderr"
    conda:
        "conda_envs/countreads_stringtie_env.yml"
    shell:
        "bash scripts/05a_getcounts_stringtie.sh"


#################################
#     MAP READS WITH BOWTIE     #
#################################

# map reads with bowtie against the predicted transcriptome
rule map_reads_bowtie:
    input:
        fastq_directory = "02_trimmed_reads/{library}_trimmed/{library}_1_paired.fastq.gz",
        fastq_2 = "02_trimmed_reads/{library}_trimmed/{library}_2_paired.fastq.gz",
        genome_cds = "00_input/mgal_genome/GCA.900618805.1_mgal_cds_noIso.fna"
    output:
        raw_counts = "03b_mapped_reads_BOWTIE/{library}.rawmapping.stats.tsv"
    log:
        stdout = "03b_mapped_reads_BOWTIE/{library}_star.stdout",
        stderr = "03b_mapped_reads_BOWTIE/{library}_star.stderr"
    params:
        fastq_directory = "02_trimmed_reads/{library}_trimmed"
    conda:
        "conda_envs/mapreads_bowtie_env.yml"
    resources:
        runtime = 1380,
        mem_mb = 200000,
        cpus_per_task = 20
    shell:
        """
        python scripts/04b_map_reads_bowtie.py \
            -d {params.fastq_directory} \
            -ref {input.genome_fna} \
            -o 03b_mapped_reads_BOWTIE
        """

# merge bowtie raw counts into one file
rule merge_bowtie_rawcounts:
    input:
        raw_counts = "03b_mapped_reads_BOWTIE/{library}.rawmapping.stats.tsv"
    output:
        merged_raw_counts = "03b_mapped_reads_BOWTIE/ALL.rawmapping.stats.csv"
    log:
        stdout = "03b_mapped_reads_BOWTIE/merge_bowtie.stdout",
        stderr = "03b_mapped_reads_BOWTIE/merge_bowtie.stderr"
    shell:
        "bash scripts/05b_merge_rawmappings_bowtie.sh"


####################
#     PLOT PCA     #
####################

# run PCA and save plots on raw counts
rule plot_pca_star:
    input:
        mapped_reads = "03a_mapped_reads_STAR/transcript_count_matrix.csv"
    output:
        pca_pdf = "04_PCA_readcounts/PCA_generawcounts.pdf",
        pca_png = "04_PCA_readcounts/PCA_generawcounts.png",
        normalised_counts = "04_PCA_readcounts/gene_counts_vstTransformed_norm.tsv"
    log:
        stdout = "04_PCA_readcounts/plot_pca.stdout",
        stderr = "04_PCA_readcounts/plot_pca.stderr"
    conda:
        "conda_envs/R_env.yml"
    shell:
        "Rscript scripts/06_plotPCA_readcounts.R"


###########################
#     RUN DE ANALYSIS     #
###########################

# run the time-series DE analysis
# this is a checkpoint rule, as we then need to check how many modules maSigPro identified
checkpoint run_maSigPro:
    input:
        mapped_reads = "03a_mapped_reads_STAR/transcript_count_matrix.csv",
        normalised_counts = "04_PCA_readcounts/gene_counts_vstTransformed_norm.tsv"
    output:
        maSigPro_Rdata = "05_masigpro_analysis/DE_timeseries.Rdata",
        maSigPro_modules_png = "05_masigpro_analysis/vstnorm_alltimepoints_maSigPro_clusters.png",
        maSigPro_modules_pdf = "05_masigpro_analysis/vstnorm_alltimepoints_maSigPro_clusters.pdf",
        modules_dir = directory("05_masigpro_analysis/01_genes_per_module")
    log:
        stdout = "05_masigpro_analysis/masigpro.stdout",
        stderr = "05_masigpro_analysis/masigpro.stderr"
    conda:
        "conda_envs/R_env.yml"
    shell:
        "Rscript scripts/07_DE_timeseries.R > {log.stdout} 2> {log.stderr}"

# list the modules identified by maSigPro
def get_modules(wildcards):
    checkpoint_output = checkpoints.run_maSigPro.get().output.modules_dir
    return expand("05_masigpro_analysis/01_genes_per_module/module_{mSigPro_modules}_genes.ls",
                  mSigPro_modules = glob_wildcards(checkpoint_output + "/module_{mSigPro_modules}_genes.ls").mSigPro_modules)


##############################
#     GO TERM ENRICHMENT     #
##############################

# perform GO enrichment on maSigPro modules
rule perform_go_enrich:
    input:
        mSigPro_modules = get_modules
    output:
        go_modules_prefix = "05_masigpro_analysis/01_genes_per_module/module_{mSigPro_modules}_genes_GOterms_"
    shell:
        """
        Rscript scripts/08_perform_GOenrich.R \
            05_masigpro_analysis/01_genes_per_module/gene_universe_GOterms.tsv \
            {input.mSigPro_modules} \
            {output.go_modules_prefix}
        """
