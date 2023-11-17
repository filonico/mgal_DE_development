#!/bin/env python3

# Given a directory containing paired and trimmed fastq files, this script map reads against a reference transcriptome.
# REQUIRED SOFTWARES: bowtie2
#
# Note that the structure of the input directory should be as follow:
# AccNo/
# ├── AccNo_1_paired.fastq.gz
# └── AccNo_2_paired.fastq.gz
# This format can be easily obtained using the script https://github.com/filonico/misc_scripts/blob/main/trim_reads.py
#
# The script creates an output directory structured as follow:
# ./
# └── your_output_dir/
#     ├── AccNo.mapped.bam
#     ├── AccNo.mapped.log
#     ├── AccNo.mapped.sorted.filtered.bam
#     ├── AccNo.mapped.sorted.filtered.bam.bai
#     └── AccNo.rawmapping.stats.tsv
#
#
# Written by:   Filippo Nicolini
# Last updated: 17/11/2023
#
#------------------------------------------------------------------


import subprocess, argparse, sys, os


############################################
#     Defining arguments of the script     #
############################################

# Initialise the parser class
parser = argparse.ArgumentParser(description = "Trim reads using Trimmomatic and perform quality check.")

# Define some options/arguments/parameters
parser.add_argument("-d", "--input_dir", required = True, help = "Directory containing trimmed paired fastq files to map. Note that the structure of input directory should be as follow: input_dir/{input_dir_1.fastq.gz, input_dir_2.fastq.gz}")
parser.add_argument("-ref", "--reference_transcriptome", required = True, help = "Reference transcriptome used to map reads.")
parser.add_argument("-o", "--output_dir", help = "Name of the output directory.", default = "04_mappings")

# This line checks if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

# Collect the inputted arguments into a dictionary
args = parser.parse_args()


###########################################
#     Defining functions to map reads     #
###########################################

# Function to index the reference transcriptome
def index_transcriptome(transcriptome):
    try:
        bowtiebuild_process = subprocess.run(f"bowtie2-build {transcriptome} {transcriptome}",
                                             shell = True)
    except:
        print("An error occured:", err.stderr)

def map_reads(input_directory, acc, indexed_transcriptome, output_directory):
    try:
        bowtie_process = subprocess.run(f"bowtie2 -x {indexed_transcriptome} "
                                        f"-1 {input_directory}/{acc}_1_paired.fastq.gz -2 {input_directory}/{acc}_2_paired.fastq.gz "
                                        f"-S {output_directory}/{acc}.mapped.sam --no-discordant -p 30 "
                                        f"2> {output_directory}/{acc}.mapped.log",
                                        shell = True,
                                        capture_output = True,
                                        text = True)
        
    except subprocess.CalledProcessError as err:
        print("An error occured:", err.stderr)

def from_sam_to_rawcounts(sam_file, acc, indexed_transcriptome, output_directory):
    output_acc = output_directory + "/" + acc

    try:
        # Generate bam file from sam file
        generate_bam_process = subprocess.run(f"samtools view -b {sam_file} > {output_acc}.mapped.bam",
                                              shell = True,
                                              capture_output = True,
                                              text = True)
        
        # Remove sam file
        subprocess.run(f"rm {output_acc}.mapped.sam", shell = True)

        # Sort and filter bam file
        sort_n_filter_process = subprocess.run(f"samtools sort -@ 30 {output_acc}.mapped.bam | "
                                               f"samtools view -b -t {indexed_transcriptome} -F 4 -h -@ 30 "
                                               f"> {output_acc}.mapped.sorted.filtered.bam",
                                               shell = True,
                                               capture_output = True,
                                               text = True)
        
        # Get raw counts statistics out of bam file
        get_rowcounts_process = subprocess.run(f"samtools index {output_acc}.mapped.sorted.filtered.bam && "
                                               f"samtools idxstats {output_acc}.mapped.sorted.filtered.bam "
                                               f"> {output_acc}.rawmapping.stats.tsv",
                                               shell = True,
                                               capture_output = True,
                                               text = True)
    
    except subprocess.CalledProcessError as err:
        print("An error occured:", err.stderr)



#####################
#     Map reads     #
#####################

# Define a variable to store the accession number of your run
if not os.path.basename(args.input_dir):
    ACC = os.path.basename(args.input_dir[:-1])
else:
    ACC = os.path.basename(args.input_dir)


# Create output directory
if not os.path.isdir(args.output_dir):
    print()
    print(f"Creating output directory in {args.output_dir}/")
    subprocess.run(f"mkdir {args.output_dir}", shell = True)


# Index the reference transcriptome
if not any(file.endswith(".bt2") for file in os.listdir(os.path.dirname(args.reference_transcriptome))):
    print("Indexing reference transcriptome...")
    index_transcriptome(args.reference_transcriptome)


# Map reads
print()
print(f"-- {ACC} --")
print("  Mapping reads...")
map_reads(args.input_dir, ACC, args.reference_transcriptome, args.output_dir)


# Convert to bam and get raw count statistics
SAM_file = args.output_dir + "/" + ACC + ".mapped.sam"

print("  Converting sam to bam and getting raw counts...")
from_sam_to_rawcounts(SAM_file, ACC, args.reference_transcriptome, args.output_dir)

print("Done")
print()
