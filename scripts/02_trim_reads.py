#!/bin/env python3

# Given a directory containing paired fastq files, this script trim reads using trimmomatic.
# REQUIRED SOFTWARES: trimmomatic, fastqc
#
# Note that the structure of the input directory should be as follow:
# AccNo/
# ├── AccNo_1.fastq.gz
# └── AccNo_2.fastq.gz
# This format can be easily obtained using the script https://github.com/filonico/rphi_DE/blob/main/scripts/01_download_reads.py
#
#
# The script creates an output directory structured as follow:
# ./
# └── your_output_dir/
#     ├── 01_fastqc/
#     |   └── {results of fastqc analysis}
#     ├── SRRXXXXXX1_trimmed/
#     |   ├── SRRXXXXXX1_1_paired.fastq.gz
#     |   ├── SRRXXXXXX1_1_unpaired.fastq.gz
#     |   ├── SRRXXXXXX1_2_paired.fastq.gz
#     |   └── SRRXXXXXX1_2_unpaired.fastq.gz
#     ├── SRRXXXXXX2_trimmed/
#     |   ├── SRRXXXXXX2_1_paired.fastq.gz
#     |   ├── SRRXXXXXX2_1_unpaired.fastq.gz
#     |   ├── SRRXXXXXX2_2_paired.fastq.gz
#     |   └── SRRXXXXXX2_2_unpaired.fastq.gz
#     ...
#     └── SRRXXXXXXN_trimmed/
#         ├── SRRXXXXXXN_1_paired.fastq.gz
#         ├── SRRXXXXXXN_1_unpaired.fastq.gz
#         ├── SRRXXXXXXN_2_paired.fastq.gz
#         └── SRRXXXXXXN_2_unpaired.fastq.gz
#
#
# Written by:   Filippo Nicolini
# Last updated: 09/10/2023
#
#------------------------------------------------------------------


import subprocess, argparse, sys, os


############################################
#     Defining arguments of the script     #
############################################

# Initialise the parser class
parser = argparse.ArgumentParser(description = "Trim reads using Trimmomatic and perform quality check.")

# Define some options/arguments/parameters
parser.add_argument("-d", "--input_dir", required = True, help = "Directory containing paired fastq files to trim. Note that the structure of input directory should be as follow: input_dir/{input_dir_1.fastq.gz, input_dir_2.fastq.gz}")
parser.add_argument("-adapt", "--illumina_adapters", required = True, help = "File containing Illumina adapters.")
parser.add_argument("-o", "--output_dir", help = "Name of the output directory.", default = "02_trimmed_reads")

# This line checks if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

# Collect the inputted arguments into a dictionary
args = parser.parse_args()


##########################################################################
#     Defining functions to download reads and perform quality check     #
##########################################################################

# Function to trim reads, given a directory containing paired fastq files
def trim_reads(input_directory, acc, trim_output_dir):
    subprocess.run(f"mkdir {trim_output_dir}", shell = True)

    input_acc = input_directory + "/" + acc
    trim_output_acc = trim_output_dir + "/" + acc

    try:
        trimmomatic_process = subprocess.run(f"trimmomatic PE -threads 15 -phred33 "
                                             f"{input_acc}*1.fastq.gz {input_acc}*2.fastq.gz "
                                             f"{trim_output_acc}_1_paired.fastq.gz {trim_output_acc}_1_unpaired.fastq.gz "
                                             f"{trim_output_acc}_2_paired.fastq.gz {trim_output_acc}_2_unpaired.fastq.gz "
                                             f"ILLUMINACLIP:{args.illumina_adapters}:2:30:10 "
                                             "LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:65 ",
                                             shell = True,
                                             capture_output = True,
                                             text = True)
        
    except subprocess.CalledProcessError as err:
        print("An error occured:", err.stderr)

# Function to perform the quality check, given the fastq file
def quality_check(fastq_file):
    try:
        fastqc_process = subprocess.run(f"fastqc {fastq_file} -o {args.output_dir}/01_fastqc -f fastq",
                                        shell = True,
                                        capture_output = True,
                                        text = True)

    except subprocess.CalledProcessError as err:
        print("An error occured:", err.stderr)

###############################################
#     Trim read and perform quality check     #
###############################################

# Define a variable to store accession number of your run
if not os.path.basename(args.input_dir):
    ACC = os.path.basename(args.input_dir[:-1])
else:
    ACC = os.path.basename(args.input_dir)

# Define a variable to store the name of the output directory of trimmed reads
TRIM_OUTPUT_DIR = args.output_dir + "/" + ACC + "_trimmed"


# Create output direcotory
if not os.path.isdir(args.output_dir):
    print()
    print(f"Creating output directory in {args.output_dir}/")
    subprocess.run(f"mkdir -p {args.output_dir}/01_fastqc", shell = True)


# Trim reads
print()
print(f"-- {ACC} --")
print("  Trimming reads...")
trim_reads(args.input_dir, ACC, TRIM_OUTPUT_DIR)


# Quality check
FASTQ_file = args.output_dir + "/" + ACC + "_trimmed/*_paired.fastq.gz"
print("  Checking read quality...")
quality_check(FASTQ_file)


print("Done")
print()
