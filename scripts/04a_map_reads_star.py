#!/bin/env python3

# Given a directory containing paired and trimmed fastq files, this script maps reads against a reference genome, also modeling intron splits.
# REQUIRED SOFTWARES: STAR
#
# Note that the structure of the input directory should be as follow:
# AccNo/
# ├── AccNo_1_paired.fastq.gz
# └── AccNo_2_paired.fastq.gz
# This format can be easily obtained using the script scripts/02_trim_reads.py
#
# The script creates an output directory structured as follow:
# ./
# └── your_output_dir/
#     ├── AccNo_trimmed_Aligned.out.bam
#     ├── AccNo_trimmed_Aligned.sortedByCoord.out.bam
#     ├── AccNo_trimmed_Log.final.out
#     ├── AccNo_trimmed_Log.out
#     ├── AccNo_trimmed_Log.progress.out
#     ├── AccNo_trimmed_SJ.out.tab
#     └── AccNo_trimmed__STARgenome
#
#
# Written by:   Filippo Nicolini
# Last updated: 27/11/2023
#
#------------------------------------------------------------------


import subprocess, argparse, sys, os


############################################
#     Defining arguments of the script     #
############################################

# Initialise the parser class
parser = argparse.ArgumentParser(description = "Map reads against a reference genome using STAR, than index and sort the resulting bam file.")

# Define some options/arguments/parameters
parser.add_argument("-d", "--input_dir",
                    help = "Directory containing trimmed paired fastq files to map. Note that the structure of input directory should be as follow: input_dir/{input_dir_1.fastq.gz, input_dir_2.fastq.gz}",
                    required = True)

parser.add_argument("-i", "--genome_index_directory",
                    help = "Directory where to store genome indexes.",
                    required = True)

parser.add_argument("-r", "--reference_genome",
                    help = "Reference genome fasta file used to map reads.",
                    required = True)

parser.add_argument("-a", "--reference_annotation",
                    help = "Reference gtf annotation file used to map reads.",
                    required = True)

parser.add_argument("-o", "--output_dir", 
                    help = "Name of the output directory.",
                    required = True)

# This line checks if the user gave no arguments, and if so then print the help
parser.parse_args(args = None if sys.argv[1:] else ["--help"])

# Collect the inputted arguments into a dictionary
args = parser.parse_args()


###########################################
#     Defining functions to map reads     #
###########################################

# Function to index the reference genome
def index_genome(genomedir, genomefasta, genomegtf):
    try:
        starindex_process = subprocess.run("STAR --runMode genomeGenerate "
                                             "--runThreadN 15 "
                                             f"--genomeDir {genomedir} "
                                             f"--genomeFastaFiles {genomefasta} "
                                             "--sjdbGTFtagExonParentTranscript Parent "
                                             "--sjdbGTFtagExonParentGene ID "
                                             f"--sjdbGTFfile {genomegtf}",
                                             shell = True,
                                             capture_output = True)
        
        starindex_process.check_returncode()

    except subprocess.CalledProcessError as err:
        print("An error occured:", err.stderr)

def map_reads(genomedir, readfiles, genomegtf, output_prefix):

    try:
        star_process = subprocess.run("STAR --runMode alignReads "
                                "--runThreadN 15 "
                                f"--genomeDir {genomedir} "
                                f"--readFilesIn {readfiles} "
#                                "--readFilesCommand zcat "
                                f"--sjdbGTFfile {genomegtf} "
				"--sjdbGTFtagExonParentTranscript Parent "
				"--sjdbGTFtagExonParentGene ID "
                                f"--outFileNamePrefix {output_prefix} "
                                "--outSAMtype BAM Unsorted SortedByCoordinate",
                                shell = True,
                                capture_output = True,
                                text = True)
        
        star_process.check_returncode()

    except subprocess.CalledProcessError as err:
        print("An error occured:", err)



#---------------------------------------------------------------------------------

###################################################
#     Create directories and index the genome     #
###################################################

print()

# Create output directory and index directory
if not os.path.isdir(args.output_dir):
    print(f"Creating output directory in {args.output_dir}/")
    subprocess.run(f"mkdir {args.output_dir}", shell = True)

if not os.path.isdir(args.genome_index_directory):
    print(f"Creating directory for indexed genome in {args.genome_index_directory}/")
    subprocess.run(f"mkdir {args.genome_index_directory}", shell = True)


# Index the reference genome
genomeparameter = args.genome_index_directory + "/genomeParameters.txt"

if not os.path.isfile(genomeparameter):
    print()

    file_to_use = []
    
    # Check if fasta or gtf files of the reference genome are gzipped: in this case, unzip them 
    for file in args.reference_genome,args.reference_annotation:
        
        if file.endswith(".gz"):
            print(f"Gunzipping {file}...")
            subprocess.run(f"gunzip {file}", shell = True)
            file_to_use.append(file[:-3])
        else:
            print(f"{file} already unzipped.")
            file_to_use.append(file)

    print()
    print("Indexing reference genome...")
    index_genome(args.genome_index_directory, file_to_use[0], file_to_use[1])

else:
    print()
    print("Reference genome already indexed.")

#####################
#     Map reads     #
#####################

# Define a variable to store the accession number of your run
if not os.path.basename(args.input_dir):
    ACC = os.path.basename(args.input_dir[:-1])
else:
    ACC = os.path.basename(args.input_dir)

print()
print(f"-- {ACC} --")

if args.input_dir.endswith("/"):
    READFILES = args.input_dir + "*_paired.fastq"
else:
    READFILES = args.input_dir + "/*_paired.fastq"

OUTPUT_PREFIX = args.output_dir + "/" + ACC + "_"

print("  Mapping reads...")
map_reads(args.genome_index_directory, READFILES, args.reference_annotation, OUTPUT_PREFIX)

print("Done")
print()
