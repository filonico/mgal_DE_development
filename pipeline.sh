#!/bin/bash

##############################
#     Download raw reads     #
##############################

# REQUIRES: conda_env/ncbi_env.yml
python3 scripts/01_download_reads.py -i 00_input/SRA_toDownload.ls -o 01_rawreads
