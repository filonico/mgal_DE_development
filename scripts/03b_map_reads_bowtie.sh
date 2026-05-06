#!/bin/bash

source ~/miniforge3/bin/activate mapreads_bowtie_env 

for i in 02_trimmed_reads/SRR*; do
	
	python scripts/04b_map_reads_bowtie.py \
		-d $i \
		-ref 00_input/mgal_genome/GCA.900618805.1_mgal_cds_noIso.fna \
		-o 03b_mapped_reads_BOWTIE

done
