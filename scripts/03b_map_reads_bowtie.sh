#!/bin/bash

for i in 02_trimmed_reads/SRR*; do
	
	python scripts/04a_map_reads_bowtie.py \
		-d $i \
		-ref 00_input/mgal_genome/GCA.900618805.1_mgal_cds_noIso.fna \
		-o 03b_mapped_reads_BOWTIE

done
