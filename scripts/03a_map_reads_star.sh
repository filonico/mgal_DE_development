#!/bin/bash

for i in 02_trimmed_reads/SRR*; do
	
	python scripts/04_map_reads_star.py \
		-d $i \
		-i 03a_mapped_reads_STAR/01_genome_index \
		-r 00_input/mgal_genome/GCA.900618805.1_mgal_genomic.fna \
		-a 00_input/mgal_genome/GCA.900618805.1_mgal_genomic.gtf \
		-o 03a_mapped_reads_STAR

done
