#!/bin/bash

for i in 02_trimmed_reads/SRR253874{30..39}*; do
	
	python3 scripts/04a_map_reads_star.py \
		-d "$i" \
		-i 03a_mapped_reads_STAR/01_genome_index \
		-r 00_input/mgal_genome/GCA.900618805.1_mgal_genomic.fna \
		-a 00_input/mgal_genome/GCA.900618805.1_mgal_genomic_noIso.gff \
		-o 03a_mapped_reads_STAR

done
