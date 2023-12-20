#!/bin/bash

for i in 03_mapped_reads/*sortedByCoord.out.bam; do

	ACC="$(basename $i | awk -F "_" '{print $1}')";

	echo Running stringtie on "$i"...
	
	# run stringtie
	bash scripts/06_run_stringtie.sh $i 00_input/mgal_genome/GCA.900618805.1_mgal_genomic.gtf 03_mapped_reads/"$ACC"_stringtie.gtf 2> 03_mapped_reads/"$ACC"_stringtie.log

	# create input file for prepDE.py script
	echo -e "$ACC"'\t'03_mapped_reads/"$ACC"_stringtie.gtf >> 03_mapped_reads/input_prepDE.tsv

	echo -e Done'\n'

done

echo "Extracting read counts for DE analysis using prepDE.py..."

# run prepDE.py to extract read counts
prepDE.py -i 03_mapped_reads/input_prepDE.tsv -l 99 -g 03_mapped_reads/gene_count_matrix.csv -t 03_mapped_reads/transcript_count_matrix.csv
