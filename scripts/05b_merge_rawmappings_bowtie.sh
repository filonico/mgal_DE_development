#!/bin/bash

if [ -e 03b_mapped_reads_BOWTIE/conditions.tsv ]; then
    
    rm 03b_mapped_reads_BOWTIE/conditions.tsv

fi


for i in 03b_mapped_reads_BOWTIE/SRR*tsv; do

    TMP="03b_mapped_reads_BOWTIE/TMP" &&
    ACC="$(basename $i | awk -F "_" '{print $1}')" &&
    # COND="$(grep $ACC 00_input_files/tech_replicates.tsv | awk '{print $NF}')" &&

    # create a table with the list of experiment conditions in the reading order (we will use this for DE analysis with NOIseq)
    # NB: this file is the same as 00_input/tech_replicates.tsv, but experiments are in the order that is expected for the DE analysis (i.e., the same of the columns in the rawcount table) 
    # echo -e $ACC'\t'$COND >> 03b_mapped_reads_BOWTIE/conditions.tsv &&

    if [ ! -e $TMP ]; then

        # select just the column of mapped reads
        awk '{print $1","$3}' $i | sed -E "1i gene,$ACC" | head -n -1 > $TMP &&
        continue

    else

        if [ -e 03b_mapped_reads_BOWTIE/ALL.rawmapping.stats.csv ]; then

            rm 03b_mapped_reads_BOWTIE/ALL.rawmapping.stats.csv

        fi

        awk '{print $1","$3}' $i | sed -E "1i gene,$ACC" | head -n -1 | join -j 1 -t "," $TMP - > 03b_mapped_reads_BOWTIE/ALL.rawmapping.stats.csv &&
        cp 03b_mapped_reads_BOWTIE/ALL.rawmapping.stats.csv $TMP

    fi

done

rm 03b_mapped_reads_BOWTIE/TMP
