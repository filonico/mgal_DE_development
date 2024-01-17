#!/bin/bash

# check if directories and files already exist
if [ ! -d 04_mapping_counts/ ]; then
	mkdir 04_mapping_counts
fi

if [ -f 04_mapping_counts/input_prepDE.tsv ]; then
	rm -f 04_mapping_counts/input_prepDE.tsv
fi

# first round of stringtie against reference gff annotation file
echo -e '\n'First round of StringTie \(vs reference\)...

OUTDIR="04_mapping_counts"
REF_GFF="00_input/mgal_genome/GCA.900618805.1_mgal_genomic.gff"

for i in 03_mapped_reads/*sortedByCoord.out.bam; do

	ACC="$(basename $i | awk -F "_" '{print $1}')"

	# run stringtie
	stringtie "$i" -G "$REF_GFF" -o "$OUTDIR"/"$ACC"_stringtie_round1.gtf -p 15 2> "$OUTDIR"/"$ACC"_stringtie_round1.log

	echo -e '\t'"$ACC": done

done

echo -e Done'\n'

# merging stringtie gtf files into a new gtf annotation file
echo -e Merging created gtf files and re-running StringTie...'\n'

ls "$OUTDIR"/*gtf > 04_mapping_counts/input_stringtiemerge.ls

stringtie --merge -G 00_input/mgal_genome/GCA.900618805.1_mgal_genomic.gff -p 15 -o "$OUTDIR"/GCA.900618805.1_mgal_genomic_stringtiemerged.gtf 04_mapping_counts/input_stringtiemerge.ls

# second round of stringtie against reference gff annotation file
echo Second round of StringTie \(vs merged gtf\)...

MERGED_GTF="$OUTDIR/GCA.900618805.1_mgal_genomic_stringtiemerged.gtf"

for i in 03_mapped_reads/*sortedByCoord.out.bam; do

        ACC="$(basename $i | awk -F "_" '{print $1}')"

        # run stringtie
	stringtie "$i" -e -G "$MERGED_GTF" -o "$OUTDIR"/"$ACC"_stringtie_round2.gtf -p 15 2> "$OUTDIR"/"$ACC"_stringtie_round2.log

        echo -e '\t'"$ACC": done

        # create input file for prepDE.py script
        echo -e "$ACC"'\t'"$OUTDIR"/"$ACC"_stringtie_round2.gtf >> "$OUTDIR"/input_prepDE.tsv

done

echo -e Done'\n'

# run prepDE.py to extract read counts
echo "Extracting read counts for DE analysis using prepDE.py..."

prepDE.py -i "$OUTDIR"/input_prepDE.tsv -l 99 -g "$OUTDIR"/gene_count_matrix.csv -t "$OUTDIR"/transcript_count_matrix.csv
