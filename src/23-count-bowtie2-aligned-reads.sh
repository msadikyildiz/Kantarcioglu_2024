#!/bin/bash

PROJECT_PATH="/project/greencenter/Toprak_lab/shared/TolC-Mutagenesis"
EXP_ID="20230403"

INPUTDIR="$PROJECT_PATH/data/$EXP_ID/alignment/bowtie2-single"
OUTPUT="$PROJECT_PATH/dump/$EXP_ID/${EXP_ID}_aligned_read_counts_single.txt"

printf "" > $OUTPUT
for SAMPLE_PATH in $INPUTDIR/*.bam; do
    SAMPLE=${SAMPLE_PATH##*/}
    SAMPLE=${SAMPLE%.bam}
    printf "${SAMPLE}\t" >> $OUTPUT
    samtools idxstats $SAMPLE_PATH | grep "TolC-oxb14" | cut -f 3 | cat - >> $OUTPUT
done