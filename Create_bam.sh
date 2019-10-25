#!/bin/bash

REF="control"
PATH_TO_PICARD=$1
PATH_TO_GATK=$2
FILE=$3
OUT_DIR=$4
TMP_DIR=$5
pp=$6

cd $TMP_DIR

if [ $pp == "p" ]
then
	bwa mem -t 10 ${OUT_DIR}/${REF}.fa ${OUT_DIR}/${FILE}_1.fastq ${OUT_DIR}/${FILE}_2.fastq \
    -R '@RG\tID:foo\tLB:bar\tPL:illumina\tSM:Sample1\tPU:L001' | \
    samtools view -Su | samtools sort | samtools rmdup - ${FILE}.sr.rm.bam
else
	bwa mem -t 10 ${OUT_DIR}/${REF}.fa ${OUT_DIR}/${FILE}.fastq \
    -R '@RG\tID:foo\tLB:bar\tPL:illumina\tSM:Sample1\tPU:L001' | \
    samtools view -Su | samtools sort | samtools rmdup -s - ${FILE}.sr.rm.bam
fi

samtools index ${FILE}.sr.rm.bam

java -jar ${PATH_TO_GATK}/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-R ${OUT_DIR}/${REF}.fa -I ${FILE}.sr.rm.bam -o ${FILE}_IndelRealigner.intervals

java -jar ${PATH_TO_GATK}/GenomeAnalysisTK.jar --filter_bases_not_stored -T IndelRealigner \
-R ${OUT_DIR}/${REF}.fa -I ${FILE}.sr.rm.bam -targetIntervals ${FILE}_IndelRealigner.intervals \
-o ${OUT_DIR}/${FILE}.final.bam

samtools flagstat ${OUT_DIR}/${FILE}.final.bam