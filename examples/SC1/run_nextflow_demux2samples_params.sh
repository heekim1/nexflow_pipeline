#!/bin/sh

module load jdk/11.0.2
#module load nextflow_latest/20.10.0
NEXTFLOW_VER=$(realpath ../../nextflow_ver/20.10.0)
OUTPUT_DIR=$(realpath ./output_demux2samples_10k_params)
mkdir -p $OUTPUT_DIR

$NEXTFLOW_VER/nextflow run ../../nextflow/workflows/demux2samples.nf \
-profile cluster \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-work-dir $OUTPUT_DIR/work \
--runFastq "/run_fastq_10k.fastq.gz" \
--regexpattern "/pocc_patterns_072722.txt" \
--reference  "/GRCh38.fa" \
--regions "/resource/regions3.bed" \
--chunkSize 2500 \
--splitByChrom true \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt


