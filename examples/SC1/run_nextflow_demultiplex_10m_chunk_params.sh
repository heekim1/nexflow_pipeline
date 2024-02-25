#!/bin/sh

module load jdk/11.0.2
#module load nextflow_latest/20.10.0
NEXTFLOW_VER=$(realpath ../../nextflow_ver/20.10.0)
OUTPUT_DIR=$(realpath ./output_demux_10m_chunk_params)
mkdir -p $OUTPUT_DIR

$NEXTFLOW_VER/nextflow run ../../nextflow/workflows/demultiplex.nf \
-profile cluster \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-work-dir $OUTPUT_DIR/work \
--runFastq "/sc1/groups/naa/DATA/SBX/fastq/poccrun/run_post_qc_sequences.fastq.gz" \
--regexpattern "/sc1/groups/naa/DATA/SBX/pipeline/resource/pocc_patterns_072722.txt" \
--chunkSize 10000000 \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt


