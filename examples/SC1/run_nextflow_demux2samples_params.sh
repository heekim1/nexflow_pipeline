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
--runFastq "/sc1/groups/naa/DATA/SBX/fastq/poccrun/run_fastq_10k.fastq.gz" \
--regexpattern "/sc1/groups/naa/DATA/SBX/pipeline/resource/pocc_patterns_072722.txt" \
--reference  "/sc1/groups/bfx-red/data/datainsights/reference_genome/human/hg38decoy/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa" \
--regions "/sc1/groups/naa/DATA/SBX/pipeline/resource/regions3.bed" \
--chunkSize 2500 \
--splitByChrom true \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt


