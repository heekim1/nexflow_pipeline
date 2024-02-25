#!/bin/sh

module load jdk/11.0.2
#module load nextflow_latest/20.10.0
NEXTFLOW_VER=$(realpath ../../nextflow_ver/20.10.0)
OUTPUT_DIR=$(realpath ./output_spliBam_10K_params)
mkdir -p $OUTPUT_DIR

$NEXTFLOW_VER/nextflow run ../../nextflow/workflows/splitBam_test.nf \
-profile cluster \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-work-dir $OUTPUT_DIR/work \
--sampleName "cfDNA6" \
--bam "/sc1/groups/onco/Analysis/bioinfo_analyses/kimh89/SBX/bfx-ngs-AVC-225/pipeline/examples/SC1/output_single_sample_10k_params/analysis/cfDNA6/align/cfDNA6.sorted.bam" \
--bai "/sc1/groups/onco/Analysis/bioinfo_analyses/kimh89/SBX/bfx-ngs-AVC-225/pipeline/examples/SC1/output_single_sample_10k_params/analysis/cfDNA6/align/cfDNA6.sorted.bam.bai" \
--chromPrefix "chr" \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt


