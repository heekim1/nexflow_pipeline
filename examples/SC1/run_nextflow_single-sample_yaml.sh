#!/bin/sh

module load jdk/11.0.2
#module load nextflow_latest/20.10.0
NEXTFLOW_VER=$(realpath ../../nextflow_ver/20.10.0)
OUTPUT_DIR=$(realpath ./output_single_sample_10k_yaml)
mkdir -p $OUTPUT_DIR

$NEXTFLOW_VER/nextflow run ../../nextflow/workflows/single-sample.nf \
-profile local \
-params-file ../../nextflow/yaml/SC1/params_single_sample_10k_sc1.yaml \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-w $OUTPUT_DIR/work \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt

