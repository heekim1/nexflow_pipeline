#!/bin/sh

module load jdk/11.0.2
#module load nextflow_latest/20.10.0
NEXTFLOW_VER=$(realpath ../../nextflow_ver/20.10.0)
OUTPUT_DIR=$(realpath ./output_single_sample_10k_params)
mkdir -p $OUTPUT_DIR

$NEXTFLOW_VER/nextflow run ../../nextflow/workflows/single-sample.nf \
-profile local \
-params-file ../../nextflow/yaml/SC1/params_single_sample_10k_sc1.yaml \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-w $OUTPUT_DIR/work \
--sampleNmae "cfDNA6" \
--sampleFastq "/sc1/groups/onco/Analysis/bioinfo_analyses/kimh89/SBX/bfx-ngs-AVC-225/pipeline/examples/SC1/output_demux_10k_params/analysis/demux/cfDNA6/cfDNA6.demuxed.fq.gz" \
--reference  "/sc1/groups/bfx-red/data/datainsights/reference_genome/human/hg38decoy/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa" \
--regions "/sc1/groups/naa/DATA/SBX/pipeline/resource/regions3.bed" \
--splitByChrom true \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt

