#!/bin/bash

bash wrapper_params.sh \
-p /sc1/groups/onco/Analysis/bioinfo_analyses/kimh89/SBX/bfx-ngs-AVC-225/pipeline/nextflow/params/SC1/params_demux2samples_10m_chunk_sc1.txt \
-seg single-sample \
-d /sc1/groups/onco/Analysis/bioinfo_analyses/kimh89/SBX/bfx-ngs-AVC-225/pipeline/examples/SC1/output_demux_10m_chunk_params/analysis/demux \
-o /sc1/groups/onco/Analysis/bioinfo_analyses/kimh89/SBX/bfx-ngs-AVC-225/pipeline/examples/SC1/output_demux_10m_chunk_params \
-e cluster
