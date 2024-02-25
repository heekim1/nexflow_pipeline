# Draft for Pipeline Structure and Usage
- **[projectDir](./nextflow)** -- directory where the main workflow script is located.
- **[modules](./nextflow/modules)** -- directory where module scripts are located. Reuseable modules are in here. A module can consists of process definitions or functions.
- **parameterfile** -- file containing parameters to pass to a script using the `-params-file` or `-c` option. This will be a yaml file containing the following information:
  
  - [params](./nextflow/params) : Parameters and configurations that will be accessible in the pipeline script through the option `-c`
  - [yaml](./nextflow/yamal)  : to define parameters that will be accessible in the pipeline script through the option `-params-file`
  
- **[nextflow.config](./nextflow/workflows/nextflow.config)** -- file containing configuration settings for a workflow; it is stored in same dir as a workflow script.

  - `nextflow.config` includes profiles to predefine multiple configurations for different execution platforms.
  - `nextflow.config` includes the withLabel selectors that allow the configuration of all processes annotated with a label directive
  
  - Example

  ```
  profiles {
      aws {
         TBD
      }
    
      cluster {
        docker.enabled = true
        process.executor = 'sge'
        process.queue = 'all.q'
        process.penv = 'smp'
        process.clusterOptions = { "-l h_vmem=${task.ext.vmem} -S /bin/bash -l docker_version=new -P onco" }
        docker.runOptions = "-u=\$UID --rm -v /sc1:/sc1"
      }
    
      local {
        docker.enabled = true
        process.executor = 'local'
        process.clusterOptions = ''
        docker.runOptions = "-u=\$UID --rm -v /sc1:/sc1"
      }
  }
  
  process {
      withName:  demux{
         memory = '32 GB'
         cpus = 8
         container = "${params.demux_sbx_index_toolsVersion}"
         ext {
              vmem = '4G'
         }
      }
  }
  
  process {
      withName:  cpuDemux{
         memory = '32 GB'
         cpus = 8
         container = "${params.cpuDemuxVersion}"
         ext {
              vmem = '4G'
         }
      }
  }
  ```

- **[versions.config](./nextflow/workflows/versions.config)** -- file to define the docker container versions.
  
## SBX Pipeline Framework 

```
├── modules                           (High Level Process Definition)
│   ├── align.nf
│   ├── consensusGen.nf
│   ├── demux.nf
│   ├── postAln.nf
│   └── umiCluster.nf
├── nextflow.config                  (Profiles and Process Configuration )
├── demultiplex.nf                   (Workflow for Demux )
├── single-sample.nf                 (Workflow for Samples including Aln, Cluster, Consensus)
├── demux2samples.nf                 (Workflow for both Demux and Sample)
└── versions.config                  (Docker Container Versions)

```

## Workflow

A workflow is a component which encloses the invocation of one or more processes and operators:

### Demultiplex (demultiplex.nf)
```
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { demux } from '../modules/demux' params( params )
include { cpuDemux } from '../modules/cpuDemux' params( params )

workflow {

    runFastq = channel.fromPath(params.runFastq)
    regexpattern = channel.fromPath(params.regexpattern)
    if (params.demuxTool == "sbx_index_tool")
        sampleFastqFiles = demux(runFastq, regexpattern)
    else
        sampleFastqFiles = cpuDemux(runFastq, regexpattern)

}

```
### Single-Sample (single-sample.nf)
```
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { align } from '../modules/align' params( params )
include { umiCluster } from '../modules/umiCluster' params( params )
include { consensusGen } from '../modules/consensusGen' params( params )
include { postAln } from '../modules/postAln' params( params )

workflow {

    sampleName = params.sampleName
    sampleFastq = channel.fromPath(params.sampleFastq)
    reference = channel.fromPath(params.reference)
    regions = channel.fromPath(params.regions)

    (bam, bai) = align(sampleName, sampleFastq, reference)

    (clusteredFastq, clusteredTxt) = umiCluster(sampleName, bam, bai, reference, regions)

    (consensusFastq) = consensusGen(sampleName, clusteredFastq, clusteredTxt)

    (postAlnBam, postAlnBai) = postAln(sampleName, consensusFastq, reference)

}
```

### dexmux2samples (demultiplex + singl-sample)
```
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { demux } from '../modules/demux' params( params )
include { cpuDemux } from '../modules/cpuDemux' params( params )
include { align } from '../modules/align' params( params )
include { umiCluster } from '../modules/umiCluster' params( params )
include { consensusGen } from '../modules/consensusGen' params( params )
include { postAln } from '../modules/postAln' params( params )

workflow {

    runFastq = channel.fromPath(params.runFastq)
    regexpattern = channel.fromPath(params.regexpattern)
    if (params.demuxTool == "sbx_index_tool")
        sampleFastqFiles = demux(runFastq, regexpattern)
    else
        sampleFastqFiles = cpuDemux(runFastq, regexpattern)

    sampleName = sampleFastqFiles.flatten() | map { fastq -> fastq.SimpleName }
    reference = channel.value(params.reference)
    regions = channel.value(params.regions)
    (bam, bai) = align(sampleName, sampleFastqFiles.flatten(), reference)

    sampleName = bam.flatten() | map { bam -> bam.SimpleName }
    (clusteredFastq, clusteredTxt) = umiCluster(sampleName, bam, bai, reference, regions)

    sampleName = clusteredFastq.flatten() | map { fastq -> fastq.SimpleName }
    (consensusFastq) = consensusGen(sampleName, clusteredFastq, clusteredTxt)

    sampleName = consensusFastq.flatten() | map { fastq -> fastq.SimpleName }
    (postAlnBam, postAlnBai) = postAln(sampleName, consensusFastq, reference)

}
```

## Run SBX workflows (demultiplex.nf, single-sample.nf, and demux2samples.nf)

Demultiplex - demultiplex.nf
Single-Sample - single-sample.nf
Demultiplex + Samples - demux2samples.nf

### Demultiplex

#### Using params
##### Command
```
$NEXTFLOW_VER/nextflow run ../nextflow/workflows/demultiplex.nf \
-profile local \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-work-dir $OUTPUT_DIR/work \
--runFastq "/sc1/groups/naa/DATA/SBX/fastq/poccrun/run_fastq_10k.fastq.gz" \
--regexpattern "/cloud_0722_pocc/pocc_patterns_072722.txt" \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt

or

$NEXTFLOW_VER/nextflow run ../nextflow/workflows/demultiplex.nf \
-profile local \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-work-dir $OUTPUT_DIR/work \
-c params_demultiplexing_10k.txt
```
##### Params (for the option -c)
```
params.runFastq = "./run_fastq_10k.fastq.gz"
params.regexpattern = "./pocc_patterns_072722.txt"
params.finalDir = "./examples/out_10k"
```

#### Using Yaml
##### Command
```
$NEXTFLOW_VER/nextflow run ../nextflow/workflows/demultiplex.nf \
-profile local \
-params-file ../nextflow/yaml/SC1/params_demux_10k_sc1.yaml \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-w $OUTPUT_DIR/work \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt
```
##### Yaml (for the option -params-file )
```
runFastq: './run_fastq_10k.fastq.gz'
regexpattern: '/cloud_0722_pocc/pocc_patterns_072722.txt'
finalDir: './output_all'

```
### Single-Sample
#### Using Params
##### Command
```
$NEXTFLOW_VER/nextflow run ../nextflow/workflows/single-sample.nf \
-profile local \
-params-file ../nextflow/yaml/SC1/params_single_sample_10k_sc1.yaml \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-w $OUTPUT_DIR/work \
--sampleNmae "cfDNA6" \
--sampleFastq "./cfDNA6.run_fastq_10k.umis.demux.fq.gz" \
--reference  "./GRCh38.fa" \
--regions "/NAA0061/regions3.bed" \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt

or

$NEXTFLOW_VER/nextflow run ../nextflow/workflows/single-sample.nf \
-profile local \
-params-file ../nextflow/yaml/SC1/params_single_sample_10k_sc1.yaml \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-w $OUTPUT_DIR/work \
-c params_single-sample_10k.txt.txt
```

##### Params (for the option -c)
```
params.runFastq = "./run_fastq_10k.fastq.gz"
params.regexpattern = "./pipeline/examples/params/pocc_patterns_072722.txt"
params.finalDir = "./pipeline/examples/out_10k"

params.sampleName="cfDNA6"
params.sampleFastq = "./cfDNA6.run_fastq_10k.umis.demux.fq.gz"
params.reference = "./GRCh38.fa"
params.bed = "/NAA0061/regions3.bed"
```
#### Using Yaml
##### Command
```
$NEXTFLOW_VER/nextflow run ../nextflow/workflows/single-sample.nf \
-profile local \
-params-file ../nextflow/yaml/SC1/params_single_sample_10k_sc1.yaml \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-w $OUTPUT_DIR/work \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt

```
##### Yaml (for the option -params-file )
```
runFastq: '/sc1/groups/naa/DATA/SBX/fastq/poccrun/run_fastq_10k.fastq.gz'
regexpattern: './pocc_patterns_072722.txt'
finalDir: './output_all'
sampleName: 'cfDNA6'
sampleFastq: './cfDNA6.run_fastq_10k.umis.demux.fq.gz'
reference: './GRCh38.fa'
regions: '/NAA0061/regions3.bed'
```

### Demux2Samples
#### Using Params
##### Command
```
$NEXTFLOW_VER/nextflow run ../nextflow/workflows/demux2samples.nf \
-profile local \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-work-dir $OUTPUT_DIR/work \
--runFastq "/sc1/groups/naa/DATA/SBX/fastq/poccrun/run_fastq_10k.fastq.gz" \
--regexpattern "/sc1/groups/bfx-red/projects/Carlo/pipeline_runs/users/shannog2/cloud_0722_pocc/pocc_patterns_072722.txt" \
--reference  "/sc1/groups/bfx-red/data/datainsights/reference_genome/human/hg38decoy/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa" \
--regions "/sc1/groups/naa/People/padhukab/Analyses/NAA/NAA0061/regions3.bed" \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt

or

$NEXTFLOW_VER/nextflow run ../nextflow/workflows/demux2samples.nf \
-profile cluster \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-w $OUTPUT_DIR/work \
-c params_single-sample_10k.txt.txt

```
##### params (for the option -c)
```
params.runFastq = "/run_fastq_10k.fastq.gz"
params.regexpattern = "pocc_patterns_072722.txt"
params.finalDir = "./pipeline/examples/out_10k"

params.sampleName="cfDNA6"
params.sampleFastq = "./cfDNA6.run_fastq_10k.umis.demux.fq.gz"
params.reference = "/GRCh38.fa"
params.bed = "/NAA0061/regions3.bed"
```
#### Using Yaml
##### Command
```
$NEXTFLOW_VER/nextflow run ../nextflow/workflows/demux2samples.nf \
-profile local \
-params-file ../nextflow/yaml/SC1/params_demux_10k_sc1.yaml \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-w $OUTPUT_DIR/work \
--finalDir $OUTPUT_DIR UTPUT_DIR/log.txt | tee -a  $OUTPUT_DIR/log.txt

```
##### Yaml (for the option -params-file)
```
runFastq: '/run_fastq_10k.fastq.gz'
regexpattern: '/pocc_patterns_072722.txt'
finalDir: './output_all'
reference: '/GRCh38.fa' \
regions: '/NAA0061/regions3.bed'
```

### Wrapper 
1. Provide easy-to-use command line-based tool to utilize demultiplex.nf and single-sample.nf. 
2. Provide organized log directories for each sample which helps users troubleshoot.
#### Usage
```
(base) [kimh89@lb100hmdev examples]$ bash wrapper_params.sh -h

The following have been reloaded with a version change:
  1) jdk/8u92-b14 => jdk/11.0.2


Program: Wrapper(CLI) for the bfx-ngs pipeline
Usage:
        bash run_nextflow_demux_single-sample.sh [options]
Options: -p|--params         parameter file
         -o|--output-dir     output directory
         -e|--executor       executor (local|cluster|aws)
         -seg|--segment      semgment (demux: demux only | single-sample: sample process only). (Default=12)
         -d|--demux-dir      demultiplexed directory only if -seg=single-sample. (Default=/analysis/demux)
         -h|--help

```
#### Command for demultiplexing+single-sample
```
bash wrapper_params.sh \
-p /pipeline/nextflow/params/SC1/params_demux2samples_10k_sc1.txt \
-o /pipeline/examples/output_wrapper \
-e local
```
#### Params for demuxed single-sample(s)
```
bash wrapper_params.sh \
-p /params_demux2samples_10k_sc1.txt \
-o /pipeline/examples/output_samples \
-e local \
-d /pipeline/examples/output_wrapper/analysis/demux \
-seg single-sample
```

### Single Process run (Unit Test)

1. Create a workflow, for example, umiCluster.nf. 

```
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { umiCluster } from '../modules/umiCluster' params( params )

workflow {

    sampleName = params.sampleName
    bam = channel.fromPath(params.bam)
    bai = channel.fromPath(params.bai)
    reference = channel.fromPath(params.reference)
    regions = channel.fromPath(params.regions)

    (clusteredFastq, clusteredTxt) = umiCluster(sampleName, bam, bai, reference, regions)

}
```
2. Run with params
```
$NEXTFLOW_VER/nextflow run ../nextflow/workflows/umiCluster.nf \
-profile local \
-ansi-log false \
-with-report $OUTPUT_DIR/report.html \
-with-timeline $OUTPUT_DIR/timeline.html \
-with-trace $OUTPUT_DIR/trace.txt \
-work-dir $OUTPUT_DIR/work \
--sampleName "cfDNA6" \
--bam "/cfDNA6.sorted.bam" \
--bai "/cfDNA6.sorted.bam.bai" \
--reference  "GRCh38.fa" \
--regions "/NAA0061/regions3.bed" \
--finalDir $OUTPUT_DIR | tee -a  $OUTPUT_DIR/log.txt

```


## Repos

https://github.com/Roche-DIA-RDS-CSI/bfx-ngs
https://github.com/Roche-DIA-RDS-CSI/sbx-algo-rnd-cloud-pipeline
