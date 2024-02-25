nextflow.enable.dsl=2

process consensusGen {
  tag { (chrom == '') ? sampleName : sampleName + "-" + chrom }

  publishDir "${params.finalDir}/analysis/${sampleName}/consensusGen", mode: 'copy'

  input:
      val(sampleName)
      val(chrom)
      path(clusteredFastq)
      path(clusteredTxt)

  output:
      path "*.consensus.fastq"

  script:
  """
      cobbsalad \
      --input ${clusteredFastq}  \
      --consensus_output_file '${sampleName}.consensus.fastq' \
      -t ${task.cpus} \
      --mode ics \
      --algorithm 7 \
      --min_pass 1 \
      --discard_front 0 \
      --discard_back 0 \
      --adapter '' \
      --trim_consensus_left 2 \
      --trim_consensus_right 2

      if [[ "${params.splitByChrom}" == "true" && "${chrom}" != "" ]]; then
          mv ${sampleName}.consensus.fastq ${sampleName}.${chrom}.consensus.fastq
      fi
      
  """
}
