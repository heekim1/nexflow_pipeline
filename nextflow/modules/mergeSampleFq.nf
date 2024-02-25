nextflow.enable.dsl=2

process mergeSampleFq {
  tag { "mergeSampleFq-" + sampleName }

  publishDir "${params.finalDir}/analysis", mode: 'copy'

  input:
      val(sampleName)
      path("*")

  output:
      path "demux/${sampleName}/*"

  script:
  """
      mkdir -p demux/${sampleName}
      cat `find .  -maxdepth 1 -name "*.fq.gz"` > demux/${sampleName}/${sampleName}.demuxed.fq.gz 
  """
}
