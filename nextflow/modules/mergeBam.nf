nextflow.enable.dsl=2

process mergeBam {
  tag { "mergeBam-" + sampleName }

  publishDir "${params.finalDir}/analysis/${sampleName}/mergeBam", mode: 'copy'

  input:
      val(sampleName)
      path("*")

  output:
      path "${sampleName}.merged.bam"
      path "${sampleName}.merged.bam.bai"

  script:
  """
      samtools merge -c ${sampleName}.merged.bam *.bam
      samtools index ${sampleName}.merged.bam
  """
}
