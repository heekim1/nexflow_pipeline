nextflow.enable.dsl=2

process align {
  tag { (chrom == '') ? sampleName : sampleName + "-" + chrom }

  publishDir "${params.finalDir}/analysis/${sampleName}/align", mode: 'copy'

  input:
      val(sampleName)
      val(chrom)
      path(sampleFastq)
      val(reference)

  output:
      path "*.sorted.bam"
      path "*.sorted.bam.bai"

  script:
  """
      bwa mem \
      -R "@RG\\tID:${sampleName}\\tSM:bwa" \
      -v 1 \
      -t ${task.cpus} \
      '${reference}' \
      '${sampleFastq}' | \
      samtools view  -F256 -bh - > '${sampleName}.bam'

      samtools sort '${sampleName}.bam' > '${sampleName}.sorted.bam'
      samtools index '${sampleName}.sorted.bam'

  """
}
