nextflow.enable.dsl=2

process splitBam {
  tag { sampleName }

  //publishDir "${params.finalDir}/analysis/${sampleName}/splitBam", mode: 'copy'

  input:
      val(sampleName)
      path(bam)
      path(bai)

  output:
      path "${sampleName}.*.sorted.bam"
      path "${sampleName}.*.sorted.bam.bai"
      
  script:
  """
      ##- split chr1-22,X,Y and M into chromosom.sorted.bam
      for f in `samtools idxstats ${bam} | grep ${params.chromPrefix} | awk '\$1 ~ /${params.chromPrefix}[0-9XYM]+\$/ {print \$1}'`;
      do
          samtools view -bh ${bam} \$f -o ${sampleName}.\$f.sorted.bam; 
          samtools index ${sampleName}.\$f.sorted.bam; 
      done

      ##- retrieve all scaffolds, which is others than chr1-22,X,Y and M, and merge into scaffolds.sorted.bam.
      samtools idxstats ${bam} | grep ${params.chromPrefix} | awk '\$1 !~ /${params.chromPrefix}[0-9XYM]+\$/ {print \$1"\t"0"\t"\$2}' > scaffolds.bed
      samtools view -hb ${bam} -L scaffolds.bed -o ${sampleName}.scaffolds.sorted.bam
      samtools index ${sampleName}.scaffolds.sorted.bam
      
  """
}
