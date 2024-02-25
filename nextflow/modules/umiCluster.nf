nextflow.enable.dsl=2

process umiCluster {
  tag { (chrom == '') ? sampleName : sampleName + "-" + chrom }

  publishDir "${params.finalDir}/analysis/${sampleName}/umiCluster", mode: 'copy'

  input:
      val(sampleName)
      val(chrom)
      path(bam)
      path(bai)
      val(reference)
      path(bed)

  output:
      path "*.clusters.fastq"
      path "*.clusters.txt"

  script:
  """
      uidclusterer \
      -t ${task.cpus} \
      --min_barcode_family_size 1 \
      --max_barcode_family_size 50 \
      --no_uid_clustering 1 \
      --consensus_chunk_size=20 \
      --positional_dedup_cutoff=2 \
      --alignment_flag_filter 2308 \
      --min_mapq 4 \
      --anchor_by_left \
      --disable_consensus \
      --trim \
      --umi_order 2 \
      --cutoff_for_good_cluster 5 \
      --umi_type randomer \
      --umi_strand 3 \
      --umi_distance_type hamming \
      -m '${bed}' \
      --ref ${reference} \
      -b  '${bam}' \
      --clustering_result_file '${sampleName}.clusters.txt' \
      --out_bam_file '${sampleName}.clusters.bam' > '${sampleName}.clusters.fastq'

      if [[ "${params.splitByChrom}" == "true" && "${chrom}" != "" ]]; then
          mv ${sampleName}.clusters.txt ${sampleName}.${chrom}.clusters.txt
          mv ${sampleName}.clusters.fastq ${sampleName}.${chrom}.clusters.fastq
      fi
  """
}
