nextflow.enable.dsl=2

process cpuDemux {
  tag { "cpuDemux" + chunkNum }

  memory "${params.parameters.cpuDemux.memory}"
  cpus "${params.parameters.cpuDemux.cpus}"
  clusterOptions { "${params.parameters.cpuDemux.clusterOptions}" }
  container "${params.cpuDemuxVersion}"

  publishDir "${params.finalDir}/analysis", mode: 'copy'

  input:
      val(chunkNum)
      path(runFastq)
      path(regexpattern)

  output:
      path 'demux/**/*' 

  script:
  """
      mkdir -p demux

      cpu_demux.py \
      -i '${runFastq}' \
      -o ./demux \
      -p '${regexpattern}'

      input="${regexpattern}"
      while IFS= read -r line
      do
        if [[ \$line =~ ^#.* ]]; then
            sample_name=\${line:1}
            mkdir -p demux/\${sample_name}
            mv demux/*\${sample_name}.* demux/\${sample_name}/
        fi
      done < "\$input"
  """
}
