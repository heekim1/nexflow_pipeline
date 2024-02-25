nextflow.enable.dsl=2

process demux {
  tag { "demux-" + chunkNum }

  input:
      val(chunkNum)
      path(runFastq)
      path(regexpattern)

  output:
      path "demux/**/*"

  script:
  """
      mkdir -p demux

      sbx_index_tools.py \
      -i  '${runFastq}' \
      --output_dir ./demux \
      -p '${regexpattern}' \
      --demux both \
      --match BESTMATCH \
      --expected_5prime_elements runway sid5 umi5 anchor5 \
      --expected_3prime_elements insert anchor3 umi3 sid3 end \
      --annotate_elements sid5 sid3 umi5 umi3 \
      --min_insert_length 3

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
