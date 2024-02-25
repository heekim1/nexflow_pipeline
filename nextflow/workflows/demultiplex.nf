#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { demux } from '../modules/demux' params( params )
include { cpuDemux } from '../modules/cpuDemux' params( params )
include { mergeSampleFq } from '../modules/mergeSampleFq' params( params )

workflow {

    runFastq = channel.fromPath(params.runFastq)
                      .splitFastq( by: params.chunkSize, file: true)
    runFastq.view()
    // tokenized fastq name ends with base_name.<chunk_number>.fastq. 
    // "-2" is the index of the chunk_number.
    chunkNum = runFastq.flatten() | map { f -> f.Name.tokenize('.')[-2] }
    regexpattern = channel.value(params.regexpattern)
    if (params.demuxTool == "sbx_index_tool")
        sampleFastq = demux(chunkNum, runFastq, regexpattern)
    else
        sampleFastq = cpuDemux(chunkNum, runFastq, regexpattern)

    sampleFastq.flatten()
               .map { fastq -> [fastq.simpleName, fastq] }.groupTuple()
               .multiMap { s, f ->
                          sampleName: s
                          chunks: f
                }.set{ result }

    demuxedFastqFile = mergeSampleFq(result.sampleName, result.chunks)
}

workflow.onComplete {
    println "Demultiplexing workflow started at: $workflow.start"
    println "Demultiplexing workflow completed at: $workflow.complete"
    println "Demultiplexing workflow duration: $workflow.duration"
    println "Exit status: $workflow.exitStatus"
    println "Error message: $workflow.errorMessage"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Demultiplexing workflow execution stopped with the following message: ${workflow.errorMessage}"
    println "Exit status of the task that caused the workflow execution to fail: ${workflow.exitStatus}"
    println "Detailed error of the task that caused the workflow execution to fail: ${workflow.exitStatus}"
}
