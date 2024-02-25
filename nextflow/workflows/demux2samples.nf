#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { demux } from '../modules/demux' params( params )
include { cpuDemux } from '../modules/cpuDemux' params( params )
include { mergeSampleFq } from '../modules/mergeSampleFq' params( params )
include { align } from '../modules/align' params( params )
include { splitBam } from '../modules/splitBam' params( params )
include { umiCluster } from '../modules/umiCluster' params( params )
include { consensusGen } from '../modules/consensusGen' params( params )
include { postAln } from '../modules/postAln' params( params )
include { mergeBam } from '../modules/mergeBam' params( params )

def getSampleName(obj){
    return obj.flatten() | map { f -> f.SimpleName }
}

def getChrom(obj){
    if(params.splitByChrom)
        return obj.flatten().map { f -> f.Name.tokenize('.')[1] }
    else
        return ""
}

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
        sampleFastq = cpuDemux(runFastq, regexpattern)

    sampleFastq.flatten()
               .map { fastq -> [fastq.simpleName, fastq] }.groupTuple()
               .multiMap { s, f ->
                          sampleName: s
                          chunks: f
                }.set{ result }

    demuxedFastqFile = mergeSampleFq(result.sampleName, result.chunks)

    sampleName = demuxedFastqFile.flatten() | map { fastq -> fastq.SimpleName }
    chrom = ""
    reference = channel.value(params.reference)
    regions = channel.value(params.regions)
    (bam, bai) = align(getSampleName(demuxedFastqFile), chrom, demuxedFastqFile.flatten(), reference)

    if( params.splitByChrom) {
        (bam, bai) = splitBam(sampleName, bam, bai)
    }
    
    (clusteredFastq, clusteredTxt) = umiCluster(getSampleName(bam), getChrom(bam), bam.flatten(), bai.flatten(), reference, regions)

    (consensusFastq) = consensusGen(getSampleName(clusteredFastq), getChrom(clusteredFastq), clusteredFastq, clusteredTxt)

    (postAlnBam, postAlnBai) = postAln(getSampleName(consensusFastq), getChrom(consensusFastq), consensusFastq, reference)

    // map.groupTuple() returns list of ["sampleName","Sample's full path"]
    // multipMap returns two lists. One is a list of the sampleName. The other is a list of the full paths. 
    if( params.splitByChrom ) {
        postAlnBam.flatten()
                   .map{ f -> [f.simpleName, f] }.groupTuple()
                   .multiMap { s, f ->
                               sampleName:s
                               chunks:f
                   }.set{ result }

        (postAlnBam, postAlnBai) = mergeBam(result.sampleName, result.chunks)
    }
 
}

workflow.onComplete {
    println "demux2samples workflow started at: $workflow.start"
    println "demux2samples workflow completed at: $workflow.complete"
    println "demux2samples workflow duration: $workflow.duration"
    println "Exit status: $workflow.exitStatus"
    println "Error message: $workflow.errorMessage"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "demux2samples workflow execution stopped with the following message: ${workflow.errorMessage}"
    println "Exit status of the task that caused the workflow execution to fail: ${workflow.exitStatus}"
    println "Detailed error of the task that caused the workflow execution to fail: ${workflow.exitStatus}"
}
