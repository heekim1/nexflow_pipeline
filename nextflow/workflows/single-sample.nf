#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { align } from '../modules/align' params( params )
include { splitBam } from '../modules/splitBam' params( params )
include { umiCluster } from '../modules/umiCluster' params( params )
include { consensusGen } from '../modules/consensusGen' params( params )
include { postAln } from '../modules/postAln' params( params )
include { mergeBam } from '../modules/mergeBam' params( params )

def getChrom(obj){
    if(params.splitByChrom == true)
        return obj.flatten().map { f -> f.Name.tokenize('.')[1] }
    else
        return ""
}

workflow {

    sampleName = params.sampleName
    chrom = ""
    sampleFastq = channel.fromPath(params.sampleFastq)
    reference = channel.value(params.reference)
    regions = channel.value(params.regions)

    (bam, bai) = align(sampleName, chrom, sampleFastq, reference)

    if( params.splitByChrom ) {
        (bam, bai) = splitBam(sampleName, bam, bai)
    }
    (clusteredFastq, clusteredTxt) = umiCluster(sampleName, getChrom(bam), bam.flatten(), bai.flatten(), reference, regions)

    (consensusFastq) = consensusGen(sampleName, getChrom(clusteredFastq), clusteredFastq, clusteredTxt)

    (postAlnBam, postAlnBai) = postAln(sampleName, getChrom(consensusFastq), consensusFastq, reference)

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
    println "Single-sample workflow started at: $workflow.start"
    println "Single-sample  workflow completed at: $workflow.complete"
    println "Single-sample workflow duration: $workflow.duration"
    println "Exit status: $workflow.exitStatus"
    println "Error message: $workflow.errorMessage"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Single-sample workflow execution stopped with the following message: ${workflow.errorMessage}"
    println "Exit status of the task that caused the workflow execution to fail: ${workflow.exitStatus}"
    println "Detailed error of the task that caused the workflow execution to fail: ${workflow.exitStatus}"
}
