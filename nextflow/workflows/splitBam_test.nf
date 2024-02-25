#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { splitBam } from '../modules/splitBam' params( params )

workflow {

    sampleName = params.sampleName
    bam = channel.fromPath(params.bam)
    bai = channel.fromPath(params.bai)

    (bam, bai) = splitBam(sampleName, bam, bai)
    bam.view()
    chrom = bam.flatten().map { f -> f.Name.tokenize('.')[1] }
    chrom.view()

}


workflow.onComplete {
    println "Exit status: $workflow.exitStatus"
    println "Error message: $workflow.errorMessage"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Exit status of the task that caused the workflow execution to fail: ${workflow.exitStatus}"
    println "Detailed error of the task that caused the workflow execution to fail: ${workflow.exitStatus}"
}
