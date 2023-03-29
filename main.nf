#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

include { lisa_validation } from './workflows/lisa_validation.nf'

workflow {
    lisa_validation()
}