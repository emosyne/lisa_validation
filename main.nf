#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

include { lisa_percohort_devel } from './workflows/lisa_percohort_devel.nf'

workflow {
    lisa_percohort_devel()
}