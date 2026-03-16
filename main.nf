#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HEADTEXT } from './modules/local/headtext.nf'

workflow {

    input_ch = Channel.fromPath(params.input)

    HEADTEXT(input_ch)
}
