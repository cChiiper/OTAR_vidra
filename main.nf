#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HEADTEXT } from './modules/local/headtext.nf'

workflow {

    if (!params.input) {
        error "Missing required parameter: --input"
    }

    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .map { input_file ->
            tuple([id: input_file.baseName], input_file)
        }

    ch_headtext_script = Channel.value(file("${projectDir}/tools/headtext.py"))

    HEADTEXT(ch_input, ch_headtext_script)
}
