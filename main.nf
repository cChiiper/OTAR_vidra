#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RUN_BAYESIAN_ANALYSIS_PER_GENE } from './modules/shared/run_bayesian_analysis_per_gene.nf'

def joinPath(String base, String child) {
    "${base.toString().replaceAll('/+$', '')}/${child.toString().replaceAll('^/+', '')}"
}

def partitionGene(def pathLike) {
    pathLike
        .toString()
        .replaceAll('/+$', '')
        .tokenize('/')
        .last()
        .replaceFirst('^as_gene=', '')
}

workflow {
    if (!params.bucket_uri) {
        error "Missing required parameter: --bucket_uri"
    }

    if (!params.analysis_ready_dir) {
        error "Missing required parameter: --analysis_ready_dir"
    }

    if (!params.annotations_dir) {
        error "Missing required parameter: --annotations_dir"
    }

    def analysisRoot = joinPath(params.bucket_uri, params.analysis_ready_dir)
    def annotationsFile = joinPath(
        joinPath(params.bucket_uri, params.annotations_dir),
        params.annotations_output_name
    )

    genePartitions = Channel
        .fromPath("${analysisRoot}/as_gene=*", checkIfExists: true, type: 'dir')
        .map { partitionPath -> tuple(partitionGene(partitionPath), partitionPath) }

    if (params.genes) {
        requestedGenes = Channel
            .fromPath(params.genes, checkIfExists: true)
            .splitText()
            .map { it.trim() }
            .filter { it }
            .unique()
            .map { gene -> tuple(gene, true) }

        genePartitions = genePartitions
            .join(requestedGenes)
            .map { gene, partitionPath, _ -> tuple(gene, partitionPath) }
    }

    genePartitions = genePartitions.ifEmpty {
        error "No gene partitions found under ${analysisRoot}"
    }

    annotations = Channel.fromPath(annotationsFile, checkIfExists: true).first()

    RUN_BAYESIAN_ANALYSIS_PER_GENE(genePartitions, annotations)
}
