process RUN_BAYESIAN_ANALYSIS_PER_GENE {
    label 'process_single'

    tag "${gene}"

    input:
    tuple val(gene), path(gene_partition)
    path annotations_file

    output:
    path("run_bayesian_analysis.done"), emit: done

    script:
    def resultsDir = "${params.bucket_uri}/${params.bayes_results_dir}"
    """
    python3 "${params.bayes_script}" \
        --analysis_dir "${gene_partition}" \
        --annotations_file "${annotations_file}" \
        --output_dir "${resultsDir}" \
        --gene "${gene}" \
        --h1 ${params.h1}

    touch run_bayesian_analysis.done
    """
}
