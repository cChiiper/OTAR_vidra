process HEADTEXT {
    label 'process_single'

    tag "${meta.id}"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.id}.head.txt"), emit: txt

    script:
    """
    python3 "${params.headtext_script}" \
        --input "${input_file}" \
        --output "${meta.id}.head.txt" \
        --n_lines ${params.n_lines}
    """
}
