process HEADTEXT {

    tag "$input_file"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path input_file

    output:
    path "head.txt"

    script:
    """
    python ${projectDir}/tools/headtext.py \
        --input "${input_file}" \
        --output head.txt \
        --n_lines ${params.n_lines}
    """
}
