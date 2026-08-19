process GET_DF {
    label 'symlink_output'
    tag "${sample_id}"
    label 'process_single'

    container params.base_container

    input:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path(dummy), path(cen_info), path(breakpoint_functions)
        val tcn_bool

    output:
        tuple val(sample_id), path("${sample_id}_reads_final.csv.gz"), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path(breakpoint_functions), emit: master  // Note reads here is now the getdf output

    script:
    """
    echo 'ffff'
    Rscript ${projectDir}/modules/fg_modules/get_df/assemble_full_df.R \
        --segs "${segs}" \
        --dummy "${dummy}" \
        --output "${sample_id}_reads_final.csv.gz" \
        --tcn_bool "${tcn_bool}" \
        --breakpoint_functions "${breakpoint_functions}" \
        --hscn_path "${hscn}" \
        --cen_info "${cen_info}"
    """
}