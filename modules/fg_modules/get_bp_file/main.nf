process GET_BP_FILE {
    // errorStrategy 'ignore'
    label 'copy_output'
    tag "${sample_id}"
    label 'process_high'

    container params.base_container


    input:
        tuple val(sample_id), path(reads_2), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path(segs_bp_anno), path(segs_bp_anno_A), path(segs_bp_anno_B), path(breakpoint_functions)
    output:
        tuple val(sample_id), path(reads_2), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path(segs_bp_anno), path(segs_bp_anno_A), path(segs_bp_anno_B), path("${sample_id}_breakpoint_file.csv.gz"), path("${sample_id}_breakpoint_file_A.csv.gz"), path("${sample_id}_breakpoint_file_B.csv.gz"), path(breakpoint_functions), emit: master

    script:
    """
    Rscript ${projectDir}/modules/fg_modules/get_bp_file/get_bp_file.R "${segs_bp_anno}"   "${sample_id}" "${breakpoint_functions}"
    Rscript ${projectDir}/modules/fg_modules/get_bp_file/get_bp_file.R "${segs_bp_anno_A}" "${sample_id}" "${breakpoint_functions}" "A"
    Rscript ${projectDir}/modules/fg_modules/get_bp_file/get_bp_file.R "${segs_bp_anno_B}" "${sample_id}" "${breakpoint_functions}" "B"
    """
}
