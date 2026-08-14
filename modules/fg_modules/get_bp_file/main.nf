process GET_BP_FILE {
    // errorStrategy 'ignore'
    label 'copy_output'
    tag "${library_id}-${sample_id}-${exp_con}"
    label 'process_high'

    conda "${moduleDir}/r44.yml"


    input:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads_2), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path(segs_bp_anno), path(segs_bp_anno_A), path(segs_bp_anno_B), path(breakpoint_functions)
    output:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads_2), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path(segs_bp_anno), path(segs_bp_anno_A), path(segs_bp_anno_B), path("${library_id}-${sample_id}_breakpoint_file.csv.gz"), path("${library_id}-${sample_id}_breakpoint_file_A.csv.gz"), path("${library_id}-${sample_id}_breakpoint_file_B.csv.gz"), path(breakpoint_functions), emit: master

    script:
    """
    Rscript ${projectDir}/modules/fg_modules/get_bp_file/get_bp_file.R "${segs_bp_anno}"   "${library_id}" "${sample_id}" "${breakpoint_functions}"
    Rscript ${projectDir}/modules/fg_modules/get_bp_file/get_bp_file.R "${segs_bp_anno_A}" "${library_id}" "${sample_id}" "${breakpoint_functions}" "A"
    Rscript ${projectDir}/modules/fg_modules/get_bp_file/get_bp_file.R "${segs_bp_anno_B}" "${library_id}" "${sample_id}" "${breakpoint_functions}" "B"
    """
}
