process ANNOTATE_BPS {
    // errorStrategy 'ignore'

    label 'copy_output'
    tag "${library_id}-${sample_id}-${exp_con}"
    label 'process_single'
    label 'process_long'
    conda "${moduleDir}/r44.yml"


    input:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads_2), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path(breakpoint_functions)
        val tcn_bool
    output:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads_2), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path("${library_id}_segs_bp_anno.csv.gz"), path("${library_id}_segs_bp_anno_A.csv.gz"), path("${library_id}_segs_bp_anno_B.csv.gz"), path(breakpoint_functions), emit: master

    script:
    """
    echo 'fasdfafsdf'
    Rscript ${projectDir}/modules/fg_modules/annotate_bps/annotate_bps.R "${reads_2}" "${library_id}" "${tcn_bool}" "${breakpoint_functions}"
    """
}
