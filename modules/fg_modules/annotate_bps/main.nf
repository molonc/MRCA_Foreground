process ANNOTATE_BPS {
    // errorStrategy 'ignore'

    label 'copy_output'
    tag "${sample_id}"
    label 'process_single'
    label 'process_long'
    container params.base_container


    input:
        tuple val(sample_id), path(reads_2), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path(breakpoint_functions)
        val tcn_bool
    output:
        tuple val(sample_id), path(reads_2), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path("${sample_id}_segs_bp_anno.csv.gz"), path("${sample_id}_segs_bp_anno_A.csv.gz"), path("${sample_id}_segs_bp_anno_B.csv.gz"), path(breakpoint_functions), emit: master

    script:
    """
    echo 'fasdfafsdf'
    Rscript ${projectDir}/modules/fg_modules/annotate_bps/annotate_bps.R "${reads_2}" "${sample_id}" "${tcn_bool}" "${breakpoint_functions}"
    """
}
