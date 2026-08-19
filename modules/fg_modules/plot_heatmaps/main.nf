
process PLOT_HEATMAPS {
    label 'copy_output'
    tag "${sample_id}"
    label 'process_medium'

    container params.base_container

    input:
        tuple val(sample_id), path(hdp_segs_file), path(hdp_segs_A), path(hdp_segs_B), path(reads_2), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path(segs_bp_anno), path(segs_bp_anno_A), path(segs_bp_anno_B), path(bp_file), path(bp_file_A), path(bp_file_B), path(breakpoint_functions), path(ref_bins)

    output:
        path "CNV_ABSOLUTE_ALLELE_RATIO_${sample_id}.png", optional: true
        path "CNV_ABSOLUTE_ALLELE_PHASE_${sample_id}.png", optional: true
        path "CNV_ABSOLUTE_${sample_id}.png"
        path "CNV_FOREGROUND_${sample_id}.png"
        path "CNV_FOREGROUND_A_${sample_id}.png", optional: true
        path "CNV_FOREGROUND_B_${sample_id}.png", optional: true
        path "CNV_FOREGROUND_EDGES_${sample_id}.png"
        path "CNV_FOREGROUND_EDGE_STATUS_A_${sample_id}.png", optional: true
        path "CNV_FOREGROUND_EDGE_STATUS_B_${sample_id}.png", optional: true
        path "CNV_FOREGROUND_EDGES_A_${sample_id}.png", optional: true
        path "CNV_FOREGROUND_EDGES_B_${sample_id}.png", optional: true
        path "breakpoint_cells/CNV_FOREGROUND_EDGES_*"
    script:
    """
    echo ''
    Rscript ${projectDir}/modules/fg_modules/plot_heatmaps/plot_heatmaps.R \
        --reads ${reads_2} \
        --metrics ${metrics} \
        --tree ${tree} \
        --breakpoints ${bp_file} \
        --annotated_segs ${segs_bp_anno} \
        --library_id ${sample_id} \
        --sample_id ${sample_id} \
        --hdp_segs_A ${hdp_segs_A} \
        --hdp_segs_B ${hdp_segs_B} \
        --segs_bp_anno_A ${segs_bp_anno_A} \
        --segs_bp_anno_B ${segs_bp_anno_B} \
        --ref_bins ${ref_bins}

    """
}

