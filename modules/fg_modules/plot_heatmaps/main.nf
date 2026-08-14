
process PLOT_HEATMAPS {
    label 'copy_output'
    tag "${library_id}-${sample_id}-${exp_con}"
    label 'process_medium'

    conda "${moduleDir}/r44.yml"

    input:
        tuple val(sample_id), val(library_id), val(exp_con), path(hdp_segs_file), path(hdp_segs_A), path(hdp_segs_B), path(reads_2), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path(segs_bp_anno), path(segs_bp_anno_A), path(segs_bp_anno_B), path(bp_file), path(bp_file_A), path(bp_file_B), path(breakpoint_functions), path(ref_bins)

    output:
        path "CNV_ABSOLUTE_ALLELE_RATIO_${library_id}.png", optional: true
        path "CNV_ABSOLUTE_ALLELE_PHASE_${library_id}.png", optional: true
        path "CNV_ABSOLUTE_${library_id}.png"
        path "CNV_FOREGROUND_${library_id}.png"
        path "CNV_FOREGROUND_A_${library_id}.png", optional: true
        path "CNV_FOREGROUND_B_${library_id}.png", optional: true
        path "CNV_FOREGROUND_EDGES_${library_id}.png"
        path "CNV_FOREGROUND_EDGE_STATUS_A_${library_id}.png", optional: true
        path "CNV_FOREGROUND_EDGE_STATUS_B_${library_id}.png", optional: true
        path "CNV_FOREGROUND_EDGES_A_${library_id}.png", optional: true
        path "CNV_FOREGROUND_EDGES_B_${library_id}.png", optional: true
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
        --library_id ${library_id} \
        --sample_id ${sample_id} \
        --hdp_segs_A ${hdp_segs_A} \
        --hdp_segs_B ${hdp_segs_B} \
        --segs_bp_anno_A ${segs_bp_anno_A} \
        --segs_bp_anno_B ${segs_bp_anno_B} \
        --ref_bins ${ref_bins}

    """
}

