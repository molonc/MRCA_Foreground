process NND_PUNISHMENT_CALC {
    label 'symlink_output'
    tag "${sample_id}"
    label 'process_single'

    container params.base_container

    input:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn)
    output:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn), path("${sample_id}_all_cells_nnd.csv"), emit: master

        path "${sample_id}_all_cells_nnd.csv"
        path "${sample_id}_beta_parameters.csv"
        path "${sample_id}_divergent_cells_heatmap.png"
        path "${sample_id}_divergent_cells_nnd.csv"
        path "${sample_id}_nnd_beta_qqplot.png"
        path "${sample_id}_nnd_density.png"
        path "${sample_id}_single_cell_copy_divergent_cells.png"

    script:
    """
    Rscript ${projectDir}/modules/fg_modules/nnd/nnd_beta_divergence_punishment.R \
        --hscn ${hscn} \
        --sample_id ${sample_id} \
        --threshold_percentile 0.99 \
        --distance_method manhattan \
        --plot_heatmap
    """
}