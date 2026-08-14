process NND_PUNISHMENT_CALC {
    label 'symlink_output'
    tag "${library_id}-${sample_id}-${exp_con}"
    label 'process_single'

    conda "${moduleDir}/r44.yml"

    input:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads), path(alleles), path(metrics), path(hscn)
    output:
        tuple val(sample_id), val(library_id), val(exp_con),path(reads), path(alleles), path(metrics), path(hscn), path("${library_id}_all_cells_nnd.csv"), emit: master

        path "${library_id}_all_cells_nnd.csv"
        path "${library_id}_beta_parameters.csv"
        path "${library_id}_divergent_cells_heatmap.png"
        path "${library_id}_divergent_cells_nnd.csv"
        path "${library_id}_nnd_beta_qqplot.png"
        path "${library_id}_nnd_density.png"
        path "${library_id}_single_cell_copy_divergent_cells.png"

    script:
    """
    Rscript ${projectDir}/modules/fg_modules/nnd/nnd_beta_divergence_punishment.R \
        --hscn ${hscn} \
        --sample_id ${library_id} \
        --threshold_percentile 0.99 \
        --distance_method manhattan \
        --plot_heatmap
    """
}