
process FILTER_PUNISHED {
    label 'symlink_output'
    tag "${library_id}-${sample_id}-${exp_con}"
    label 'process_single'

    conda "${moduleDir}/r44.yml"
    input:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads), path(alleles), path(metrics), path(hscn), path(all_nnd), path(punish)

    output:
        tuple val(sample_id), val(library_id), val(exp_con),  path(reads), path(alleles), path(metrics), path("${library_id}_hscn_filtered.csv.gz"), path(all_nnd), path(punish), emit: master

    script:
    """
    Rscript ${projectDir}/modules/fg_modules/filter_cells/filter_cells.R \
        --hscn ${hscn} \
        --punishment_data ${punish} \
        --divergence_data ${all_nnd} \
        --output ${library_id}_hscn_filtered.csv.gz \
        --error_threshold 50 \
        --filter_divergent ${params.filter_divergent}
    """
}