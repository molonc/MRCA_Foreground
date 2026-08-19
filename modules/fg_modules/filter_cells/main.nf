
process FILTER_PUNISHED {
    label 'symlink_output'
    tag "${sample_id}"
    label 'process_single'

    container params.base_container
    input:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn), path(all_nnd), path(punish)

    output:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path("${sample_id}_hscn_filtered.csv.gz"), path(all_nnd), path(punish), emit: master

    script:
    """
    Rscript ${projectDir}/modules/fg_modules/filter_cells/filter_cells.R \
        --hscn ${hscn} \
        --punishment_data ${punish} \
        --divergence_data ${all_nnd} \
        --output ${sample_id}_hscn_filtered.csv.gz \
        --error_threshold 50 \
        --filter_divergent ${params.filter_divergent}
    """
}