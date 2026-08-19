process BAF_PUNISHMENT_CALC {
    label 'symlink_output'
    tag "${sample_id}"
    label 'process_single'

    container params.base_container

    input:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn), path(all_nnd)
    output:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn), path(all_nnd), path("punishment_data.csv.gz"), emit: master

        path "baf_integerness_results.csv.gz"
        path "punishment_data.csv.gz"
        path "baf_punishment_boxplot.png"
        path "epsilon_density.png"

    script:
    """
    Rscript ${projectDir}/modules/fg_modules/punish_baf/BAF_punishment_score_metric.R \
        --hscn ${hscn} \
        --metrics ${metrics} \
        --output baf_integerness.csv.gz \
        --epsilon_cutoff 0.1 \
        --bin_size 5e5
    """
}