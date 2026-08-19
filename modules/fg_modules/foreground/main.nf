process FOREGROUND {
    label 'symlink_output'

    tag "${sample_id}"
    label 'process_single'

    // foreground_ancestor.R needs vroom/ape/dplyr/tidyr/tibble (base_container's
    // R env) -- not medicc2_container, which is Python-only and has no Rscript.
    container params.base_container

    input:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles)
    output:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path("${sample_id}_seg.csv.gz"), emit: master

    script:
    """
    Rscript ${moduleDir}/foreground_ancestor.R \
        ${tree} \
        ${cnprofiles} \
        ${sample_id}_seg.csv.gz
    """
}