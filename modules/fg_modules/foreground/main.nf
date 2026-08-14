process FOREGROUND {
    label 'symlink_output'

    tag "${library_id}-${sample_id}-${exp_con}"
    label 'process_single'

    conda "${moduleDir}/medicc2.yaml"

    input:
        tuple val(sample_id), val(library_id), val(exp_con),path(reads), path(alleles), path(metrics),  path(hscn), path(tree), path(cnprofiles)
    output:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads), path(alleles), path(metrics),  path(hscn), path(tree), path(cnprofiles), path("${library_id}_seg.csv.gz"), emit: master

    script:
    """
    Rscript ${moduleDir}/foreground_ancestor.R \
        ${tree} \
        ${cnprofiles} \
        ${library_id}_seg.csv.gz
    """
}