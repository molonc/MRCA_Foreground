process JITTER_CORRECT {
    label 'symlink_output'
    tag "${sample_id}"
    label 'process_single'

    container params.medicc2_container

    input:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn), path(all_nnd), path(punish), path(adata), path(cell_list), path(hg19_centro)
    output:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn), path(all_nnd), path(punish), path(adata), path("output_segments_all.tsv.gz"), emit: master
        path "centromere_bounds.tsv", emit: centromere_bounds

    script:
    """
    echo 'f'
    python ${moduleDir}/compute_segments_mspcf.py \
        ${adata} \
        ${cell_list} \
        output_segments_all.tsv.gz \
        --hg19_blacklist_path  ${hg19_centro}
    """
}