
process PREPROCESS_MEDICC2 {
    label 'symlink_output'

    tag "${sample_id}"
    label 'process_single'

    container params.medicc2_container

    input:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn), path(all_nnd), path(punish)

    output:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(hscn), path(all_nnd), path(punish), path("anndata.h5ad"), path("medicc2_input_cell_list.txt"), emit: master

    script:
    """
    echo 'f'
    python ${moduleDir}/create_signals_and_cell_files.py \
        ${hscn} \
        ${metrics} \
        --ann_data anndata.h5ad \
        --cell_name_output medicc2_input_cell_list.txt
    """
}
