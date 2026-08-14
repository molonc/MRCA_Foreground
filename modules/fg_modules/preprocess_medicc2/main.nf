
process PREPROCESS_MEDICC2 {
    label 'symlink_output'

    tag "${library_id}-${sample_id}-${exp_con}"
    label 'process_single'

    conda "${moduleDir}/medicc2.yaml"

    input:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads), path(alleles), path(metrics), path(hscn), path(all_nnd), path(punish)

    output:
        tuple val(sample_id), val(library_id),val(exp_con), path(reads), path(alleles), path(metrics), path(hscn), path(all_nnd), path(punish), path("anndata.h5ad"), path("medicc2_input_cell_list.txt"), emit: master

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
