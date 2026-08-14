process SIGNALS_PREP {
    label 'symlink_output'
    tag "${library_id}-${sample_id}-${exp_con}"
    label 'process_single'
    conda "${moduleDir}/r44.yml"

    input:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads_input), path(alleles), path(metrics)

    output:
        tuple val(sample_id), val(library_id), val(exp_con), path("${library_id}_reads.csv.gz"), path("${library_id}_alleles.csv.gz"), path("${library_id}_metrics.csv.gz"), emit: master
    script:
    """
    Rscript ${projectDir}/modules/fg_modules/signals_prep/hmmcopy_SIGNALS_prep.R --sample_id "${sample_id}" \
        --library_id "${library_id}" \
        --reads "${reads_input}" \
        --alleles "${alleles}" \
        --metrics "${metrics}" \
        --include_s_phase ${params.include_s_phase} \
        --s_phase_quality_threshold 0.5 \
        --include_low_quality_cells FALSE \
        --quality_threshold 0.75 \
        ${params.num_cells_returned != null ? "--num_cells_returned ${params.num_cells_returned}" : ""}
    """
    }

        
