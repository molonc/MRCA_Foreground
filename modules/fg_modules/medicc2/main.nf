process MEDICC2 {
    label 'symlink_output'
    tag "${library_id}-${sample_id}-${exp_con}"
    label 'process_high'

    conda "${moduleDir}/medicc2.yaml"

    input:
        tuple val(sample_id), val(library_id), val(exp_con),path(reads), path(alleles), path(metrics), path(hscn), path(all_nnd), path(punish), path(adata), path(jitter_segs)
        val tcn_bool
    output:
        tuple val(sample_id), val(library_id),  val(exp_con),path(reads), path(alleles), path(metrics), path(hscn), path("medicc2_output_allele_wgd/output_segments_all.tsv_final_tree.new"),  path("medicc2_output_allele_wgd/output_segments_all.tsv_final_cn_profiles.tsv"), emit: master

    script:
    """
    if [[ "${tcn_bool}" == "true" ]]; then
        medicc2 ${jitter_segs} medicc2_output_allele_wgd \
            --n-cores "${task.cpus}" \
            --total-copy-numbers \
            --input-allele-columns state
    else
        medicc2 ${jitter_segs} medicc2_output_allele_wgd \
            --n-cores "${task.cpus}" \
            --events \
            --wgd-x2
    fi
    """
}

//    --events ## allele
//    --input-allele-columns state  ## total cn
//     --total-copy-numbers \



// No allele awareess, no WGD awareness
    // medicc2 ${jitter_segs} medicc2_output_allele_wgd \
    //     --n-cores "${task.cpus}" \
    //     --total-copy-numbers \
    //     --input-allele-columns state \
    //     --no-wgd

// No allele awareess, yes WGD awareness
    // medicc2 ${jitter_segs} medicc2_output_allele_wgd \
    //     --n-cores "${task.cpus}" \
    //     --total-copy-numbers \
    //     --input-allele-columns state \

// No allele awareess, yes WGD awareness 2x. Whatever 2x is. It doesnt work in TCN
    // medicc2 ${jitter_segs} medicc2_output_allele_wgd \
    //     --n-cores "${task.cpus}" \
    //     --total-copy-numbers \
    //     --input-allele-columns state \
    //     --wgd-x2

// yes Allele awareess, no WGD awarenes
    // medicc2 ${jitter_segs} medicc2_output_allele_wgd \
    //     --n-cores "${task.cpus}" \
    //     --no-wgd \
    //     --events

// yes Allele awareess, yes WGD awareness
    // medicc2 ${jitter_segs} medicc2_output_allele_wgd \
    //     --n-cores "${task.cpus}" \
    //     --events

// yes Allele awareess, yes WGD awareness, 2x wgd
    // medicc2 ${jitter_segs} medicc2_output_allele_wgd \
    //     --n-cores "${task.cpus}" \
    //     --events \
    //     --wgd-x2
