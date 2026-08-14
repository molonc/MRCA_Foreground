process SIGNALS {
    label 'symlink_output'
    tag "${library_id}-${sample_id}-${exp_con}"
    label "process_medium"
    conda "${moduleDir}/r44.yml"
    label "error_retry"

    input:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads), path(alleles), path(metrics), path(phasing_object)

    output:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads), path(alleles), path(metrics), path("${library_id}_hscn.csv.gz"), emit: master
    script:
        if (task.attempt == 1) {
            """
            Rscript ${projectDir}/modules/fg_modules/signals/signals.R \
                --reads ${reads} \
                --allele_id ${alleles} \
                --ncores ${task.cpus} \
                --output ${library_id}_hscn.csv.gz \
                --phased_haplotypes ${phasing_object} \
                --include_chr_x TRUE
            """
        } else {
            """
            Rscript ${projectDir}/modules/fg_modules/signals/signals.R \
                --reads ${reads} \
                --allele_id ${alleles} \
                --ncores ${task.cpus} \
                --output ${library_id}_hscn.csv.gz \
                --phased_haplotypes ${phasing_object} \
                --include_chr_x FALSE
            """
        }
    }