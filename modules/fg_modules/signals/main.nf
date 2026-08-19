process SIGNALS {
    label 'symlink_output'
    tag "${sample_id}"
    label "process_medium"
    container params.base_container
    label "error_retry"

    input:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path(phasing_object)

    output:
        tuple val(sample_id), path(reads), path(alleles), path(metrics), path("${sample_id}_hscn.csv.gz"), emit: master
    script:
        if (task.attempt == 1) {
            """
            Rscript ${projectDir}/modules/fg_modules/signals/signals.R \
                --reads ${reads} \
                --allele_id ${alleles} \
                --ncores ${task.cpus} \
                --output ${sample_id}_hscn.csv.gz \
                --phased_haplotypes ${phasing_object} \
                --include_chr_x TRUE
            """
        } else {
            """
            Rscript ${projectDir}/modules/fg_modules/signals/signals.R \
                --reads ${reads} \
                --allele_id ${alleles} \
                --ncores ${task.cpus} \
                --output ${sample_id}_hscn.csv.gz \
                --phased_haplotypes ${phasing_object} \
                --include_chr_x FALSE
            """
        }
    }