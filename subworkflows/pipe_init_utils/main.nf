//
// Subworkflows for pipeline initialisation
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
        input_path
        data_dir

    main:
        sample_id_ch = Channel
            .fromPath(input_path)
            .splitCsv(sep: '\t', header: false)
            .map { row -> tuple(row[0], row[1], row[2], row[3], row[4])} //library_id, sample_id, experimental_con, qc_ax, ocunthaps_ax
        // sample_id_ch.view()

        def assembly_paths = sample_id_ch.map { library_id, sample_id, exp_con, qc, counthaps ->
            def dir = "${data_dir}/${library_id}_${sample_id}_${exp_con}"

            def reads   = file("${dir}/*hmmcopy_reads.csv.gz", checkIfExists: true)
            def allele  = file("${dir}/*allele_count.csv.gz", checkIfExists: true)
            def metrics = file("${dir}/*hmmcopy_metrics.csv.gz", checkIfExists: true)

            tuple(sample_id, library_id, exp_con, reads, allele, metrics)
        }
        // assembly_paths.view()
    emit:
        ch_fk_sheet = assembly_paths
}
