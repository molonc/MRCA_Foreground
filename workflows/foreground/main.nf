/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SIGNALS_PREP                 } from '../../modules/fg_modules/signals_prep/main.nf'
include { SIGNALS                 } from '../../modules/fg_modules/signals/main.nf'
include { NND_PUNISHMENT_CALC                 } from '../../modules/fg_modules/nnd/main.nf'
include { BAF_PUNISHMENT_CALC                 } from '../../modules/fg_modules/punish_baf/main.nf'
include { FILTER_PUNISHED                  } from '../../modules/fg_modules/filter_cells/main.nf'
include { PREPROCESS_MEDICC2                 } from '../../modules/fg_modules/preprocess_medicc2/main.nf'
include { MEDICC2                 } from '../../modules/fg_modules/medicc2/main.nf'
include { JITTER_CORRECT                } from '../../modules/fg_modules/jitter_correct/main.nf'
include { FOREGROUND                 } from '../../modules/fg_modules/foreground/main.nf'
include { GET_DF                 } from '../../modules/fg_modules/get_df/main.nf'
include { ANNOTATE_BPS                  } from '../../modules/fg_modules/annotate_bps/main.nf'
include { GET_BP_FILE                  } from '../../modules/fg_modules/get_bp_file/main.nf'
include { GET_HDP_FILE                  } from '../../modules/fg_modules/get_hdp_file/main.nf'
include { PLOT_HEATMAPS                  } from '../../modules/fg_modules/plot_heatmaps/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FOREGROUND_PIPE {
    take:
        ch_fk_sheet
    main:

        SIGNALS_PREP(
            ch_fk_sheet
        )
        phasing_object_ch = params.phasing_object ?
            Channel.fromPath(params.phasing_object)
            : Channel.value( file("${workflow.workDir}/dummy_phasing.csv").with { f ->
                f.text = ""; f
            })
        signals_in = SIGNALS_PREP.out.master.combine(phasing_object_ch)
        SIGNALS(
            signals_in
        )

        NND_PUNISHMENT_CALC(
            SIGNALS.out.master
        )

        BAF_PUNISHMENT_CALC(
            NND_PUNISHMENT_CALC.out.master
        )

        FILTER_PUNISHED(
            BAF_PUNISHMENT_CALC.out.master
        )

        PREPROCESS_MEDICC2(
            FILTER_PUNISHED.out.master
        )

        def centro_loc_ch = Channel.fromPath("${projectDir}/helper/hg19.chr2centro.json")
        JITTER_CORRECT(
            PREPROCESS_MEDICC2.out.master.combine(centro_loc_ch)
        )

        def tcn_bool = Channel.value(params.tcn) // True if medicc2 only does tcn output
        MEDICC2(
            JITTER_CORRECT.out.master,
            tcn_bool
        )
        FOREGROUND(
            MEDICC2.out.master
        )

        // Create reads data frame with new "foreground state"
        def dummy_cell = Channel.fromPath("${projectDir}/helper/dummy_cell.csv")
        def cen_info_ch = Channel.fromPath("${projectDir}/helper/centro_telo_locs.csv")
        def breakpoint_functions = Channel.fromPath("${projectDir}/helper/breakpoint_file_functions.R")
        GET_DF(
            FOREGROUND.out.master.combine(dummy_cell.combine(cen_info_ch.combine(breakpoint_functions))),
            tcn_bool
        )
        
        ANNOTATE_BPS(
            GET_DF.out.master,
            tcn_bool
        )
        GET_HDP_FILE(
            ANNOTATE_BPS.out.master
        )
        GET_BP_FILE(
            ANNOTATE_BPS.out.master
        )
        def annotated_bins  = Channel.fromPath("${projectDir}/helper/hg19_annotated_bins.csv")
        def heatmap_ch = GET_HDP_FILE.out.master.join(GET_BP_FILE.out.master, by: 0).combine(annotated_bins)
        PLOT_HEATMAPS(
            heatmap_ch
        )
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/