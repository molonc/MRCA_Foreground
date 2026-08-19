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

    The samplesheet is the single source of truth for where a library's input files
    live. There is no directory-naming convention or filename globbing here on
    purpose: different data sources (a shared data lake, a collaborator who just
    handed you three files, an old cohort re-downloaded into some other layout)
    almost never agree on a folder schema, and guessing at one is exactly what made
    this fragile before -- a folder that happened to contain two files matching
    `*allele_count.csv.gz` silently produced a broken run instead of a clear error.

    Required columns (header row, comma-separated):
        sample_id, reads, allele, metrics

    `sample_id` is the pipeline's only tracking key -- it names the results
    subdirectory (see conf/process.config) and is the tag on every process, from
    here through the rest of the pipeline. What goes into it is entirely up to you:
    it can just be a library ID, or something composite like
    `<library_id>_<sample_id>_<exp_con>` if that's what uniquely identifies a run in
    your data -- the pipeline itself has no opinion on its structure, only that it's
    non-empty and unique across the sheet.

    `reads` / `allele` / `metrics` must each be an exact path to one file -- absolute,
    or relative to the samplesheet's own location (so a samplesheet + its data can be
    kept/shared together without hardcoding a machine-specific prefix).

    Start a new sheet from `input_files/samplesheets/TEMPLATE.csv`, or use
    `helper_scripts/generate_samplesheet.py` to build/migrate one from a directory that
    already follows the old `<library_id>_<sample_id>_<exp_con>/` convention (it composes
    those three fields into one sample_id for you, matching the previous behaviour).
*/

def REQUIRED_COLUMNS = ['sample_id', 'reads', 'allele', 'metrics']

// Resolve one samplesheet cell to an existing file, erroring out with row/column
// context instead of a bare Nextflow/Groovy stack trace if it's blank or missing.
def resolveSheetPath(row, column, sheet_dir, row_desc) {
    def raw = row[column]?.trim()
    if (!raw) {
        error "Samplesheet row for ${row_desc} is missing a value in the required '${column}' column."
    }
    def resolved = raw.startsWith('/') ? file(raw, glob: false) : file("${sheet_dir}/${raw}", glob: false)
    if (!resolved.exists()) {
        error "Samplesheet row for ${row_desc}: '${column}' file does not exist: ${resolved}\n(from sheet value: '${raw}')"
    }
    resolved
}

workflow PIPELINE_INITIALISATION {

    take:
        input_path

    main:
        if (!input_path) {
            error "params.input is required: point it at a samplesheet with a header row of\n" +
                  "${REQUIRED_COLUMNS.join(',')}\n" +
                  "Copy input_files/samplesheets/TEMPLATE.csv to start one, or run\n" +
                  "helper_scripts/generate_samplesheet.py to build/migrate one from an existing\n" +
                  "<library_id>_<sample_id>_<exp_con>/ directory layout."
        }
        def sheet_file = file(input_path, checkIfExists: true)
        def sheet_dir  = sheet_file.getParent()

        // Collected into one list (rather than streamed row-by-row) so duplicate
        // sample_ids can be caught here, up front -- since sample_id now names the
        // results directory, two rows sharing one would silently clobber each
        // other's output instead of erroring.
        def assembly_paths = Channel
            .fromPath(sheet_file)
            .splitCsv(header: true)
            .toList()
            .flatMap { rows ->
                def seen = new HashSet()
                rows.collect { row ->
                    def missing_cols = REQUIRED_COLUMNS.findAll { !row.containsKey(it) }
                    if (missing_cols) {
                        error "Samplesheet ${sheet_file} is missing required column(s): ${missing_cols.join(', ')}.\n" +
                              "Expected header: ${REQUIRED_COLUMNS.join(',')}"
                    }

                    def sample_id = row.sample_id?.trim()
                    if (!sample_id) {
                        error "Samplesheet ${sheet_file} has a row with a blank sample_id: ${row}"
                    }
                    if (!seen.add(sample_id)) {
                        error "Samplesheet ${sheet_file} has more than one row with sample_id='${sample_id}'.\n" +
                              "sample_id must be unique -- it names the results directory and is this pipeline's\n" +
                              "only tracking key."
                    }
                    def row_desc = "sample_id='${sample_id}'"

                    def reads   = resolveSheetPath(row, 'reads',   sheet_dir, row_desc)
                    def allele  = resolveSheetPath(row, 'allele',  sheet_dir, row_desc)
                    def metrics = resolveSheetPath(row, 'metrics', sheet_dir, row_desc)

                    tuple(sample_id, reads, allele, metrics)
                }
            }
        // assembly_paths.view()
    emit:
        ch_fk_sheet = assembly_paths
}
