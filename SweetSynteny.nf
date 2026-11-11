#!/usr/bin/env nextflow

/*
 * SweetSynteny
 * Authors:
 * - Maria Schreiber <maria.schreiber@uni-jena.de>
 *
 * Revised workflow:
 * - channel-safe empty-input checks
 * - path-staged scientific inputs and helper scripts
 * - explicit search/neighbourhood output contracts
 * - validated clustering mode
 * - retry resources scale on retryable resource failures
 * - record run parameters in 0_metadata/run_provenance.json
 * - split hmmscan/cmscan representative FASTAs into 25-sequence chunks
 *   and scan chunks concurrently within each Nextflow task
 * - forward all clustering options of plot_context_script.py (feature_mode,
 *   min_coverage, min_cluster_cosine, HDBSCAN options, ...); unset = script default
 * - HDBSCAN ('density') on by default; density_clustering = yes | no | only
 */

import groovy.json.JsonOutput

nextflow.enable.dsl=2

// Optional clustering parameters. null = not passed on, i.e. the default of
// bin/plot_context_script.py applies. Set any of them in the params file / on the CLI.
[
    'feature_mode', 'position_weight', 'synteny_weight', 'adjacency_across_gaps',
    'collapse_overlaps', 'overlap_fraction', 'rare_unit', 'min_features_per_neighbourhood',
    'min_cluster_size', 'min_coverage', 'min_cluster_cosine', 'cohesion_stat',
    'min_shared_families', 'signature_fraction', 'assign_leftovers', 'assign_min_similarity',
    'supercluster_similarity', 'keep_small_hier_clusters', 'max_hierarchical_samples',
    'max_heatmap_features', 'density_clustering', 'hdbscan_selection', 'min_samples',
    'noise_penalty', 'max_hdbscan_candidates', 'max_silhouette_samples', 'logo_for_noise'
].each { name -> if (!params.containsKey(name)) params[name] = null }

version = '1.2.0'

// Quote arbitrary values for safe insertion into shell command arguments.
def shellQuote = { Object value ->
    if (value == null) return "''"
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

// "--flag 'value'" when the value is set, otherwise an empty string.
def optArg = { String flag, Object value ->
    value != null ? "${flag} ${shellQuote(value)}" : ''
}

def isTrue = { value ->
    value == true || value?.toString()?.equalsIgnoreCase('true')
}

def asList = { value ->
    if (value == null) return []
    if (value instanceof Collection) return value.toList()
    return value.toString().split(/[\s,]+/).findAll { it != null && !it.toString().trim().isEmpty() }
}

def hasNonCommentRow = { resultFile ->
    boolean found = false
    resultFile.withReader { reader ->
        String line
        while (!found && (line = reader.readLine()) != null) {
            def trimmed = line.trim()
            if (trimmed && !trimmed.startsWith('#')) found = true
        }
    }
    return found
}

log.info """
SweetSynteny - Nextflow Pipeline
================================
Pipeline version : ${version}
Nextflow version : ${workflow.nextflow.version}
Run name         : ${workflow.runName}
Session ID       : ${workflow.sessionId}
Revision         : ${workflow.revision ?: 'unknown'}
Launch directory : ${workflow.launchDir}
Work directory   : ${workflow.workDir}

> Results
  Result folder  : ${params.output_dir}
> Searching
  Search type    : ${params.search_types} [blastn, blastp, infernal, tblastn, from_gff]
  Bio type       : ${params.bio_type}
  Genomes folder : ${params.genomes_dir}
  Query          : ${params.query}
  Annotation     : ${params.annotation_type}
> Neighbourhoods
  Gene of interest : ${params.gene_of_interest}
  Neighbours       : ${params.neighbours}
  Included features: ${params.including_features}
  Ignore overlaps  : ${params.ignore_overlaps}
  Substring search : ${params.substring_search}
> Clustering
  Protein/ncRNA    : ${params.adjacent_gene_clustering}
  Scale            : ${params.scale}
  Feature mode     : ${params.feature_mode ?: 'default'}
  Metric / linkage : ${params.metric ?: 'default'} / ${params.linkage_method ?: 'default'}
  H-cluster cut    : ${params.cut_height_args ?: 'automatic'}
  Min coverage     : ${params.min_coverage ?: 'default'}
  Min cluster cos. : ${params.min_cluster_cosine ?: 'default'}
  HDBSCAN          : ${params.density_clustering ?: 'yes (default)'}
  Microsynteny logo: ${params.microsynteny_logo}
> Resources
  Search CPUs      : ${params.cpus}
""".stripIndent(true)

// -----------------------------------------------------------------------------
// Search
// -----------------------------------------------------------------------------
process runSearch {
    publishDir "${params.output_dir}/1_search", mode: 'copy', pattern: '*.tsv'

    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 2
    cpus { params.cpus as int }
    memory { 4.GB * task.attempt }
    time { 1.h * task.attempt }

    input:
    tuple val(id), path(genome), path(gff)
    path query_file, name: 'query_input'

    output:
    tuple val(id), path("${id}.tsv"), emit: search_results

    script:
    def searchType = params.search_types.toString().trim().toLowerCase()
    def outFile = "${id}.tsv"
    def blastOutfmt = '6 qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq'
    def evalueArg = params.evalue_threshold_blast != null ? "-evalue ${shellQuote(params.evalue_threshold_blast)}" : ''
    def matrixArg = params.blast_matrix != null ? "-matrix ${shellQuote(params.blast_matrix)}" : ''

    if (searchType == 'blastn') {
        """
        makeblastdb \\
            -in ${shellQuote(genome)} \\
            -dbtype nucl \\
            -out blastdb

        blastn \\
            -num_threads ${task.cpus} \\
            -query ${shellQuote(query_file)} \\
            -subject ${shellQuote(genome)} \\
            -out ${shellQuote(outFile)} \\
            -outfmt ${shellQuote(blastOutfmt)} \\
            ${evalueArg}
        """
    }
    else if (searchType == 'blastp') {
        """
        makeblastdb \\
            -in ${shellQuote(genome)} \\
            -dbtype prot \\
            -out blastdb

        blastp \\
            -num_threads ${task.cpus} \\
            -query ${shellQuote(query_file)} \\
            -db blastdb \\
            -out ${shellQuote(outFile)} \\
            -outfmt ${shellQuote(blastOutfmt)} \\
            ${evalueArg} \\
            ${matrixArg}
        """
    }
    else if (searchType == 'tblastn') {
        """
        makeblastdb \\
            -in ${shellQuote(genome)} \\
            -dbtype nucl \\
            -out blastdb

        tblastn \\
            -num_threads ${task.cpus} \\
            -query ${shellQuote(query_file)} \\
            -db blastdb \\
            -out ${shellQuote(outFile)} \\
            -outfmt ${shellQuote(blastOutfmt)} \\
            ${evalueArg} \\
            ${matrixArg}
        """
    }
    else if (searchType == 'infernal') {
        """
        cmsearch \\
            --cpu ${task.cpus} \\
            --tblout ${shellQuote(outFile)} \\
            ${shellQuote(query_file)} \\
            ${shellQuote(genome)}
        """
    }
    else {
        error "runSearch received unsupported search type: ${searchType}"
    }
}

// -----------------------------------------------------------------------------
// Neighbourhood extraction
// -----------------------------------------------------------------------------
process getNeighbours {
    publishDir "${params.output_dir}/2_neighbour", mode: 'copy', pattern: '*.nb.*'

    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 2
    cpus 1
    memory { 2.GB * task.attempt }
    time { 30.min * task.attempt }

    input:
    tuple val(id), path(genome), path(gff), path(search_result, name: 'search_result.tsv'), val(use_hit_file)
    path get_neighbours_script

    output:
    tuple val(id), path("${id}.nb.tsv"), emit: tsv
    path("${id}.nb.protein.mfaa"), optional: true, emit: mfaa
    path("${id}.nb.ncrna.mfna"), optional: true, emit: mfna
    path("${id}.nb.promoter.mfna"), optional: true, emit: promoter_mfna

    script:
    def features = asList(params.including_features)
    if (!features) error 'Parameter --including_features must contain at least one feature.'
    def featureArgs = features.collect { shellQuote(it) }.join(' ')

    def hitArg = isTrue(use_hit_file) ? "--hit_file ${shellQuote(search_result)}" : ''
    def subSearchArg = isTrue(params.substring_search) ? '--substring_search' : ''
    def ignoreOverArg = isTrue(params.ignore_overlaps) ? '--ignore_overlaps' : ''
    def promoterModeArg = params.promoter_mode != null ? "--promoter_mode ${shellQuote(params.promoter_mode)}" : ''
    def overlapArg = params.overlap_threshold != null ? "--overlap_threshold ${shellQuote(params.overlap_threshold)}" : ''
    def promoterArg = params.promoter != null ? "--promoter ${shellQuote(params.promoter)}" : ''
    def promoterLenArg = params.promoter_len != null ? "--promoter_len ${shellQuote(params.promoter_len)}" : ''
    def blastEvalueArg = params.evalue_threshold_blast != null ? "--evalue_threshold_blast ${shellQuote(params.evalue_threshold_blast)}" : ''
    def infernalEvalueArg = params.evalue_threshold_infernal != null ? "--evalue_threshold_infernal ${shellQuote(params.evalue_threshold_infernal)}" : ''
    def blastLenArg = params.len_threshold_blast != null ? "--len_threshold_blast ${shellQuote(params.len_threshold_blast)}" : ''
    def infernalLenArg = params.len_threshold_infernal != null ? "--len_threshold_infernal ${shellQuote(params.len_threshold_infernal)}" : ''
    def genomeIdArg = "--genome_id ${shellQuote(id)}"

    """
    python ${shellQuote(get_neighbours_script)} \\
        --bio_type ${shellQuote(params.bio_type)} \\
        --input_type ${shellQuote(params.search_types)} \\
        --fna_file ${shellQuote(genome)} \\
        --gff_file ${shellQuote(gff)} \\
        --including_features ${featureArgs} \\
        --gene_of_interest ${shellQuote(params.gene_of_interest)} \\
        --neighbours ${shellQuote(params.neighbours)} \\
        --output_path ${shellQuote("${id}.nb")} \\
        ${genomeIdArg} \\
        ${promoterArg} \\
        ${promoterLenArg} \\
        ${promoterModeArg} \\
        ${subSearchArg} \\
        ${ignoreOverArg} \\
        ${blastEvalueArg} \\
        ${infernalEvalueArg} \\
        ${blastLenArg} \\
        ${infernalLenArg} \\
        ${overlapArg} \\
        ${hitArg}
    """
}

// -----------------------------------------------------------------------------
// MMseqs clustering
// -----------------------------------------------------------------------------
process runMMseqsProtein {
    publishDir "${params.output_dir}/3_cluster", mode: 'copy'

    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 2
    cpus 4
    memory { 8.GB * task.attempt }
    time { 2.h * task.attempt }

    input:
    path protein_mfaa
    path postprocess_script

    output:
    tuple path('clusterRes.protein_mmseqs.tsv'), path('prot_clust_rep_seq.fasta'), emit: mmseqs_out

    script:
    def fastaFiles = protein_mfaa instanceof Collection ? protein_mfaa : [protein_mfaa]
    def fastaArgs = fastaFiles.collect { shellQuote(it) }.join(' ')

    """
    cat -- ${fastaArgs} > all_protein.mfaa

    # Exclude GOI by exact FASTA identifier (first header token); supports wrapped FASTA.
    awk -v target=${shellQuote(params.gene_of_interest)} \\
        'BEGIN { keep = 1 } /^>/ { id = substr(\$0, 2); sub(/[[:space:]].*/, "", id); keep = (id != target) } keep { print }' \\
        all_protein.mfaa > merged_protein.mfaa

    if ! grep -q '^>' merged_protein.mfaa; then
        echo 'No protein sequences remain after excluding the gene of interest.' >&2
        exit 2
    fi

    mmseqs easy-linclust merged_protein.mfaa prot_clust tmp \\
        --min-seq-id 0.3 \\
        --cov-mode 1 \\
        --threads ${task.cpus}

    python ${shellQuote(postprocess_script)} \\
        --input prot_clust_cluster.tsv \\
        --output clusterRes.protein_mmseqs.tsv \\
        --gene_of_interest ${shellQuote(params.gene_of_interest)}
    """
}

process runMMseqsNCRNA {
    publishDir "${params.output_dir}/3_cluster", mode: 'copy'

    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 2
    cpus 4
    memory { 8.GB * task.attempt }
    time { 2.h * task.attempt }

    input:
    path ncrna_mfna
    path postprocess_script

    output:
    tuple path('clusterRes.ncrna_mmseqs.tsv'),
          path('ncrna_clust_rep_seq.fasta'),
          emit: mmseqs_out

    script:
    def fastaFiles = ncrna_mfna instanceof Collection ? ncrna_mfna : [ncrna_mfna]
    def fastaArgs = fastaFiles.collect { shellQuote(it) }.join(' ')

    """
    cat -- ${fastaArgs} > all_ncrna.mfna

    # Exclude sis features before clustering.
    awk '
        /^>/ {
            header = substr(\$0, 2)
            id = header
            sub(/[[:space:]].*/, "", id)

            keep = 1

            # Exclude sis identifiers.
            if (id == "sis" || index(id, ":sis_sis_") > 0) {
                keep = 0
            }
        }
        keep { print }
    ' all_ncrna.mfna > merged_ncrna.mfna

    if ! grep -q '^>' merged_ncrna.mfna; then
        echo 'No ncRNA sequences remain after filtering sis features.' >&2
        exit 2
    fi

    mmseqs easy-linclust merged_ncrna.mfna ncrna_clust tmp \\
        --min-seq-id 0.8 \\
        --cov-mode 1 \\
        --threads ${task.cpus}

    python ${shellQuote(postprocess_script)} \\
        --input ncrna_clust_cluster.tsv \\
        --output clusterRes.ncrna_mmseqs.tsv \\
        --gene_of_interest ${shellQuote(params.gene_of_interest)}
    """
}

// -----------------------------------------------------------------------------
// Domain/family annotation of MMseqs representatives
// -----------------------------------------------------------------------------
process runHMMscan {
    publishDir "${params.output_dir}/3_cluster", mode: 'copy'

    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 2
    cpus 10
    memory { 32.GB * task.attempt }
    time { 7.h * task.attempt }

    input:
    tuple path(rep_fasta), path(mmseqs_map)
    path pfam_hmm
    path pfam_h3f
    path pfam_h3i
    path pfam_h3m
    path pfam_h3p
    path postprocess_script

    output:
    path 'clusterRes.hmm.protein_cluster.tsv', emit: hmm_results

    script:
    """
    set -euo pipefail

    # Confirm that the index files have been staged into this task.
    for f in Pfam-A.hmm.h3f Pfam-A.hmm.h3i Pfam-A.hmm.h3m Pfam-A.hmm.h3p; do
        if [ ! -s "\$f" ]; then
            echo "Missing or empty Pfam index file: \$f" >&2
            exit 1
        fi
    done

    mkdir -p hmm_chunks hmm_tblout

    # Split FASTA into files containing at most 25 records each.
    awk '
    /^>/ {
        seq++
        if ((seq - 1) % 25 == 0) {
            if (file != "") close(file)
            file = sprintf("hmm_chunks/chunk_%05d.fasta", int((seq - 1) / 25) + 1)
        }
    }
    seq > 0 { print > file }
    END {
        if (seq == 0) exit 2
        if (file != "") close(file)
    }
    ' ${shellQuote(rep_fasta)}

    # Run one hmmscan per chunk concurrently; xargs bounds parallelism to task.cpus.
    export PFAM_DB=${shellQuote(pfam_hmm)}
    find hmm_chunks -maxdepth 1 -type f -name 'chunk_*.fasta' -print0 |
        xargs -0 -r -n 1 -P ${task.cpus} bash -c '
            chunk="\$1"
            base=\$(basename "\$chunk" .fasta)
            hmmscan --cpu 1 \\
                --tblout "hmm_tblout/\${base}.tbl" \\
                "\$PFAM_DB" "\$chunk" \\
                > "hmm_chunks/\${base}.log" 2>&1
        ' _

    # Preserve one comment/header section and append only hit rows from other tables.
    first=1
    for tbl in hmm_tblout/chunk_*.tbl; do
        if [ "\$first" -eq 1 ]; then
            cat "\$tbl" > merged.protein.hmm.tbl
            first=0
        else
            awk '!/^#/' "\$tbl" >> merged.protein.hmm.tbl
        fi
    done

    if [ "\$first" -eq 1 ]; then
        echo 'No hmmscan output tables were generated.' >&2
        exit 1
    fi

    python ${shellQuote(postprocess_script)} \\
        --hmmscan_file merged.protein.hmm.tbl \\
        --mmseqs_map ${shellQuote(mmseqs_map)} \\
        --output_file clusterRes.hmm.protein_cluster.tsv \\
        --gene_of_interest ${shellQuote(params.gene_of_interest)}
    """
}

process runCMscan {
    publishDir "${params.output_dir}/3_cluster", mode: 'copy'

    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 2
    cpus 4
    memory { 32.GB * task.attempt }
    time { 7.h * task.attempt }

    input:
    tuple path(rep_fasta), path(mmseqs_map)
    path rfam_db
    path rfam_i1f
    path rfam_i1i
    path rfam_i1m
    path rfam_i1p
    path postprocess_script

    output:
    path 'clusterRes.ncrna_cluster.tsv', emit: cm_results

    script:
    """
    set -euo pipefail

    # Confirm that the Infernal index files have been staged into this task.
    for f in Rfam.cm.i1f Rfam.cm.i1i Rfam.cm.i1m Rfam.cm.i1p; do
        if [ ! -s "\$f" ]; then
            echo "Missing or empty Rfam index file: \$f" >&2
            exit 1
        fi
    done

    mkdir -p cm_chunks cm_tblout

    # Split FASTA into files containing at most 25 records each.
    awk '
    /^>/ {
        seq++
        if ((seq - 1) % 25 == 0) {
            if (file != "") close(file)
            file = sprintf("cm_chunks/chunk_%05d.fasta", int((seq - 1) / 25) + 1)
        }
    }
    seq > 0 { print > file }
    END {
        if (seq == 0) exit 2
        if (file != "") close(file)
    }
    ' ${shellQuote(rep_fasta)}

    # Run one cmscan per chunk concurrently; xargs bounds parallelism to task.cpus.
    export RFAM_DB=${shellQuote(rfam_db)}
    find cm_chunks -maxdepth 1 -type f -name 'chunk_*.fasta' -print0 |
        xargs -0 -r -n 1 -P ${task.cpus} bash -c '
            chunk="\$1"
            base=\$(basename "\$chunk" .fasta)
            cmscan -E 0.01 \\
                --cpu 1 \\
                --noali \\
                --tblout "cm_tblout/\${base}.tbl" \\
                "\$RFAM_DB" "\$chunk" \\
                > "cm_chunks/\${base}.log" 2>&1
        ' _

    # Preserve one comment/header section and append only hit rows from other tables.
    first=1
    for tbl in cm_tblout/chunk_*.tbl; do
        if [ "\$first" -eq 1 ]; then
            cat "\$tbl" > clusterRes.ncrna
            first=0
        else
            awk '!/^#/' "\$tbl" >> clusterRes.ncrna
        fi
    done

    if [ "\$first" -eq 1 ]; then
        echo 'No cmscan output tables were generated.' >&2
        exit 1
    fi

    python ${shellQuote(postprocess_script)} \\
        --cmscan_file clusterRes.ncrna \\
        --mmseqs_map ${shellQuote(mmseqs_map)} \\
        --output_file clusterRes.ncrna_cluster.tsv \\
        --gene_of_interest ${shellQuote(params.gene_of_interest)}
    """
}

// -----------------------------------------------------------------------------
// Merge neighbourhood/cluster files and assign colours
// -----------------------------------------------------------------------------
process clusterColoring {
    publishDir "${params.output_dir}/3_cluster",
        mode: 'copy',
        overwrite: true

    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 2
    cpus 4
    memory { 8.GB * task.attempt }
    time { 7.h * task.attempt }

    input:
    path cluster_files
    path tsv_files
    path color_script
    path plot_script
    path color_scheme_file

    output:
    path 'merged_with_color.tsv', emit: colored_tsv
    path 'summary.txt', emit: summary_txt
    path 'summary.json', emit: summary_json
    path 'summary.svg', emit: summary_svg
    path 'summary.png', emit: summary_png

    script:
    def clusterList = cluster_files instanceof Collection ? cluster_files : [cluster_files]
    def neighbourList = tsv_files instanceof Collection ? tsv_files : [tsv_files]
    def clusterArgs = clusterList.collect { shellQuote(it) }.join(' ')
    def neighbourArgs = neighbourList.collect { shellQuote(it) }.join(' ')

    """
    # Merge only the explicitly declared inputs.
    python - ${clusterArgs} -- ${neighbourArgs} <<'PY'
    import sys

    args = sys.argv[1:]
    try:
        separator = args.index('--')
    except ValueError:
        raise SystemExit('Internal error: missing file-list separator')

    cluster_files = args[:separator]
    neighbour_files = args[separator + 1:]
    if not cluster_files:
        raise SystemExit('No cluster result files were supplied')
    if not neighbour_files:
        raise SystemExit('No neighbourhood TSV files were supplied')

    def merge_files(paths, destination, drop_promoters=False):
        with open(destination, 'w', encoding='utf-8') as out_handle:
            for filename in paths:
                with open(filename, 'r', encoding='utf-8', errors='replace') as in_handle:
                    for line in in_handle:
                        fields = line.rstrip('\\r\\n').split('\\t')
                        if drop_promoters and any(field.strip().lower() == 'promoter' for field in fields):
                            continue
                        if line.endswith(('\\n', '\\r')):
                            out_handle.write(line)
                        else:
                            out_handle.write(line + '\\n')

    merge_files(cluster_files, 'clusterRes.tsv')
    merge_files(neighbour_files, 'merged_neighbours.tsv', drop_promoters=True)
    PY

    python ${shellQuote(color_script)} \\
        --cluster_file clusterRes.tsv \\
        --tsv_file merged_neighbours.tsv \\
        --output_file merged_with_color.tsv \\
        --gene_of_interest ${shellQuote(params.gene_of_interest)} \\
        --goi_type ${shellQuote(params.search_types)} \\
        --color_scheme_file '${color_scheme_file}'
    """
}

// -----------------------------------------------------------------------------
// Plot genomic context
// -----------------------------------------------------------------------------
process plottingContext {
    publishDir "${params.output_dir}/4_plot",
        mode: 'copy',
        overwrite: true

    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 2
    cpus 4
    memory { 8.GB * task.attempt }
    time { 7.h * task.attempt }

    input:
    tuple path(merged_with_color), path(name_file, name: 'name_file.txt'), val(has_name_file)
    path plot_script
    path color_scheme_file

    output:
    path 'clustering_report.txt', optional: true
    path 'clustering_parameters.json', optional: true
    path 'hierarchical/hierarchical_cut_scan.csv', optional: true      
    path 'clustering_comparison.tsv', optional: true      
    path 'neighbourhood_assignments.tsv', optional: true
    path 'excluded_neighbourhoods.tsv', optional: true
    path 'comand.txt'
    path 'summary_*.json', optional: true, emit: json_reports
    path 'density', emit: density_dir
    path 'hierarchical', emit: hierarchical_dir
    path 'density/**/*.svg', optional: true, emit: density_svgs
    path 'hierarchical/**/*.svg', optional: true, emit: hierarchical_svgs
    path 'network_graph*', optional: true
        
    script:
    def nameFileArg = isTrue(has_name_file) ? "--name_file ${shellQuote(name_file)}" : ''
    def geneNameArg = params.gene_name != null ? "--gene_name ${shellQuote(params.gene_name)}" : ''
    def scaleArg = params.scale != null ? "--scale ${shellQuote(params.scale)}" : ''
    def cutHeightArg = params.cut_height_args != null ? "--cut_height_args ${shellQuote(params.cut_height_args)}" : ''
    def logoArg = params.microsynteny_logo != null ? "--microsynteny_logo ${shellQuote(params.microsynteny_logo)}" : ''
    def svgArg = params.svg != null ? "--svg ${shellQuote(params.svg)}" : ''
    def removeCoreArg = params.remove_core_genes != null ? "--remove_core_genes ${shellQuote(params.remove_core_genes)}" : ''
    def removeRareArg = params.remove_rare_genes != null ? "--remove_rare_genes ${shellQuote(params.remove_rare_genes)}" : ''
    def coreThresholdArg = params.core_threshold != null ? "--core_threshold ${shellQuote(params.core_threshold)}" : ''
    def rareMinCountArg = params.rare_min_count != null ? "--rare_min_count ${shellQuote(params.rare_min_count)}" : ''
    def minGenePrevalence = params.min_gene_prevalence != null ? "--min_gene_prevalence ${shellQuote(params.min_gene_prevalence)}" : ''
    def paraRe = params.paralog_resolution != null ? "--paralog_resolution ${shellQuote(params.paralog_resolution)}" : ''
    def paraSim = params.paralog_sim_threshold != null ? "--paralog_sim_threshold ${shellQuote(params.paralog_sim_threshold)}" : ''
    def metricArg = params.metric != null ? "--metric ${shellQuote(params.metric)}" : ''
    def linkageMethod = params.linkage_method != null ? "--linkage_method ${shellQuote(params.linkage_method)}" : ''
    def nxThreshold = params.network_similarity_threshold != null ? "--network_similarity_threshold ${shellQuote(params.network_similarity_threshold)}" : ''

    def clusterOptArgs = [
        optArg('--feature_mode', params.feature_mode),
        optArg('--position_weight', params.position_weight),
        optArg('--synteny_weight', params.synteny_weight),
        optArg('--adjacency_across_gaps', params.adjacency_across_gaps),
        optArg('--collapse_overlaps', params.collapse_overlaps),
        optArg('--overlap_fraction', params.overlap_fraction),
        optArg('--rare_unit', params.rare_unit),
        optArg('--min_features_per_neighbourhood', params.min_features_per_neighbourhood),
        optArg('--min_cluster_size', params.min_cluster_size),
        optArg('--min_coverage', params.min_coverage),
        optArg('--min_cluster_cosine', params.min_cluster_cosine),
        optArg('--cohesion_stat', params.cohesion_stat),
        optArg('--min_shared_families', params.min_shared_families),
        optArg('--signature_fraction', params.signature_fraction),
        optArg('--assign_leftovers', params.assign_leftovers),
        optArg('--assign_min_similarity', params.assign_min_similarity),
        optArg('--supercluster_similarity', params.supercluster_similarity),
        optArg('--keep_small_hier_clusters', params.keep_small_hier_clusters),
        optArg('--max_hierarchical_samples', params.max_hierarchical_samples),
        optArg('--max_heatmap_features', params.max_heatmap_features),
        optArg('--density_clustering', params.density_clustering),
        optArg('--hdbscan_selection', params.hdbscan_selection),
        optArg('--min_samples', params.min_samples),
        optArg('--noise_penalty', params.noise_penalty),
        optArg('--max_hdbscan_candidates', params.max_hdbscan_candidates),
        optArg('--max_silhouette_samples', params.max_silhouette_samples),
        optArg('--logo_for_noise', params.logo_for_noise)
    ].findAll { it }.join(' ')

    """
    mkdir -p density hierarchical

    if [ "${params.microsynteny_logo?.toString()?.equalsIgnoreCase('yes')}" = "true" ]; then
        mkdir -p hierarchical/microsynteny_logo density/microsynteny_logo
    fi

    python ${shellQuote(plot_script)} \\
        --input_file ${shellQuote(merged_with_color)} \\
        --output_path ./ \\
        ${geneNameArg} \\
        ${scaleArg} \\
        ${cutHeightArg} \\
        --gene_of_interest ${shellQuote(params.gene_of_interest)} \\
        ${nameFileArg} \\
        ${logoArg} \\
        ${svgArg} \\
        --goi_type ${shellQuote(params.search_types)} \\
        ${removeCoreArg} \\
        ${removeRareArg} \\
        ${coreThresholdArg} \\
        ${rareMinCountArg} \\
        ${metricArg} \\
        ${minGenePrevalence} \\
        ${linkageMethod} \\
        ${paraRe} \\
        ${paraSim} \\
        ${clusterOptArgs} \\
        ${nxThreshold} \\
        --color_scheme_file ${shellQuote(color_scheme_file)}
    """
}

// -----------------------------------------------------------------------------
// Extract cluster-specific multi-FASTA files for density clustering
// -----------------------------------------------------------------------------
process extractDensityClusterFastas {
    publishDir "${params.output_dir}/4_plot/density/fasta", mode: 'copy'

    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 2
    cpus 1
    memory { 2.GB * task.attempt }
    time { 30.min * task.attempt }

    input:
    path goi_tsv
    path protein_fastas
    path ncrna_fastas
    path extract_script

    output:
    path 'cluster_*_ncrna.mfna', optional: true
    path 'cluster_*_protein.mfaa', optional: true

    script:
    """
    mkdir -p density/fasta
    python ${shellQuote(extract_script)} \\
        --input_tsv ${shellQuote(goi_tsv)} \\
        --neighbour_dir ./ \\
        --output_dir ./
    """
}

// -----------------------------------------------------------------------------
// Extract cluster-specific multi-FASTA files for hierarchical clustering
// -----------------------------------------------------------------------------
process extractHierarchicalClusterFastas {
    publishDir "${params.output_dir}/4_plot/hierarchical/fasta", mode: 'copy'

    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 2
    cpus 1
    memory { 2.GB * task.attempt }
    time { 30.min * task.attempt }

    input:
    path goi_tsv
    path protein_fastas
    path ncrna_fastas
    path extract_script

    output:
    path 'cluster_*_ncrna.mfna', optional: true
    path 'cluster_*_protein.mfaa', optional: true

    script:
    """
    mkdir -p hierarchical/fasta
    python ${shellQuote(extract_script)} \\
        --input_tsv ${shellQuote(goi_tsv)} \\
        --neighbour_dir ./ \\
        --output_dir ./
    """
}

// -----------------------------------------------------------------------------
// Record run parameters
// -----------------------------------------------------------------------------
process writeProvenance {
    publishDir "${params.output_dir}/0_metadata", mode: 'copy'

    cpus 1
    memory '256 MB'
    time '5m'

    input:
    val provenance_json

    output:
    path 'run_provenance.json'
    path 'software_versions.txt'

    script:
    """
    printf '%s\\n' ${shellQuote(provenance_json)} > run_provenance.json

    {
        echo 'blastn:'
        if command -v blastn >/dev/null 2>&1; then blastn -version 2>&1 | head -n 1 || true; else echo 'not available in this task environment'; fi
        echo 'makeblastdb:'
        if command -v makeblastdb >/dev/null 2>&1; then makeblastdb -version 2>&1 | head -n 1 || true; else echo 'not available in this task environment'; fi
        echo 'mmseqs:'
        if command -v mmseqs >/dev/null 2>&1; then mmseqs version 2>&1 | head -n 1 || true; else echo 'not available in this task environment'; fi
        echo 'hmmscan:'
        if command -v hmmscan >/dev/null 2>&1; then hmmscan -h 2>&1 | head -n 2 || true; else echo 'not available in this task environment'; fi
        echo 'cmscan:'
        if command -v cmscan >/dev/null 2>&1; then cmscan -h 2>&1 | head -n 2 || true; else echo 'not available in this task environment'; fi
        echo 'cmsearch:'
        if command -v cmsearch >/dev/null 2>&1; then cmsearch -h 2>&1 | head -n 2 || true; else echo 'not available in this task environment'; fi
        echo 'python:'
        if command -v python >/dev/null 2>&1; then python --version 2>&1 | head -n 1 || true; else echo 'not available in this task environment'; fi
    } > software_versions.txt
    """
}

// -----------------------------------------------------------------------------
// Workflow
// -----------------------------------------------------------------------------
workflow {
    def searchType = (params.search_types ?: '').toString().trim().toLowerCase()
    def validSearchTypes = ['blastn', 'blastp', 'infernal', 'tblastn', 'from_gff']
    if (!(searchType in validSearchTypes)) {
        error "Invalid search type '${params.search_types}'. Supported values: ${validSearchTypes.join(', ')}"
    }

    def skipToClustering = isTrue(params.skip_to_clustering)
    def clusteringParts = (params.adjacent_gene_clustering ?: '').toString()
        .split(/\s*,\s*/)
        .collect { it.trim().toLowerCase() }
    if (clusteringParts.size() != 2) {
        error "Parameter adjacent_gene_clustering must contain exactly two comma-separated values: protein_tool,ncrna_tool (e.g. mmseqs,cmscan)."
    }
    def proteinTool = clusteringParts[0]
    def ncrnaTool = clusteringParts[1]
    if (!(proteinTool in ['mmseqs', 'hmmscan'])) {
        error "Unsupported protein clustering/annotation tool '${proteinTool}'. Supported: mmseqs, hmmscan."
    }
    if (!(ncrnaTool in ['mmseqs', 'cmscan'])) {
        error "Unsupported ncRNA clustering/annotation tool '${ncrnaTool}'. Supported: mmseqs, cmscan."
    }

    // ---- clustering parameter validation (fail before any search is run) ----
    def oneOf = { String name, Collection allowed ->
        def v = params[name]
        if (v != null && !(v.toString().toLowerCase() in allowed))
            error "Invalid --${name} '${v}'. Supported values: ${allowed.join(', ')}"
    }
    oneOf('density_clustering', ['yes', 'no', 'only'])
    oneOf('hdbscan_selection', ['eom', 'leaf'])
    oneOf('feature_mode', ['content', 'combined', 'synteny'])
    oneOf('metric', ['cosine', 'euclidean'])
    oneOf('linkage_method', ['average', 'ward'])
    oneOf('cohesion_stat', ['worst', 'median', 'mean'])
    oneOf('rare_unit', ['genome', 'neighbourhood'])
    ['assign_leftovers', 'adjacency_across_gaps', 'collapse_overlaps', 'keep_small_hier_clusters',
     'logo_for_noise', 'remove_core_genes', 'remove_rare_genes', 'microsynteny_logo'].each { oneOf(it, ['yes', 'no']) }
    if ((params.linkage_method ?: 'average').toString().toLowerCase() == 'ward' &&
        (params.metric ?: 'cosine').toString().toLowerCase() != 'euclidean') {
        error "--linkage_method ward requires --metric euclidean (current metric: ${params.metric ?: 'cosine'})."
    }
    ['min_coverage', 'min_cluster_cosine', 'assign_min_similarity', 'signature_fraction'].each { name ->
        if (params[name] != null) {
            def x = params[name].toString().isNumber() ? params[name].toString().toDouble() : null
            if (x == null || x < 0 || x > 1) error "--${name} must be a number between 0 and 1 (got '${params[name]}')."
        }
    }
    if (params.min_gene_prevalence != null &&
        (params.remove_core_genes?.toString()?.equalsIgnoreCase('yes') || params.remove_rare_genes?.toString()?.equalsIgnoreCase('yes'))) {
        log.warn "--min_gene_prevalence is set: --remove_core_genes / --remove_rare_genes are ignored by the clustering script."
    }
    if (params.cut_height_args != null && params.cut_height_args.toString().isNumber() && params.cut_height_args.toString().toDouble() >= 1) {
        log.warn "--cut_height_args=${params.cut_height_args} is a very high cut; most neighbourhoods will end up in one cluster. Omit it (or use 0) for the automatic cut."
    }

    if (!params.output_dir) error 'Required parameter --output_dir was not set.'
    if (!params.gene_of_interest) error 'Required parameter --gene_of_interest was not set.'
    if (!params.neighbours) error 'Required parameter --neighbours was not set.'
    if (!params.cpus || (params.cpus as int) < 1) error 'Parameter --cpus must be a positive integer.'

    if (params.name_file) {
        def nameFile = file(params.name_file).toAbsolutePath()
        if (!nameFile.toFile().isFile()) error "The supplied name_file does not exist or is not a file: ${nameFile}"
    }

    if (!skipToClustering) {
        if (!params.genomes_dir) error 'Required parameter --genomes_dir was not set.'
        def genomesDir = file(params.genomes_dir).toAbsolutePath()
        if (!genomesDir.toFile().isDirectory()) error "Genomes directory does not exist: ${genomesDir}"

        if (searchType != 'from_gff') {
            if (!params.query) error "Parameter --query is required when search_types='${searchType}'."
            def queryFile = file(params.query).toAbsolutePath()
            if (!queryFile.toFile().isFile()) error "Query file does not exist or is not a file: ${queryFile}"
        }
        if (proteinTool == 'hmmscan') {
            if (!params.pfam_db) error 'Parameter --pfam_db is required when protein tool is hmmscan.'
            def db = file(params.pfam_db).toAbsolutePath()
            if (!db.toFile().isFile()) error "Pfam HMM database does not exist or is not a file: ${db}"
            ['h3f', 'h3i', 'h3m', 'h3p'].each { suffix ->
                def indexFile = file("${db}.${suffix}")
                if (!indexFile.toFile().isFile()) error "Required Pfam index file is missing: ${indexFile}"
            }
        }
        if (ncrnaTool == 'cmscan') {
            if (!params.rfam_db) error 'Parameter --rfam_db is required when ncRNA tool is cmscan.'
            def db = file(params.rfam_db).toRealPath()
            if (!db.toFile().isFile()) error "Rfam CM database does not exist or is not a file: ${db}"
            ['i1f', 'i1i', 'i1m', 'i1p'].each { suffix ->
                def indexFile = file("${db}.${suffix}")
                if (!indexFile.toFile().isFile()) error "Required Rfam index file is missing: ${indexFile}"
            }
        }
    }
    else {
        if (!params.prev_nb_dir) error 'Parameter --prev_nb_dir is required when --skip_to_clustering is true.'
        def prevDir = file(params.prev_nb_dir).toAbsolutePath()
        if (!prevDir.toFile().isDirectory()) error "Previous neighbourhood directory does not exist: ${prevDir}"
        if (proteinTool == 'hmmscan') {
            if (!params.pfam_db) error 'Parameter --pfam_db is required when protein tool is hmmscan.'
            def db = file(params.pfam_db).toAbsolutePath()
            if (!db.toFile().isFile()) error "Pfam HMM database does not exist or is not a file: ${db}"
            ['h3f', 'h3i', 'h3m', 'h3p'].each { suffix ->
                def indexFile = file("${db}.${suffix}")
                if (!indexFile.toFile().isFile()) error "Required Pfam index file is missing: ${indexFile}"
            }
        }
        if (ncrnaTool == 'cmscan') {
            if (!params.rfam_db) error 'Parameter --rfam_db is required when ncRNA tool is cmscan.'
            def db = file(params.rfam_db).toRealPath()
            if (!db.toFile().isFile()) error "Rfam CM database does not exist or is not a file: ${db}"
            ['i1f', 'i1i', 'i1m', 'i1p'].each { suffix ->
                def indexFile = file("${db}.${suffix}")
                if (!indexFile.toFile().isFile()) error "Required Rfam index file is missing: ${indexFile}"
            }
        }
    }

    def neighboursScript = file("${projectDir}/bin/get_neighbours_script.py").toAbsolutePath()
    def mmseqsPostprocessScript = file("${projectDir}/bin/postprocess_mmseqs.py").toAbsolutePath()
    def hmmscanPostprocessScript = file("${projectDir}/bin/postprocess_hmmscan.py").toAbsolutePath()
    def cmscanPostprocessScript = file("${projectDir}/bin/postprocess_cmscan.py").toAbsolutePath()
    def colorScript = file("${projectDir}/bin/color_clusters_script.py").toAbsolutePath()
    def plotScript = file("${projectDir}/bin/plot_context_script.py").toAbsolutePath()

    def requiredHelperScripts = [mmseqsPostprocessScript, colorScript, plotScript]
    if (!skipToClustering) requiredHelperScripts << neighboursScript
    if (proteinTool == 'hmmscan') requiredHelperScripts << hmmscanPostprocessScript
    if (ncrnaTool == 'cmscan') requiredHelperScripts << cmscanPostprocessScript
    requiredHelperScripts.each { scriptFile ->
        if (!scriptFile.toFile().isFile()) error "Required pipeline helper script is missing: ${scriptFile}"
    }

    def tsv_files
    def mfaa_files
    def mfna_files
    if (skipToClustering) {
        log.info "Skipping search/neighbourhood stages; loading previous results from ${file(params.prev_nb_dir).toAbsolutePath()}"
        def previousDir = file(params.prev_nb_dir).toAbsolutePath()
        def previousTsv = Channel.fromPath("${previousDir}/*.nb.tsv")
            .ifEmpty { error "No *.nb.tsv files found in --prev_nb_dir ${previousDir}" }
        def previousProtein = Channel.fromPath("${previousDir}/*.nb.protein.mfaa")
            .ifEmpty { error "No *.nb.protein.mfaa files found in --prev_nb_dir ${previousDir}; clustering requires protein neighbourhood FASTA files." }
        def previousNcrna = Channel.fromPath("${previousDir}/*.nb.ncrna.mfna")
            .ifEmpty { error "No *.nb.ncrna.mfna files found in --prev_nb_dir ${previousDir}; clustering requires ncRNA neighbourhood FASTA files." }

        tsv_files = previousTsv.collect()
        mfaa_files = previousProtein.collect()
        mfna_files = previousNcrna.collect()
    }
    else {
        def genomesDir = file(params.genomes_dir).toAbsolutePath()
        def genome_gff_pairs = Channel
            .fromPath("${genomesDir}/*", type: 'dir')
            .map { subfolder ->
                def directory = new java.io.File(subfolder.toString())
                def listedFiles = directory.listFiles()
                if (listedFiles == null) {
                    error "Cannot read genome directory '${directory}'. Check permissions and confirm this is a readable local directory."
                }

                def files = listedFiles.findAll { f -> f.isFile() }
                def fnaCandidates = files.findAll { f ->
                    def lower = f.name.toLowerCase(java.util.Locale.ROOT)
                    lower.endsWith('.fasta') || lower.endsWith('.fna') || lower.endsWith('.fa')
                }
                def gffCandidates = files.findAll { f ->
                    def lower = f.name.toLowerCase(java.util.Locale.ROOT)
                    lower.endsWith('.gff') || lower.endsWith('.gff3')
                }
                if (fnaCandidates.size() != 1 || gffCandidates.size() != 1) {
                    def entryNames = listedFiles.collect { it.name }.sort().join(', ')
                    error "Genome folder '${directory}' must contain exactly one FASTA (.fasta/.fna/.fa) and one uncompressed GFF (.gff/.gff3); found ${fnaCandidates.size()} FASTA and ${gffCandidates.size()} GFF files. Directory entries seen: [${entryNames}]."
                }
                tuple(subfolder.fileName.toString(), fnaCandidates[0].toPath(), gffCandidates[0].toPath())
            }
            .ifEmpty { error "No genome subfolders found in ${genomesDir}. Each folder must contain one FASTA and one GFF file." }

        if (searchType == 'from_gff') {
            def ready_for_neighbours = genome_gff_pairs.map { id, genome, gff ->
                tuple(id, genome, gff, genome, false)
            }
            def neighbour_results = getNeighbours(ready_for_neighbours, Channel.value(neighboursScript))

            tsv_files = neighbour_results.tsv
                .ifEmpty { error 'Neighbourhood extraction produced no required .nb.tsv output files.' }
                .map { id, path -> path }
                .collect()
            mfaa_files = neighbour_results.mfaa
                .ifEmpty { error 'Neighbourhood extraction produced no protein FASTA files; protein clustering cannot continue.' }
                .collect()
            mfna_files = neighbour_results.mfna
                .ifEmpty { error 'Neighbourhood extraction produced no ncRNA FASTA files; ncRNA clustering cannot continue.' }
                .collect()
        }
        else {
            def queryFile = file(params.query).toAbsolutePath()
            def search_results = runSearch(genome_gff_pairs, Channel.value(queryFile)).search_results
                .map { id, result -> tuple(id, result) }
                .filter { id, result ->
                    searchType == 'infernal' ? hasNonCommentRow(result) : result.size() > 0
                }
                .ifEmpty { error "No search hits found for search type '${searchType}'. Check the query, FASTA type and search parameters." }

            def ready_for_neighbours = genome_gff_pairs
                .join(search_results)
                .map { id, genome, gff, result -> tuple(id, genome, gff, result, true) }

            def neighbour_results = getNeighbours(ready_for_neighbours, Channel.value(neighboursScript))

            tsv_files = neighbour_results.tsv
                .ifEmpty { error 'Neighbourhood extraction produced no required .nb.tsv output files.' }
                .map { id, path -> path }
                .collect()
            mfaa_files = neighbour_results.mfaa
                .ifEmpty { error 'Neighbourhood extraction produced no protein FASTA files; protein clustering cannot continue.' }
                .collect()
            mfna_files = neighbour_results.mfna
                .ifEmpty { error 'Neighbourhood extraction produced no ncRNA FASTA files; ncRNA clustering cannot continue.' }
                .collect()
        }
    }

    def clustering_results_ch = Channel.empty()

    if (proteinTool == 'hmmscan') {
        def pfamHMM = file(params.pfam_db).toAbsolutePath()
        def pfamH3f = file("${pfamHMM}.h3f").toAbsolutePath()
        def pfamH3i = file("${pfamHMM}.h3i").toAbsolutePath()
        def pfamH3m = file("${pfamHMM}.h3m").toAbsolutePath()
        def pfamH3p = file("${pfamHMM}.h3p").toAbsolutePath()

        def mmseqsProt = runMMseqsProtein(mfaa_files, Channel.value(mmseqsPostprocessScript))
        def hmmInput = mmseqsProt.mmseqs_out.map { cluster_map, reps -> tuple(reps, cluster_map) }
        def hmmOut = runHMMscan(
            hmmInput,
            Channel.value(pfamHMM),
            Channel.value(pfamH3f),
            Channel.value(pfamH3i),
            Channel.value(pfamH3m),
            Channel.value(pfamH3p),
            Channel.value(hmmscanPostprocessScript)
        )
        clustering_results_ch = clustering_results_ch.mix(hmmOut.hmm_results)
    }
    else if (proteinTool == 'mmseqs') {
        def mmseqsProt = runMMseqsProtein(mfaa_files, Channel.value(mmseqsPostprocessScript))
        clustering_results_ch = clustering_results_ch.mix(mmseqsProt.mmseqs_out.map { cluster_map, reps -> cluster_map })
    }

    if (ncrnaTool == 'cmscan') {
        def mmseqsNc = runMMseqsNCRNA(mfna_files, Channel.value(mmseqsPostprocessScript))
        def cmInput = mmseqsNc.mmseqs_out.map { cluster_map, reps -> tuple(reps, cluster_map) }

        def rfamDb = file(params.rfam_db).toRealPath()
        def rfamI1f = file("${rfamDb}.i1f").toAbsolutePath()
        def rfamI1i = file("${rfamDb}.i1i").toAbsolutePath()
        def rfamI1m = file("${rfamDb}.i1m").toAbsolutePath()
        def rfamI1p = file("${rfamDb}.i1p").toAbsolutePath()

        [rfamDb, rfamI1f, rfamI1i, rfamI1m, rfamI1p].each { dbFile ->
            if (!dbFile.toFile().isFile()) error "Required Rfam database/index file is missing: ${dbFile}"
        }

        def cmOut = runCMscan(
            cmInput,
            Channel.value(rfamDb),
            Channel.value(rfamI1f),
            Channel.value(rfamI1i),
            Channel.value(rfamI1m),
            Channel.value(rfamI1p),
            Channel.value(cmscanPostprocessScript)
        )
        clustering_results_ch = clustering_results_ch.mix(cmOut.cm_results)
    }
    else if (ncrnaTool == 'mmseqs') {
        def mmseqsNc = runMMseqsNCRNA(mfna_files, Channel.value(mmseqsPostprocessScript))
        clustering_results_ch = clustering_results_ch.mix(mmseqsNc.mmseqs_out.map { cluster_map, reps -> cluster_map })
    }

    def combined_clusters_ch = clustering_results_ch.collect()
    def colorSchemeFile = file(params.color_scheme_file).toAbsolutePath()

    if (!colorSchemeFile.toFile().isFile()) {
        error "Color scheme file not found: ${colorSchemeFile}"
    }

    def colored_ch = clusterColoring(
        combined_clusters_ch,
        tsv_files,
        Channel.value(colorScript),
        Channel.value(plotScript),
        colorSchemeFile
    )

    def nameFilePath = params.name_file ? file(params.name_file).toAbsolutePath() : null
    def plotInputs = colored_ch.colored_tsv.map { coloredFile ->
        tuple(coloredFile, nameFilePath ?: coloredFile, nameFilePath != null)
    }

    def plotOutputs = plottingContext(
        plotInputs,
        Channel.value(plotScript),
        colorSchemeFile
    )

    def extractFastasScript = file("${projectDir}/bin/get_cluster_fastas.py").toAbsolutePath()
    if (!extractFastasScript.toFile().isFile()) {
        error "Required helper script is missing: ${extractFastasScript}"
    }

    // Only run an extraction when that method produced a gene_of_interest.tsv
    // (--density_clustering only has no hierarchical result; 'no' has no density result).
    extractDensityClusterFastas(
        plotOutputs.density_dir.map { dir -> file("${dir}/gene_of_interest.tsv") }.filter { it.exists() },
        mfaa_files,
        mfna_files,
        Channel.value(extractFastasScript)
    )

    extractHierarchicalClusterFastas(
        plotOutputs.hierarchical_dir.map { dir -> file("${dir}/gene_of_interest.tsv") }.filter { it.exists() },
        mfaa_files,
        mfna_files,
        Channel.value(extractFastasScript)
    )

    def provenance = [
        pipeline: 'SweetSynteny',
        pipeline_version: version,
        nextflow_version: workflow.nextflow.version.toString(),
        run_name: workflow.runName,
        session_id: workflow.sessionId.toString(),
        revision: workflow.revision ?: 'unknown',
        command_line: workflow.commandLine,
        launch_directory: workflow.launchDir.toString(),
        work_directory: workflow.workDir.toString(),
        launched_at: new Date().format("yyyy-MM-dd'T'HH:mm:ssXXX"),
        search_mode: searchType,
        clustering_mode: [protein: proteinTool, ncrna: ncrnaTool],
        inputs: [
            genomes_dir: params.genomes_dir ? file(params.genomes_dir).toAbsolutePath().toString() : null,
            query: (params.query && searchType != 'from_gff') ? file(params.query).toAbsolutePath().toString() : null,
            previous_neighbourhood_dir: (skipToClustering && params.prev_nb_dir) ? file(params.prev_nb_dir).toAbsolutePath().toString() : null,
            pfam_db: (proteinTool == 'hmmscan' && params.pfam_db) ? file(params.pfam_db).toAbsolutePath().toString() : null,
            rfam_db: (ncrnaTool == 'cmscan' && params.rfam_db) ? file(params.rfam_db).toRealPath().toString() : null,
            name_file: params.name_file ? file(params.name_file).toAbsolutePath().toString() : null
        ],
        parameters: [
            output_dir: params.output_dir,
            bio_type: params.bio_type,
            annotation_type: params.annotation_type,
            gene_of_interest: params.gene_of_interest,
            neighbours: params.neighbours,
            including_features: params.including_features,
            ignore_overlaps: params.ignore_overlaps,
            substring_search: params.substring_search,
            promoter: params.promoter,
            promoter_len: params.promoter_len,
            promoter_mode: params.promoter_mode,
            overlap_threshold: params.overlap_threshold,
            evalue_threshold_blast: params.evalue_threshold_blast,
            evalue_threshold_infernal: params.evalue_threshold_infernal,
            len_threshold_blast: params.len_threshold_blast,
            len_threshold_infernal: params.len_threshold_infernal,
            blast_matrix: params.blast_matrix,
            cpus: params.cpus,
            scale: params.scale,
            cut_height_args: params.cut_height_args,
            microsynteny_logo: params.microsynteny_logo,
            svg: params.svg,
            remove_core_genes: params.remove_core_genes,
            remove_rare_genes: params.remove_rare_genes,
            core_threshold: params.core_threshold,
            rare_min_count: params.rare_min_count,
            min_gene_prevalence: params.min_gene_prevalence,
            paralog_resolution: params.paralog_resolution,
            paralog_sim_threshold: params.paralog_sim_threshold,
            metric: params.metric,
            linkage_method: params.linkage_method,
            feature_mode: params.feature_mode,
            position_weight: params.position_weight,
            synteny_weight: params.synteny_weight,
            adjacency_across_gaps: params.adjacency_across_gaps,
            collapse_overlaps: params.collapse_overlaps,
            overlap_fraction: params.overlap_fraction,
            rare_unit: params.rare_unit,
            min_features_per_neighbourhood: params.min_features_per_neighbourhood,
            min_cluster_size: params.min_cluster_size,
            min_coverage: params.min_coverage,
            min_cluster_cosine: params.min_cluster_cosine,
            cohesion_stat: params.cohesion_stat,
            min_shared_families: params.min_shared_families,
            signature_fraction: params.signature_fraction,
            assign_leftovers: params.assign_leftovers,
            assign_min_similarity: params.assign_min_similarity,
            supercluster_similarity: params.supercluster_similarity,
            keep_small_hier_clusters: params.keep_small_hier_clusters,
            max_hierarchical_samples: params.max_hierarchical_samples,
            max_heatmap_features: params.max_heatmap_features,
            density_clustering: params.density_clustering,
            hdbscan_selection: params.hdbscan_selection,
            min_samples: params.min_samples,
            noise_penalty: params.noise_penalty,
            max_hdbscan_candidates: params.max_hdbscan_candidates,
            max_silhouette_samples: params.max_silhouette_samples,
            logo_for_noise: params.logo_for_noise,
            skip_to_clustering: skipToClustering
        ]
    ]

    writeProvenance(Channel.value(JsonOutput.prettyPrint(JsonOutput.toJson(provenance))))
}

workflow.onError {
    log.error "Pipeline failed with error: ${workflow.errorMessage}"
    log.error "Check the work directory for failed tasks: ${workflow.workDir}"
}

workflow.onComplete {
    if (workflow.success) {
        log.info 'Pipeline completed successfully!'
        log.info "Results are in: ${params.output_dir}"
    }
    else {
        log.error 'Pipeline failed. Check logs above.'
    }
    log.info "Workflow finished at: ${new Date()}"
}
