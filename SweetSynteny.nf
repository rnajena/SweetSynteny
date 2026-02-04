#!/usr/bin/env nextflow

/* 
* SweetSynteny
*
* Authors: 
* - Maria Schreiber <maria.schreiber@uni-jena.de>
*/

nextflow.enable.dsl=2

// Pipeline version
version = '1.0'

log.info """
    SweetSyntheny - NF Pipeline
    ===========================
    >Results
     Result Folder    : ${params.output_dir}
    >Parameter for Searching
     Search Type      : ${params.search_types} [blastn, blastp, infernal, tblastn]
     Bio Type         : ${params.bio_type}
     Database         : ${params.genomes_dir}
     Query            : ${params.query}
     Annotation       : ${params.annotation_type}
    >Parameter for Neighbours
     Gene of interest : ${params.gene_of_interest} [Target gene identifier]
     Neighbours : ${params.neighbours} [Neighbor range: x,y (genes) or x:y (nucleotides)]
     Include GFF features: ${params.including_features}
     Ignore overlap filter : ${params.ignore_overlaps}
     Substring search : ${params.substring_search} [Only relevant for from_gff search]
    >Clustering 
     For adjacent genes: ${params.adjacent_gene_clustering} [mmseqs,mmseqs|mmseqs,cmscan|hmmscan,mmseqs|hmmscan,cmscan]
    >Parameter for Plotting
     Scale : ${params.scale}
     For DBscan
     Cluster : ${params.cluster} [Minimal size for a cluster, default 2]
     Threshold : ${params.threshold} [Similarity threshold for clustering, default 0.3]
     Microsynteny logo : ${params.microsynteny_logo}
     For H-clustering
     Threshold : ${params.cut_height_args}
    >CPU : ${params.cpus}
    """
    .stripIndent(true)

// Process to perform sequence search using BLAST or Infernal
process runSearch {
    publishDir "${params.output_dir}/1_search", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(id), path(genome), path(gff)
    
    output:
    tuple val(id), path(genome), path(gff), path("${id}.tsv")
    
    script:
    if (params.search_types == 'blastn')
        """
        if [ ! -e "${genome}.nin" ] && [ ! -e "${genome}.00.nin" ]; then
                makeblastdb -in $genome -dbtype nucl
        fi
        blastn \\
            -num_threads ${params.cpus} \\
            -query ${params.query} \\
            -subject $genome \\
            -out ${id}.tsv \\
            -outfmt "6 qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq" \\
            -evalue ${params.evalue_threshold_blast}
        """
    else if (params.search_types == 'blastp')
        """
        if [ ! -e "${genome}.nin" ] && [ ! -e "${genome}.00.nin" ]; then
                makeblastdb -in $genome -dbtype nucl
        fi
        blastp \\
            -num_threads ${params.cpus} \\
            -query ${params.query} \\
            -db $genome \\
            -out ${id}.tsv \\
            -outfmt "6 qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq" \\
            -evalue ${params.evalue_threshold_blast}
        """
    else if (params.search_types == 'tblastn')
        """
        if [ ! -e "${genome}.nin" ] && [ ! -e "${genome}.00.nin" ]; then
                makeblastdb -in $genome -dbtype nucl
        fi
        tblastn \\
            -num_threads ${params.cpus} \\
            -query ${params.query} \\
            -db $genome \\
            -out ${id}.tsv \\
            -outfmt "6 qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq" \\
            -evalue ${params.evalue_threshold_blast}
        """
    else if (params.search_types == 'infernal')
        """
        cmsearch \\
            --cpu ${params.cpus} \\
            --tblout ${id}.tsv \\
            ${params.query} \\
            $genome
        """
    else
        error "Invalid search type: ${params.search_types}"
}

// Process to identify neighboring genes
process getNeighbours {
    //publishDir "${params.output_dir}/2_neighbour", mode: 'copy', pattern: '*_{neighbours_output.tsv,neighbours_output.ncrna.mfna,neighbours_output.protein.mfaa}'
    publishDir "${params.output_dir}/2_neighbour/", mode: 'copy'

    input:
        tuple val(id), 
        path(genome), 
        path(gff), 
        val(search_result)
    
    output:
        val(id) // To keep track of the ID
        path("${id}.nb.tsv"),          optional: true, emit: tsv
        path("${id}.nb.protein.mfaa"), optional: true, emit: mfaa
        path("${id}.nb.ncrna.mfna"),   optional: true, emit: mfna

    script:
    
    // Convert the Groovy list [gene, ncRNA] into a space-separated string for the shell
    def feature_list = params.including_features.join(' ')

    // Logic: include hit file flag only if search_result is a path/file
    def hit_input = (search_result instanceof Path) ? "--hit_file $search_result" : ""

    // Boolean flags: expand to the whole flag or empty string
    def subSearch = (params.substring_search == true || params.substring_search.toString() == "True") ? "--substring_search" : ""
    def ignoreOver = (params.ignore_overlaps == true || params.ignore_overlaps.toString() == "True") ? "--ignore_overlaps" : ""

    // Always include promoter_mode and overlap_threshold (they have defaults in nextflow.config)
    def promoterMode = "--promoter_mode ${params.promoter_mode}"
    def overlapArg   = "--overlap_threshold ${params.overlap_threshold}"

    """
    mkdir -p "${params.output_dir}/2_neighbour"

    python ${projectDir}/bin/get_neighbours_script.py \\
        --bio_type ${params.bio_type} \\
        --input_type ${params.search_types} \\
        --fna_file $genome \\
        --gff_file $gff \\
        --including_features ${feature_list} \\
        --gene_of_interest ${params.gene_of_interest} \\
        --neighbours ${params.neighbours} \\
        --output_path ${id}.nb \\
        --promoter ${params.promoter} \\
        --promoter_len ${params.promoter_len} \\
        ${promoterMode} \\
        ${subSearch} \\
        ${ignoreOver} \\
        --evalue_threshold_blast ${params.evalue_threshold_blast} \\
        --evalue_threshold_infernal ${params.evalue_threshold_infernal} \\
        --len_threshold_blast ${params.len_threshold_blast} \\
        --len_threshold_infernal ${params.len_threshold_infernal} \\
        --overlap_threshold ${params.overlap_threshold} \\
        ${hit_input}
    """
}

// MMseqs clustering for Proteins
process runMMseqsProtein {
    publishDir "${params.output_dir}/3_cluster", mode: 'copy'
    input: 
        path protein_mfaa
    output: 
        path "clusterRes.protein_mmseqs.tsv", emit: results

    script:
    """
    cat ${protein_mfaa} > merged_protein.mfaa
    
    # MMseqs creates prot_clust_cluster.tsv, prot_clust_rep_seq.fasta, and prot_clust_all_seqs.fasta
    mmseqs easy-linclust merged_protein.mfaa prot_clust tmp --min-seq-id 0.3 --cov-mode 1
    
    # We only pass the _cluster.tsv to the python script
    python ${projectDir}/bin/postprocess_mmseqs.py \
        --input prot_clust_cluster.tsv \
        --output clusterRes.protein_mmseqs.tsv \
        --gene_of_interest ${params.gene_of_interest}
    """
}

// MMseqs clustering for ncRNA
process runMMseqsNCRNA {
    publishDir "${params.output_dir}/3_cluster", mode: 'copy'
    input: 
        path ncrna_mfna
    output: 
        path "clusterRes.ncrna_mmseqs.tsv", emit: results

    script:
    """
    cat ${ncrna_mfna} > merged_ncrna.mfna
    
    # MMseqs creates ncrna_clust_cluster.tsv and others
    mmseqs easy-linclust merged_ncrna.mfna ncrna_clust tmp --min-seq-id 0.8 --cov-mode 1
    
    python ${projectDir}/bin/postprocess_mmseqs.py \
        --input ncrna_clust_cluster.tsv \
        --output clusterRes.ncrna_mmseqs.tsv \
        --gene_of_interest ${params.gene_of_interest}
    """
}

process runHMMscan {
    publishDir "${params.output_dir}/3_cluster", mode: 'copy'
    input:
        path protein_mfaa
    output:
        path "clusterRes.hmm.protein_cluster.tsv", emit: hmm_results

    script:
    """
    cat ${protein_mfaa} > merged_protein.mfaa
    # Use the merged file 'merged_protein.mfaa' here:
    hmmscan --tblout merged.protein.hmm.tbl --cpu 10 ${params.pfam_db} merged_protein.mfaa
    
    python ${projectDir}/bin/postprocess_hmmscan.py \\
        --hmmscan_file merged.protein.hmm.tbl \\
        --output_file clusterRes.hmm.protein_cluster.tsv \\
        --gene_of_interest ${params.gene_of_interest}
    """
}

process runCMscan {
    publishDir "${params.output_dir}/3_cluster", mode: 'copy'
    input:
        path mfna_files // Input variable name
    output:
        path "clusterRes.ncrna_cluster.tsv", emit: cm_results

    script:
    """
    cat ${mfna_files} > merged_ncrna.mfna
    # Use the merged file 'merged_ncrna.mfna' here:
    cmscan -E 0.01 --cpu 10 --noali --tblout clusterRes.ncrna ${params.rfam_db} merged_ncrna.mfna

    python ${projectDir}/bin/postprocess_cmscan.py \\
        --cmscan_file clusterRes.ncrna \\
        --output_file clusterRes.ncrna_cluster.tsv \\
        --gene_of_interest ${params.gene_of_interest}
    """
}

// Process for clustering and coloring results
process clusterColoring {
    publishDir "${params.output_dir}/3_cluster", mode: 'copy'

    input:
        path cluster_files // This contains [clusterRes.hmm.protein_cluster.tsv, clusterRes.ncrna_cluster.tsv]
        path tsv_files     // This contains all the *.nb.tsv files

    output:
        path 'merged_with_color.tsv', emit: colored_tsv
        path "summary.png"
        path "summary.svg"

    script:
    """
    # 1. Merge all neighborhood data staged in the local directory
    # Using 'cat *.nb.tsv' works because Nextflow put them all here
    cat *.nb.tsv | grep -v 'promoter' > merged_neighbours.tsv
    
    # 2. Merge the cluster result files
    # Since 'cluster_files' is a list of two files, we cat them together
    cat ${cluster_files} > clusterRes.tsv

    # 3. Run the coloring script
    python ${projectDir}/bin/color_clusters_script.py \\
        --cluster_file clusterRes.tsv \\
        --tsv_file merged_neighbours.tsv \\
        --output_file merged_with_color.tsv \\
        --gene_of_interest ${params.gene_of_interest}
    """
}

// Process for plotting genomic context
process plottingContext {
    publishDir "${params.output_dir}/4_plot", mode: 'copy'
    
    input:
    path merged_with_color

    output:
    path "dbscan"
    path "ward_cosinus"
    path "ward_euclidean"
    path "*_korrelations_heatmap.png"
    path "*_distribution.png"
    path "clustering_report.txt"

    script:
    """
    mkdir -p dbscan
    mkdir -p ward_cosinus
    mkdir -p ward_euclidean

    python ${projectDir}/bin/plot_context_script.py \\
        --input_file $merged_with_color \\
        --output_path ./ \\
        --scale ${params.scale} \\
        --cluster ${params.cluster} \\
        --threshold ${params.threshold} \\
        --cut_height_args ${params.cut_height_args} \\
        --gene_of_interest ${params.gene_of_interest} \\
        --name_file ${params.name_file} \\
        --method ${params.method} > clustering_report.txt 2>&1
    """
}

workflow {
    // Define the three main data channels
    def tsv_files, mfaa_files, mfna_files

    if (params.skip_to_clustering) {
        log.info "Skipping search/neighbour steps. Loading files from: ${params.prev_nb_dir}"
        
        // Load existing files from the previous output directory
        tsv_files  = Channel.fromPath("${params.prev_nb_dir}/*.nb.tsv").collect()
        mfaa_files = Channel.fromPath("${params.prev_nb_dir}/*.nb.protein.mfaa").collect()
        mfna_files = Channel.fromPath("${params.prev_nb_dir}/*.nb.ncrna.mfna").collect()

    } else {
        // --- ORIGINAL SEARCH/NEIGHBOUR LOGIC ---
        genome_gff_pairs = Channel
            .fromPath("${params.genomes_dir}/*", type: 'dir')
            .map { subfolder -> 
                def fna = subfolder.listFiles().find { it.name.endsWith('.fasta') || it.name.endsWith('.fna') }
                def gff = subfolder.listFiles().find { it.name.endsWith('.gff') }
                if (fna && gff) [subfolder.name, fna, gff]
            }
            .filter { it != null }

        if (params.search_types == 'from_gff') {
            ready_for_neighbours = genome_gff_pairs.map { id, fna, gff -> 
                tuple(id, fna, gff, params.gene_of_interest) 
            }
        } else {
            search_results = runSearch(genome_gff_pairs)
            ready_for_neighbours = search_results.filter { id, g, gf, res -> res.size() > 0 }
        }

        neighbour_results = getNeighbours(ready_for_neighbours)
        
        tsv_files  = neighbour_results.tsv.collect()
        mfaa_files = neighbour_results.mfaa.collect()
        mfna_files = neighbour_results.mfna.collect()
    }

    // --- CLUSTERING LOGIC (Remains the same, but uses the channels defined above) ---
    def tools = params.adjacent_gene_clustering.split(',')
    def protein_tool = tools[0]
    def ncrna_tool   = tools[1]
    def clustering_results_ch = Channel.empty()

    if (protein_tool == 'hmmscan') {
        clustering_results_ch = clustering_results_ch.mix(runHMMscan(mfaa_files).hmm_results)
    } else if (protein_tool == 'mmseqs') {
        clustering_results_ch = clustering_results_ch.mix(runMMseqsProtein(mfaa_files).results)
    }

    if (ncrna_tool == 'cmscan') {
        clustering_results_ch = clustering_results_ch.mix(runCMscan(mfna_files).cm_results)
    } else if (ncrna_tool == 'mmseqs') {
        clustering_results_ch = clustering_results_ch.mix(runMMseqsNCRNA(mfna_files).results)
    }

    combined_clusters_ch = clustering_results_ch.collect()
    colored_ch = clusterColoring(combined_clusters_ch, tsv_files)
    plottingContext(colored_ch.colored_tsv)
}
log.info "Pipeline completed at: $workflow.complete"
