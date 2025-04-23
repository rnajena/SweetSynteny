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
     Search Type      : ${params.types} [blastn, blastp, infernal, tblastn]
     Database         : ${params.genomes_dir}
     Query            : ${params.query}
    >Parameter for Neighbours
     Gene of interest : ${params.gene_of_interest} [Target gene identifie]
     Neighbours : ${params.neighbours} [Neighbor range: x,y (genes) or x:y (nucleotides)]
    >Parameter for Plotting
     Clustering level : ${params.cluster_level} [sequence_level or in futur seconday_level]
     Scale : ${params.scale}
     Plotting ending : ${params.plotting} [.png or .svg]
     Cluster : ${params.cluster} [Minimal size for a cluster, default 2]
     Threshold : ${params.threshold} [Similarity threshold for clustering, default 0.3]
    """
    .stripIndent(true)

// Process to perform sequence search using BLAST or Infern
process runSearch {
    publishDir "${params.output_dir}/search_results", mode: 'copy', pattern: '*_search_result.tsv'

    input:
    tuple val(id), path(genome), path(gff)
    
    output:
    tuple val(id), path(genome), path(gff), path("${id}_search_result.tsv")
    
    script:
    if (params.types == 'blastn')
        """
        makeblastdb -in $genome -dbtype nucl
        blastn \\
            -num_threads $task.cpus \\
            -query ${params.query} \\
            -subject $genome \\
            -out ${id}_search_result.tsv \\
            -outfmt "6 qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq" \\
            -evalue 0.01
        """
    else if (params.types == 'blastp')
        """
        makeblastdb -in $genome -dbtype nucl
        blastp \\
            -num_threads $task.cpus \\
            -query ${params.query} \\
            -db $genome \\
            -out ${id}_search_result.tsv \\
            -outfmt "6 qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq" \\
            -evalue 0.01
        """
    else if (params.types == 'tblastn')
        """
        makeblastdb -in $genome -dbtype nucl
        tblastn \\
            -num_threads $task.cpus \\
            -query ${params.query} \\
            -db $genome \\
            -out ${id}_search_result.tsv \\
            -outfmt "6 qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq" \\
            -evalue 0.01
        """
    else if (params.types == 'infernal')
        """
        cmsearch \\
            --cpu $task.cpus \\
            --tblout ${id}_search_result.tsv \\
            ${params.query} \\
            $genome
        """
    else
        error "Invalid search type: ${params.types}"
}

// Process to identify neighboring genes
process getNeighbours {
    publishDir "${params.output_dir}/neighbour_results", mode: 'copy', pattern: '*_{neighbours_output.tsv,neighbours_output.srna.mfna,neighbours_output.protein.mfaa}'

    input:
    tuple val(id), path(genome), path(gff), path(search_result)
    
    output:
    tuple val(id), path("${id}_neighbours_output.tsv"), path("${id}_neighbours_output.protein.mfaa"), path("${id}_neighbours_output.srna.mfna")
    
    script:
    """
    # TODO mkdir neighbour_results
    python ${projectDir}/bin/get_neighbours_script.py \\
        --hit_file $search_result \\
        --input_type ${params.types} \\
        --fna_file $genome \\
        --gff_file $gff \\
        --gene_of_interest ${params.gene_of_interest} \\
        --neighbours ${params.neighbours} \\
        --output_path ${id}_neighbours_output
    """
}

// Process for clustering and coloring results
process clusterColoring {
    publishDir "${params.output_dir}/clustered_results", mode: 'copy', pattern: 'merged_neighbours*'
    publishDir "${params.output_dir}/clustered_results", mode: 'copy', pattern: 'merged_with_color.tsv'

    input:
    path tsv_files
    path mfna_files
    path mfaa_files

    output:
    path "merged_neighbours.tsv"
    path "merged_neighbours.protein.mfaa"
    path "merged_neighbours.srna.mfna"
    path "clusterRes_cluster.tsv"
    path "merged_with_color.tsv", emit: colored_tsv

    script:
    """
    # Merge results
    cat $tsv_files > merged_neighbours.tsv
    cat $mfna_files > merged_neighbours.srna.mfna
    cat $mfaa_files > merged_neighbours.protein.mfaa

    # Cluster proteins
    mmseqs createdb \\
        --dbtype 0 merged_neighbours.protein.mfaa \\
        merged_db.protein
    mkdir tmp.protein
    mmseqs easy-linclust \\
        merged_neighbours.protein.mfaa \\
        clusterRes.protein \\
        tmp.protein \\
        --min-seq-id 0.5 \\
        -c 0.8 \\
        --cov-mode 1 \\
        --cluster-mode 2

    # Cluster ncRNAs
    if [ "${params.cluster_level}" == "sequence_level" ]; then
        mmseqs createdb \\
            --dbtype 2 \\
            merged_neighbours.srna.mfna \\
            merged_db.srna
        mkdir tmp.srna
        mmseqs easy-linclust \\
            merged_neighbours.srna.mfna \\
            clusterRes.srna tmp.srna \\
            --min-seq-id 0.5 \\
            -c 0.8 \\
            --cov-mode 0 \\
            --cluster-mode 2 -k 6 
    elif [ "${params.cluster_level}" == "secondary_level" ]; then
        echo "RFAM clustering not implemented yet" # > clusterRes.srna
        # TODO
    else
        echo "Invalid clustering type: ${params.cluster_level}" >&2
        exit 1
    fi

    cat *_cluster.tsv > clusterRes_cluster.tsv

    python ${projectDir}/bin/color_clusters_script.py \\
        --cluster_file clusterRes_cluster.tsv \\
        --tsv_file merged_neighbours.tsv \\
        --output_file merged_with_color.tsv
    """
}

// Process for plotting genomic context
process plottingContext {
    publishDir "${params.output_dir}/plot_results", mode: 'copy', pattern: 'plot*'

    input:
    path merged_with_color

    output:
    path "plot*"

    script:
    """
    python ${projectDir}/bin/plot_context_script.py \\
        --input_file $merged_with_color \\
        --output_path plot \\
        --scale ${params.scale} \\
        --output_ending ${params.plotting} \\
        --cluster ${params.cluster} \\
        --threshold ${params.threshold}
    """
}

workflow {
    // Create channel of genome-GFF pairs
    genome_gff_pairs = Channel
        .fromPath("${params.genomes_dir}/*", type: 'dir')
        .map { subfolder -> 
            def fna = subfolder.listFiles().find { it.name.endsWith('.fasta') }
            def gff = subfolder.listFiles().find { it.name.endsWith('.gff') }
            tuple(subfolder.name, fna, gff)
        }
        .filter { it[1] != null && it[2] != null }
    // Run search process
    search_results = runSearch(genome_gff_pairs)
    // Filter out empty results
    non_empty_results = search_results.filter { id, genome, gff, result ->
        result.size() > 0
    }
    // Get neighboring genes
    neighbour_results = getNeighbours(non_empty_results)
    // Collect TSV and MFNA files separately
    tsv_files = neighbour_results.map { it[1] }.collect()
    mfna_files = neighbour_results.map { it[2] }.collect()
    mfaa_files = neighbour_results.map { it[3] }.collect()
    // Perform clustering and coloring
    clusterColoring(tsv_files, mfna_files, mfaa_files)
    // Generate context plot
    plottingContext(clusterColoring.out.colored_tsv)
}   
log.info "Pipeline completed at: $workflow.complete"
