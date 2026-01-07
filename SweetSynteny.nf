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
     Annotation       : ${params.annotation_type}
    >Parameter for Neighbours
     Gene of interest : ${params.gene_of_interest} [Target gene identifier]
     Neighbours : ${params.neighbours} [Neighbor range: x,y (genes) or x:y (nucleotides)]
     Include GFF features: ${params.including_features}
    >Parameter for Plotting
     Scale : ${params.scale}
     Cluster : ${params.cluster} [Minimal size for a cluster, default 2]
     Threshold : ${params.threshold} [Similarity threshold for clustering, default 0.3]
     Microsynteny logo : ${params.microsynteny_logo}
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
    if (params.types == 'blastn')
        """
        if [ ! -e "${genome}.nin" ] && [ ! -e "${genome}.00.nin" ]; then
                makeblastdb -in $genome -dbtype nucl
        fi
        blastn \\
            -num_threads $task.cpus \\
            -query ${params.query} \\
            -subject $genome \\
            -out ${id}.tsv \\
            -outfmt "6 qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq" \\
            -evalue 0.01
        """
    else if (params.types == 'blastp')
        """
        if [ ! -e "${genome}.nin" ] && [ ! -e "${genome}.00.nin" ]; then
                makeblastdb -in $genome -dbtype nucl
        fi
        blastp \\
            -num_threads $task.cpus \\
            -query ${params.query} \\
            -db $genome \\
            -out ${id}.tsv \\
            -outfmt "6 qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq" \\
            -evalue 0.01
        """
    else if (params.types == 'tblastn')
        """
        if [ ! -e "${genome}.nin" ] && [ ! -e "${genome}.00.nin" ]; then
                makeblastdb -in $genome -dbtype nucl
        fi
        tblastn \\
            -num_threads $task.cpus \\
            -query ${params.query} \\
            -db $genome \\
            -out ${id}.tsv \\
            -outfmt "6 qseqid sseqid bitscore evalue pident length mismatch gapopen qstart qend qlen sstart send sstrand slen qseq sseq" \\
            -evalue 0.01
        """
    else if (params.types == 'infernal')
        """
        cmsearch \\
            --cpu $task.cpus \\
            --tblout ${id}.tsv \\
            ${params.query} \\
            $genome
        """
    else
        error "Invalid search type: ${params.types}"
}

// Process to identify neighboring genes
process getNeighbours {
    //publishDir "${params.output_dir}/2_neighbour", mode: 'copy', pattern: '*_{neighbours_output.tsv,neighbours_output.srna.mfna,neighbours_output.protein.mfaa}'
    publishDir "${params.output_dir}/2_neighbour/", mode: 'copy'

    input:
        tuple val(id), 
        path(genome), 
        path(gff), 
        val(search_result)
    
    output:
        tuple val(id), 
        path("${id}.nb.tsv"),           
        path("${id}.nb.protein.mfaa"),  
        path("${id}.nb.srna.mfna")      

    script:
    
    // Convert the Groovy list [gene, ncRNA] into a space-separated string for the shell
    def feature_list = params.including_features.join(' ')

    // Logic: if search_result is a path/file, get its name; otherwise it's just a string
    def hit_input = (search_result instanceof Path) ? "--hit_file $search_result" : ""
    
    """
    mkdir -p "${params.output_dir}/2_neighbour"

    if [ "${params.types}" != "from_gff" ]; then
        python ${projectDir}/bin/get_neighbours_script.py \\
            --hit_file $search_result \\
            --input_type ${params.types} \\
            --fna_file $genome \\
            --gff_file $gff \\
            --gene_of_interest ${params.gene_of_interest} \\
            --neighbours ${params.neighbours} \\
            --output_path ${id}.nb \\
            --promoter ${params.promoter} \\
            --promoter_len ${params.promoter_len} \\
            --including_features ${feature_list} \\
            --evalue_threshold_blast ${params.evalue_threshold_blast} \\
            --evalue_threshold_infernal ${params.evalue_threshold_infernal} \\
            --len_threshold_blast ${params.len_threshold_blast} \\
            --len_threshold_infernal ${params.len_threshold_infernal} \\
            --from_gff_feature ${params.from_gff_feature}
    else
        python ${projectDir}/bin/get_neighbours_script.py \\
            --input_type ${params.types} \\
            --fna_file $genome \\
            --gff_file $gff \\
            --gene_of_interest ${params.gene_of_interest} \\
            --neighbours ${params.neighbours} \\
            --output_path ${id}.nb \\
            --promoter ${params.promoter} \\
            --promoter_len ${params.promoter_len} \\
            --including_features ${feature_list} \\
            --evalue_threshold_blast ${params.evalue_threshold_blast} \\
            --evalue_threshold_infernal ${params.evalue_threshold_infernal} \\
            --len_threshold_blast ${params.len_threshold_blast} \\
            --len_threshold_infernal ${params.len_threshold_infernal} \\
            --from_gff_feature ${params.from_gff_feature}
    fi
    """
}

// Process for clustering and coloring results
process clusterColoring {
    publishDir "${params.output_dir}/3_cluster", mode: 'copy'
    publishDir "${params.output_dir}/3_cluster", mode: 'copy'

    input:
        path(tsv_files)
        path(mfna_files)
        path(mfaa_files)

    output:
        path 'merged_neighbours.tsv'
        path 'merged_protein.mfaa'
        path 'merged_srna.mfna'
        path 'clusterRes.tsv'
        path 'merged_with_color.tsv', emit: colored_tsv
        path "summary.png"
        path "summary.svg"

    script:
    """
    mkdir -p "${params.output_dir}/3_cluster"

    # 1. Merge all neighborhood data
    cat *.nb.tsv | grep -v 'promoter' > merged_neighbours.tsv
    cat *.nb.protein.mfaa > merged_protein.mfaa
    cat *.nb.srna.mfna > merged_srna.mfna

    #mmseqs easy-linclust \\
    #    merged_neighbours.protein.mfaa \\
    #    clusterRes.protein \\
    #    tmp.protein \\
    #    --min-seq-id 0.5 \\
    #    -c 0.8 \\
    #    --cov-mode 1 \\
    #    --cluster-mode 2

    hmmscan \\
        --tblout merged.protein.hmm.tbl \\
        --domtblout merged.protein.hmm.domtbl \\
        --cpu 4 \\
        ${params.pfam_db} \\
        merged_protein.mfaa

    python /home/we93kif/maria_projects/SweetSynteny/bin/postprocess_hmmscan.py \\
        --hmmscan_file merged.protein.hmm.tbl \\
        --output_file clusterRes.hmm.protein_cluster.tsv \\
        --gene_of_interest ${params.gene_of_interest} \\

    cmscan -E 0.01 --cpu 10 --noali --tblout clusterRes.srna \\
        ${params.rfam_db} merged_srna.mfna

    python ${projectDir}/bin/postprocess_cmscan.py \\
        --cmscan_file clusterRes.srna \\
        --output_file clusterRes.srna_cluster.tsv \\
        --gene_of_interest ${params.gene_of_interest}
    
    cat clusterRes.srna_cluster.tsv clusterRes.hmm.protein_cluster.tsv > clusterRes.tsv

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
    path "p_value.png"
    path "p_value.svg"
    path "korrelations_heatmap.png"
    path "korrelations_heatmap.svg"

    script:
    """
    mkdir -p "${params.output_dir}/4_plot/dbscan"
    mkdir -p "${params.output_dir}/4_plot/ward_cosinus"
    mkdir -p "${params.output_dir}/4_plot/ward_euclidean"

    python ${projectDir}/bin/plot_context_script.py \\
        --input_file $merged_with_color \\
        --output_path ./ \\
        --scale ${params.scale} \\
        --cluster ${params.cluster} \\
        --threshold ${params.threshold} \\
        --cut_height_args ${params.cut_height_args} \\
        --gene_of_interest ${params.gene_of_interest} \\
        --name_file ${params.name_file}
    """
}

workflow {
    // Create channel of genome-GFF pairs
    genome_gff_pairs = Channel
        .fromPath("${params.genomes_dir}/*", type: 'dir')
        .map { subfolder -> 
            def fna = subfolder.listFiles().find { it.name.endsWith('.fasta') || it.name.endsWith('.fna') }
            def gff = subfolder.listFiles().find { it.name.endsWith('.gff') }
            if (fna && gff) [subfolder.name, fna, gff]
        }
        .filter { it != null }

    if (params.types == 'from_gff') {
            // Skip search: Prepare the channel to match the input format for getNeighbours
            ready_for_neighbours = genome_gff_pairs.map { id, fna, gff -> 
                tuple(id, fna, gff, params.gene_of_interest) 
            }
    } else {
            // Run search: Process through runSearch and filter
            search_results = runSearch(genome_gff_pairs)
            
            ready_for_neighbours = search_results.filter { id, genome, gff, result ->
                result.size() > 0
            }
    }

    // 3. Final Step: Get neighboring genes
    neighbour_results = getNeighbours(ready_for_neighbours)
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
