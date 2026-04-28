# WORK IN PROGRESS
This pipeline is work-in-progress, you might find bugs, some are known, while others remain undiscovered. Before getting desperate, please check out the Issues that are already opened and discussed. We encourage the community to contribute by reporting any issues they encounter on GitHub. Feel free to reach out to me via email (maria.schreiber@uni-jena.de) or open an issue directly. It's important to note that I cannot be held responsible for any results obtained using SweetSynteny or any conclusions drawn from them.

***
# SweetSynteny - Unraveling Microsynteny Patterns
SweetSynteny identifies microsynteny patterns (conservation of gene order and orientation within small genomic regions) across genomes by combining sequence-based clustering, gene-context analysis, and visualization. It supports both protein and non-coding RNA searches and produces phylogenetic, cluster, and neighborhood summaries.

Microsynteny provides crucial insights into evolutionary relationships and functional conservation. 

Key features of SweetSynteny:
- Flexible input:
    1. different number of organisms (from bacteria to eukaryotes)
    2. different searches (`cmsearch` for sRNA or `blast` for protein or with the gff files `from_gff`)
- Data filtering (E-value, hit length)
- Sequence-driven clustering and color-pattern Microsynteny clustering
    1. on sequence / structur level: `mmseq easy lineclust` or `cmscan` or `hmmscan` (see table)
    3. pca for dimension reduction
    2. on microsynteny level: on global level: microsynteny cluster by ward or dbscan 
- Comprehensive results:
    1. phylogenetic trees using `dendrogram` build by scipy.cluster.hierarchy or scatterplot
    2. statistical summaries of adjacent genes and genome location
    3. microsynteny plots
    4. statistics on the similarity of the microsynteny locations, e.g. cosinus similarity
    5. Optional: get gene of interest sequence and its promoter sequence (default: 100 nt upstream or up to the next adjacent gene)
- Implementation: Nextflow

| Contig:Counter | Gene Name          | Start  | Stop   | Strand| Bio_type       | Color   |
|-----------------|--------------------|--------|--------|-------|----------------|---------|
| NZ_CP013002.1:0 | gene-AQ619_RS00960 | 215167 | 216307 | sense | protein_coding | #FFFFFF |

So, as you can see, with SweetSynteny, your Microsynteny analysis will be, well... sweet!

***

## Graphical Workflow

![Workflow graph](/fig/workflow.2.png)

## Dependencies and installation

The pipeline is written in Nextflow. In order to run `SweetSynteny`, I recommend creating a conda environment dedicated for NextFlow.
1. Install [miniconda](https://docs.conda.io/projects/miniconda/en/latest/) or [conda]()
2. Create a conda environment and install NextFlow within this environment and install everything else.
    ```bash
    mamba create -n env_name
    conda activate env_name
    mamba install -c conda-forge -c bioconda   nextflow openjdk   \
        infernal blast mmseqs2   \
        matplotlib pandas platformdirs pytest requests seaborn numpy scipy scikit-learn
    ```
3. sugar
   ```
   pip install rnajena-sugar
   ```
5. Clone the github repository for the latest version of `SweetSynteny`
   ```bash
   nextflow pull rnajena/SweetSynteny
   ```
6.  Done!

## Input layout

Please choose 2 or more genomes you want to search and save them here. Prepare a genome directory containing one subdirectory per genome. Each genome folder should include a GFF annotation file and a FASTA genome file.
- And use following structure:
  
    └── genomes_dir
  
        ├── genome1_dir
  
        │    ├── db.gff
  
        │    └── db.fna
  
        ├── genome2_dir
  
        .    ├── db.gff
  
        .    └── db.fna
    ...

## Usage

Let us briefly go over the most important parameters and options. 

```bash
nextflow run SweetSynteny.nf -params-file /path/to/para.json -c nextflow.config
```

<samp>search_types infernal|blastn|blastp|tblastn </samp>

- For protein(s) we recommended a (m)fasta of amino acid sequences and tblastn
- For sRNA(s) we recommend a corresponding CM from RFAM or self-built\
- You have the choice

<samp>bio_type ncRNA|protein </samp>

<samp>genomes_dir FOLDER </samp> 

<samp>annotation_type .gff | other_types_for_futur </samp>

<samp>query .cm | .fna </samp>

- Path to CM or FASTA of the gene of interest

<samp>output_dir FOLDER </samp>

- Path to output folder

<samp>gene_of_interest string </samp>

- Name of the gene of interest

<samp>adjacent_gene_clustering hmmscan,cmscan | mmseq,cmscan | hmmscan,mmseq | mmseq,mmseq </samp>

- Choose clustering for adjacent genes

<samp>neighbours x:y | x-y </samp>

- Set numbers of neighbours ( : ) or number of nucleotides ( - )
- x and y should be Integer numbers
- It is also possible for e.g. riboswitches to write 0,4 and only focus on the downstream genes.

<samp>scale yes | no </samp> 

- Choose if you want to scaled and aligned the microsynteny plots

<samp>cluster >2 </samp>

- Choose minimal cluster size for `DBscan` clustering

<samp>threshold 0-1 </samp>

- Select a similarity threshold for clustering

<samp> cut_height_args float </samp>

- Cutting threshold for ward clustering
 
<samp>pfam_db : /path/to/result/folder/Pfam-A.hmm</samp>

- please, download the pfam db and call ... 

<samp>rfam_db : /path/to/result/folder/Rfam.cm</samp>

- please, download the rfam db and call ... 

<samp>name_file : ""</samp>

- Path to genome name file
- It should look like this:

| strain          | contig               | organism_name                                      |
|-----------------|----------------------|----------------------------------------------------|
| GCF_000731315.1 | NZ_HG938354.1        | Neorhizobium galegae bv. orientalis str. HAMBI 540 |
| GCF_000731315.1 | NZ_HG938353.1        | Neorhizobium galegae bv. orientalis str. HAMBI 540 |
| GCF_042657465.1 | NZ_JBHSLC010000080.1 | Azospirillum himalayense                           |
| GCF_042657465.1 | NZ_JBHSLC010000008.1 | Azospirillum himalayense                           |
| GCF_042657465.1 | NZ_JBHSLC010000094.1 | Azospirillum himalayense                           |

<samp>ignore_overlaps : "True"|"False"</samp>
- When you know you search hit overlaps with another annotation, but not more than 75%

<samp>substring_search : "True"|"False"</samp>
- Only when `from_gff`
- If you want to search for all SRPs but in you gff file you find SRP, bacterial_SRP, etc

<samp>cpu : int</samp>

### Use a config file.

See example `para.json`

### Output interpretation

SweetSynteny generates:

- cluster and microsynteny summary tables
- gene neighborhood tables
- phylogenetic trees and scatterplots
- statistical summaries of adjacent genes and genome locations
- `summary.png`, `summary.svg`, and `summary.json`

### Other tools
<details><summary>Click here for all citations</summary>

  * SUGAR:
    * `Eulenfeld, Tom. "Sugar: A Python framework for bioinformatics." Journal of Open Source Software 10.111 (2025): 8122.`

  * BLAST:
    * `Korf, Ian, Mark Yandell, and Joseph Bedell. Blast. " O'Reilly Media, Inc.", 2003.`
      
  * INFERNAL:
     * `Nawrocki, Eric P., Diana L. Kolbe, and Sean R. Eddy. "Infernal 1.0: inference of RNA alignments." Bioinformatics 25.10 (2009): 1335-1337.`
       
  * MMSeqs2:
    * `Steinegger, M., Söding, J. "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets". Nat Biotechnol 35, 1026–1028 (2017)`
      
  * ETE3:
    * `Huerta-Cepas, Jaime, François Serra, and Peer Bork. "ETE 3: reconstruction, analysis, and visualization of phylogenomic data." Molecular biology and evolution 33.6 (2016): 1635-1638.`

  * DNA Features Viewer
    * `Edinburgh Genome Foundry by Zulko. https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer`      
</details>

## Cite us
If you use SweetSynteny for your analysis, please cite our github repository.

```bibtex
@software{Maria_Schreiber_SweetSynteny,
author = {Maria Schreiber, Emanuel Barth, Manja Marz},
license = {MIT},
title = {{SweetSynteny}},
url = {https://github.com/rnajena/SweetSynteny}
}
```
