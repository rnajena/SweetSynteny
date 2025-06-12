# WORK IN PROGRESS
This pipeline is work-in-progress, you might fing bugs, some are known, while others remain undiscovered. Before getting desperate, please check out the Issues that are already opened and discussed. We encourage the community to contribute by reporting any issues they encounter on GitHub. Feel free to reach out to me via email (maria.schreiber@uni-jena.de) or open an issue directly. It's important to note that I cannot be held responsible for any results obtained using SweetSynteny or any conclusions drawn from them.

***
# SweetSynteny - Unraveling Microsynteny Patterns
Microsynteny, the conservation of gene order and orientation within small genomic regions across different species, provides crucial insights into evolutionary relationships and functional conservation. 

Key features of SweetSynteny:
- Flexible input:
    1. different number of organisms (from bacteria to eukaryotes)
    2. different searches (`cmsearch` for sRNA or `blast` for protein)
- Sequence-driven clustering and color-pattern Microsynteny clustering
    1. on sequence / structur level (-> see Table): `mmseq easy lineclust` or `cmscan`
    2. on microsynteny level: `dendrogram`
    3. on global level: umap comparing microsynteny cluster
- Comprehensive results:
    1. phylogenetic trees using `dendrogram` build by scipy.cluster.hierarchy
    2. statistical summaries
    3. microsynteny plots
    4. statistics on the similarity of the microsynteny locations, e.g. cosinus similarity
    5. Optional: get gene of interest sequence and its promoter sequence (default: 100 nt upstream)

| Conitig:Counter | Gene Name          | Start  | Stop   | Strand| Bio_type       | Color   |
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
    conda create -n nextflow -c bioconda nextflow
    conda activate nextflow
    conda install bioconda::infernal
    conda install bioconda::blast
    conda install bioconda::mmseq
    conda install -c conda-forge matplotlib pandas platformdirs pytest requests seaborn
    ```
3. sugar
   ```
   pip install dna_features_viewer
   pip install rnajena-sugar
   ```
5. Clone the github repository for the latest version of `SweetSynteny`
   ```bash
   nextflow pull rnajena/SweetSynteny
   ```
6.  Done!

## Usage
Let us briefly go over the most important parameters and options. 

<samp>types infernal|blastn|blastp|tblastn </samp>

- For protein(s) we recommended a (m)fasta of amino acid sequences and tblastn
- For sRNA(s) we recommend a corresponding CM from RFAM or self-built\
- You have the choice

<samp>genomes_dir FOLDER </samp> 

- Please choose 2 or more genomes you want to search and save them here.
- And use following structure:
  
    └── genomes_dir
  
        ├── genome1_dir
  
        │    ├── db.gff
  
        │    └── db.fna
  
        ├── genome2_dir
  
        .    ├── db.gff
  
        .    └── db.fna
    ...

<samp>query .cm | .fna </samp>

- Path to CM or FASTA of the gen of interest

<samp>output_dir FOLDER </samp>

- Path to output folder

<samp>gene_of_interest string </samp>

- Name of the gene of interest

<samp>cluster_level sequence_level | structur_level </samp>

- Chose clustering for sRNAs

<samp>neighbours x:y | x-y </samp>

- Set numbers of neighbours (:) or number of nucleotides (-)
- x and y should be Integer numbers

<samp>scale yes | no </samp> 

- Chose if you want to scale the microsynteny plots

<samp>plotting png | svg </samp>

- Select which output format you prefer for the microsynteny plots

<samp>cluster >2 </samp>

- Chose minimal cluster size for `DBscan` clustering

<samp>threshold 0-1 </samp>

- Select a similarity threshold for clustering

### Use a config file.

See example `para.json`

### Running the pipeline
`nextflow run SweetSynteny.nf -params-file /SweetSynteny/para.json`

### Other tools
<details><summary>Click here for all citations</summary>
    
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
