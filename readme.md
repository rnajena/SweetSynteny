# WORK IN PROGRESS - DISCLAIMER
This pipeline is work-in-progress, you might fing bugs, some are known, while others remain undiscovered. Before getting desperate, please check out the Issues that are already opened and discussed. We encourage the community to contribute by reporting any issues they encounter on GitHub. Feel free to reach out to me via email or open an issue directly. It's important to note that I cannot be held responsible for any results obtained using SweetSynteny or any conclusions drawn from them.

***
# SweetSynteny - Overview
- Searching with `blastn` or `tblastn` or `cmsearch`
- Gets neighours of you hit and saves them in a tsv file
- Clustering with 
    - on sequence / structur level (-> see Table): `mmseq easy lineclust` or `cmscan` [TODO] 
    - on microsynteny level: `DBscan` or `Hierarchical clustering` [TODO]

| Conitig:Counter | Gene Name          | Start  | Stop   | Strand| Bio_type       | Color   |
|-----------------|--------------------|--------|--------|-------|----------------|---------|
| NZ_CP013002.1:0 | gene-AQ619_RS00960 | 215167 | 216307 | sense | protein_coding | #FFFFFF |

- Generates micorsynteny plots and gives you statistics on the similarity of the microsynteny locations
- Generates trees using `ete3` [TODO]
***

## Graphical Workflow

![Workflow graph](/fig/workflow.png)

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

<samp>genomes_dir FOLDER </samp>      

<samp>query .cm | .fna </samp>

<samp>output_dir FOLDER </samp>

<samp>gene_of_interest string </samp>

<samp>cluster_level sequence_level | sequence_level </samp>

<samp>neighbours x:y | x-y </samp>

<samp>scale yes | no </samp> 

<samp>plotting png | svg </samp>

<samp>cluster >2 </samp>

<samp>threshold 0-1 </samp>

### Use a config file.

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
