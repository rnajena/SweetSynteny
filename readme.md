# WORK IN PROGRESS - DISCLAIMER
This pipeline is work-in-progress, you might fing bugs, some are known, while others remain undiscovered. Before getting desperate, please check out the Issues that are already opened and discussed. We encourage the community to contribute by reporting any issues they encounter on GitHub. Feel free to reach out to me via email or open an issue directly. It's important to note that I cannot be held responsible for any results obtained using SweetSynteny or any conclusions drawn from them.

# TODO
- bin folder for python scripts
- modules folder for nf

***
# SweetSynteny - Overview
- Searching with `blastn` or `tblastn` or `cmsearch`
- Gets neighours of you hit and saves them in a tsv file
- Clustering with 
    - on sequence / structur level: `mmseq easy lineclust` or `cmscan` [TODO] -> see Table
    - on microsynteny level: `DBscan` or `h`

| Conitig:Counter | Gene Name          | Start  | Stop   | Ori   | Bio_type       | Color   |
|-----------------|--------------------|--------|--------|-------|----------------|---------|
| NZ_CP013002.1:0 | gene-AQ619_RS00960 | 215167 | 216307 | sense | protein_coding | #FFFFFF |

- Generates micorsynteny plots
- Compares ...
- Generates trees using `ete3`
***

## Graphical Workflow

![Workflow graph](/fig/workflow.png)

## Dependencies and installation
The pipeline is written in Nextflow. In order to run `SweetSynteny`, I recommend creating a conda environment dedicated for NextFlow.
1. Install [miniconda](https://docs.conda.io/projects/miniconda/en/latest/) or [conda]()
2. Next, make sure that conda is part of your $PATH variable, which is usually the case.
3. Create a conda environment and install NextFlow within this environment and install everything else.
    ```bash
    conda create -n nextflow -c bioconda nextflow
    conda activate nextflow
    conda install bioconda::infernal
    conda install bioconda::blast
    conda install bioconda::mmseq
    conda install -c conda-forge matplotlib pandas platformdirs pytest requests seaborn
    pip install rnajena-sugar
    ```
4. Clone the github repository for the latest version of `SweetSynteny`
    ```bash
    nextflow pull rnajena/SweetSynteny
    ```
5.  Done!

## Usage
Let us briefly go over the most important parameters and options. 
--hit_file
--input_type
--gff_file
--fna_file
--gene_of_interest
--neighbours
--output_path

### Use a config file.

### Running the pipeline
`nextflow run SweetSynteny.nf -params-file /home/we93kif/maria_projects/SweetSynteny/para.json`


### Other tools
<details><summary>Click here for all citations</summary>
  * BLAST
  * INFERNAL
  * MMSeqs2:
    * `Steinegger, M., Söding, J. "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets". Nat Biotechnol 35, 1026–1028 (2017)`
  * DBscan
  * ETE3
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




"""
python /home/we93kif/maria_projects/SweetSynteny/bin/plot_context_script.py \
--scale no \
--input_file /data/fass5/projects/ProjectFroehlich/data/result/merged_and_clustered_results_crfa5_proteinVSsrna/merged_with_color.tsv \
--output_path /home/we93kif/maria_projects/SweetSynteny/test \
--output_ending png
"""