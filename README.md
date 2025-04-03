# Analysis of miRNA targets and signaling pathways involved using multiMiR

A simple script to demonstrate the capabilities of R for searching miRNA targets across multiple databases, and performing functional enrichment analysis to identify associated signaling pathways.

## Usage

### 1. Download the repo from GitHub

Use the git clone command to create a local copy of the repository.

```shell
git clone https://github.com/villena-francis/mirnas_exercise
```

### 2. Set up the environment

This pipeline requires a [conda package manager](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). We recommend using [Miniforge](https://github.com/conda-forge/miniforge), a lightweight conda alternative that prioritizes compatibility with conda-forge channel.

To generate an environment with R and all the necessary packages to run the script use:

```shell
# Navigate to the visor-simulations directory
cd /path/to/mirnas_exercise
# Create the environment
conda env create -f env/r_mirnas.yaml
```

You will also need to add a list of miRNAs on which you want to run the analysis, in the `data` directory as a text file named `mirnas.txt`. Additionally, you will need a list of key terms to refine the pathways in relation to your topic of interest, in the same directory as a file named `context_keywords.txt` (you can use the example files to check the format and test the script).



### 3. Run the script

Make sure you are still in `mirnas_exercise` directory, activate the 
Conda environment **r_mirnas**, and execute the R script:

```shell
# Activate the environment
conda activate r_mirnas
# Run the script
Rscript mirnas.R
```

When you are done, your working directory should look something like this:

```
.
├── data
│   ├── context_keywords.txt
│   └── mirnas.txt
├── env
│   └── r_mirnas.yaml
├── mirnas.R
├── README.md
└── results
    ├── all_mirnas_complete_targets.csv
    ├── analysis_summary.txt
    ├── athero_relevant_terms.txt
    ├── common_targets.txt
    ├── GO_athero_barplot.pdf
    ├── GO_athero_network.pdf
    ├── GO_enrichment_results.csv
    ├── KEGG_barplot.pdf
    ├── KEGG_network.pdf
    └── miRNA_targets_venn.pdf
```
