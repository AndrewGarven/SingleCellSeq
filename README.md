# Single Cell Sequencing Data Analysis with Scanpy

## Table of Contents

1. [Description and Acknowledgements](#description)
2. [Libraries](#libraries)
3. [Finding Public Data](#finding-public-data)
4. [Data Download](#data-download)
5. [Quality Control: FastQC](#fastqc)
6. [Read Processing and Alignment: CellRanger](#read-processing-and-alignment-cellranger)
7. [Technical Artifact Removal: CellBender](#technical-artifact-removal-cellbender)
8. [Scanpy Single-Cell Sequencing](#scanpy-single-cell-sequencing)
9. [Scanpy Loading Single-Cell RNA-seq Count Data](#loading-single-cell-rna-seq-count-data-scanpy)
10. [Scanpy Pre-Processing Single-Cell RNA-seq Count Data](#pre-processing-single-cell-rna-seq-count-data-scanpy)
11. [Scrublet Doublet Prediction](#doublet-prediction-scrublet)
12. [Normalization](#normalization)
13. [Dimensional Reduction](#dimensional-reduction)
14. [PCA analysis](#principle-component-analysis)
15. [Harmony Dataset Integration](#harmony-dataset-integration)

## Description and Acknowledgements

The "Single Cell Sequencing Data Analysis with Scanpy" project aims to evaluate publicly available single cell sequencing data from the Sequence Read Archive (SRA) database using the powerful `scanpy` library. This project will encompass the entire data analysis workflow, from data discovery and retrieval to the final evaluation using `scanpy`.

This work is largely based on work done by Dr. Hamid Ghaedi [here](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/README.md)

## Libraries

- [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki)
- [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
- [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ov)
- [Apptainer](https://apptainer.org/docs/user/main/docker_and_oci.html)
- [Docker](https://www.docker.com/)
- [CellBender](https://cellbender.readthedocs.io/en/latest/)
## Finding Public Data

The Sequence Read Archive ([SRA](https://www.ncbi.nlm.nih.gov/sra)) is a publically accessible database maintained by the National Center for Biotechnology Information ([NCBI](https://www.ncbi.nlm.nih.gov)) containing vast archives of bioinformatic sequencing data. As such, it is a great place to start looking for publically available sequencing data that may be useful to answer your research question. In my case, I would like to find a **single-cell sequencing** dataset of **human bladder cancer samples**.  

Thankfully, the SRA database can be searched for 'key terms' using the `esearch` command from the [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) command line tool. Entrez Direct can be downloaded by running the following commands (as outlined [here](https://www.ncbi.nlm.nih.gov/books/NBK179288/)) :

1. `sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"`
2. `echo "export PATH=\$HOME/edirect:\$PATH" >> $HOME/.bash_profile`
3. `export PATH=${HOME}/edirect:${PATH}`

** After Entrez Direct is downloaded approriately, you will need to restart your terminal session

With Entrez Direct downloaded you can run the following command replacing `<Key Terms>` with your own. This function will query the SRA database looking for datasets matching your search criteria and store those results in a .csv file to be viewed in Excel (or perfered editor). 

`esearch -db sra -query '<Key Terms>' | efetch -format runinfo > sra_results.csv`

Note: like all search queries, there are effective modes of setting search criteria or `<Key Terms>` in the SRA ([documentation here](https://www.ncbi.nlm.nih.gov/sra/docs/srasearch/)). For example, to find the single cell sequencing reads from all human bladder cancer samples in the SRA my search may look like the following:

`<Key Terms>` = ((((Homo sapiens[Organism]) AND Single Cell Sequencing)) AND Bladder Cancer)

When I run the above `esearch` command using the above `<Key Terms>` I get the following spread sheet: 

![SRA_results](images/SRA_results.png)

This spread sheet contains an abundance of information regarding all data samples in the SRA that match your search criteria. I would recommend that you sort the spreadsheet by the `"BioProject"` column, and begin manually searching the `"BioProject"` of interest in the SRA database. For example, the `"BioProject"` of most interest to me has the ID `"PRJNA558456"` so I would search the [NCBI](https://www.ncbi.nlm.nih.gov) database as follows:

![BioProjectSearch](images/BioProjectSearch.png)

** Note: please read all the documentation and attached publications to verify the data is as you expect it. If it is we can begin to download all of the samples contained within the selected `"BioProject"`. To do so, I like to download a RunInfo table similar our previously generated `esearch` table above; however, this table will only contain information from our selected `"BioProject"` this can be done by running the following (Replacing `<BioProject>` with your selected Bioproject ID): 

`esearch -db sra -query <BioProject>  | efetch --format runinfo > runinfo.csv`

When I run this script for my BioProject `PRJNA558456` I get the following **runinfo.txt** file: 

![RunInfoBioProject](images/BioProjectRunInfo.png)

## Data Download

As seen [here](https://www.biostars.org/p/359441/#360008) using this runinfo file, we can download SRA data in parallel (GNU) using the [SRA toolkit](https://github.com/ncbi/sra-tools/wiki). 

to download the SRA toolkit please follow the appropriate instructions [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)

Once downloaded you will need to configure your SRA toolkit using the following command 

`vdb-config -i`

this will open an interactive screen you can navigate with your arrow keys and enter.

![vdb-config](images/vdb-config.png)

navigate to `Cache` tab and input a path to an empty directory that will accept all the fastq files from your selected BioProject. In my case, this path is:  `lustre06/project/6065374/garvena/singlecellseq/Data/InputData/Lai`

![vdb-config-path](images/vdb-config-path.png)


Once SRA-toolkit is configured you can 'bulk download' all fastq files from the bioproject of interest using the following command: (please change `<number of samples>` and `<path to runinfo.txt>` to match your data download)


```bash
#!/bin/bash

module load sra-toolkit/3.0.0

parallel --verbose -j 8 curl {} ::: $(tail -n +2 /home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/InputData/test/runinfo.csv | cut -d ',' -f10) >> sra_dump.log
wait
exit
```

## Quality Control: FastQC

Running Quality Control on fastq files downloaded from a public repository is always a good idea. for this I typically use the [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) module. to download fastqc please follow instructions outlined [here](https://raw.githubusercontent.com/s-andrews/FastQC/master/INSTALL.txt)


Once installed you can run the following code changing:
1.  `<path/to/fastq>` with the path to the directory containing the subdirectories with fastq files downloaded above
2. `<path/to/fastQCoutput>` with the path you wish to download the fastqc output
3. `<NumberOfThreads>` with an integer regarding the number of Threads available to multi-thread this process. 


```bash
module load fastqc

fastqc <path/to/fastq>/*/*.fastq.gz -o <path/to/fastQCoutput>/fastQC -f fastq -t <numberOfThreads>
```

fastqc will return 2 files for each sample:

1. ID_fastqc.html
2. ID_fastqc.zip

open the ID_fastqc.html file in your prefered browser by running one of the following commands
1. Google Chrome 

```bash
google-chrome ID_fastqc.html
```
2. Firefox 

```bash
firefox ID_fastqc.html
```

3. Simply locate the file on your device and double click to open in default browser.

This will open up a fastqc report containing a variety of QC metrics for each sample as seen below

![fastqc](images/fastqc.png) 

please navigate through all the summary statistics provided on the right hand side including: Basic Statistics, Per base sequence quality, etc.

## Read Processing and Alignment: CellRanger

Assuming all downloaded fastq files pass quality control checks, we are ready to align and quantify the single-cell sequencing RNA transcripts. This process will be done using [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (Note: CellRanger requires 10X Genomics data as input) 

CellRanger has very strict input requirements, and as such we need to download additional data and further preprocess our data to meet these requirements:

1. [Download CellRanger](#download-cell-ranger)
2. [Rename Fastq Files](#rename-fastq-files)
3. [Download Reference Genome](#download-reference-genome)
4. [Find Chemistry Type](#find-chemistry-type)

### Download Cell Ranger

cellranger can be downloaded by following the instructions provided by 10X genomics [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_in).

in brief:

1. you download the binarys by running

```bash
#!/bin/bash
curl -o cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1690612992&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2OTA2MTI5OTJ9fX1dfQ__&Signature=jqvvJ~lrWzNspxFCdtzBBXz9SIfKS6G31AlqGmj0RMXXoRoPkRfow4StdJkKSQHdillHZX6QEJThl5UaPqiz9dwfVQl~fs2iJA9jVuyvTXzs1VSuboKtNXzESUvbsfEWj3Xt0DykD-V6zeQnbprM397ENlkBvbeUxjRK-HtmMdUCb9K4BlcoXZbB0X7MwrlQNYp1CWxPXT4fIb0Q2aqAPuFGU5DfEyMXzK~hNuOmeQeMye6uF6GPaRHWFN53A6axWsIf55PHmMVQIsLPtGLVkVKMWTOzslvW8mA45lkbVFmlvqZJwfX-zUCWigmnstDUJdKQMBcR27KjxeLbG2aE7Q__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```

unpack those binaries in an appropriate location by running

```bash
#!/bin/bash
tar -xzvf cellranger-7.1.0.tar.gz
```

then prepend your cell ranger directory to your path by running:

*** note change `<path/to/cellranger>` with the location of the binary download above

```bash
#!/bin/bash
export PATH=<path/to/cellranger>/cellranger-7.1.0:$PATH
```

### Rename Fastq Files

As highlighted by [Hamid](https://github.com/hamidghaedi/scRNA_seq-analysis/blob/main/README.md) and a 10X blog post [here](https://kb.10xgenomics.com/hc/en-us/articles/115003802691), natively downloaded SRA files using the fastq-dump have an incompatable file name with cellranger as shown below: 

![cellrangername](images/cellrangername.png)

therefore, we need to rename our `SRR*_1.fastq.gz` file to something like `SRR*_S1_L00*_*_001.fastq.gz` this can be done using the following script. you will need to change `<path/to/input/directory>` to the directory contain the raw fastq files downloaded previously. 

```bash
#!/bin/bash

# Set the path to the parent directory containing the subdirectories (directories with data)
data_dir=<path/to/inputdata/directory>

# Change to the data directory
cd "$data_dir"

# Loop through the subdirectories
for subdir in */; do
    echo "Processing files in $subdir"

    # Change to the subdirectory
    cd "$subdir"

    # Loop through the files and rename them
    for file in *.fastq.gz; do
        if [[ $file == *_S1_L*_R1_* ]]; then
            new_name="${file/_S1_L*_R1_/_R1_}"
        elif [[ $file == *_S1_L*_R2_* ]]; then
            new_name="${file/_S1_L*_R2_/_R2_}"
        else
            echo "Skipping file: $file (Unexpected filename format)"
            continue
        fi

        # Rename the file
        mv "$file" "$new_name"
        echo "Renamed: $file -> $new_name"
    done

    # Change back to the parent directory for the next subdirectory
    cd "$data_dir"
done
```

If executed correctly your data structure should look like the following:

![cellrangerdatastructure](images/cellrangerdatastructure.png)

### Download Reference Genome

As cellranger can be run for a variety of model organisms including: mice (Mus musculus), humans (Homo Sapiens), etc. it is important that we provide cellranger with an appropriate reference genome. All stable reference genomes are provided [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) by 10X genomics, please select the most up-to-date reference. For me (August 2023) I can download the most recent GRCh38 Human reference genome by running the following

```bash 
#!/bin/bash
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```

and unpacked by running the following:

```bash 
#!/bin/bash
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
```

### Find Chemistry Type

This final step requires reading the publication that the data was derived from to evaluate the sequencing chemistry used in the generation of your dataset (often found in the supplemental data). In my case this was found to be `SC3Pv3`; however, this will differ for every dataset (documentation [here](https://kb.10xgenomics.com/hc/en-us/articles/115003764132-How-does-Cell-Ranger-auto-detect-chemistry-)). 

alternatively, you can set this parameter to `auto` and cellranger will attempt to resolve the chemistry for you. 

### Running Cell Ranger

we are now ready to run cellranger on our samples. to do so run the following command inputting the appropriate: 

1. `<path/to/raw/fastq>`
2. `<path/to/ref-genome>`
3. `<chemistry/type>`

*** note this process requires a significant amount of RAM (over 300Gb for my samples) for a long period of time (over 32 hours for my samples)

```bash 
#!/bin/bash
module load cellranger

for d in <path/to/raw/fastq>/*

do
cellranger count --id=$(basename "$d") \
                 --transcriptome=<path/to/ref-genome>/refdata-gex-GRCh38-2020-A \
                 --fastqs=$d \
                 --chemistry=<chemistry/type>
done
```

when cellranger has run correctly on a sample you will receive a message saying:

`Pipestance completed successfully!` 

upon recieveing this message for each sample, you should recieve the output from cellranger in the form of a directory with the same structure as indicated below: 

![cellrangeroutput](images/cellrangeroutput.png)

with the `/outs` directory within each sample directory you should have the files of most interest including directories containing `filtered_feature_bc_matrix` and `raw_feature_bc_matrix`.

Inside the `filtered_feature_bc_matrix` and `raw_feature_bc_matrix` directories, you may find the following files:

1. `matrix.mtx`: This is the main file containing the digital gene-barcode matrix in a sparse format. It represents the number of UMIs (unique molecular identifiers) detected for each gene in each cell.
2. `features.tsv`: This file provides information about the gene features (genes) included in the matrix. It typically includes the gene names and other relevant metadata.
3. `barcodes.tsv`: This file contains information about the cell barcodes used to uniquely identify each cell in the matrix.

These files are what are going to be leveraged for Scanpy Analysis.

*** note: there are a variety of useful outputs for a whole number of things that may be covered at a later date. for now, feel free to evaluate these outputs on your own (documentation [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct))

## Technical Artifact Removal: CellBender

Single-cell sequencing assays produce a significant amount of background noise counts primarily derived from cell-free droplets containing RNA fragments. This can result in a variety of errors including batch effects, incorrect differential expression, etc. As a result, we need to remove this background 'noise' from each of our samples. This can be done using [CellBender](https://github.com/broadinstitute/CellBender) ([publication](https://www.biorxiv.org/content/10.1101/791699v2), [documention](https://cellbender.readthedocs.io/en/latest/)) CellBenders main functionality `remove-background` is a module that removes counts due to ambient RNA molecules and random barcode swapping from (raw) UMI-based scRNA-seq gene-by-cell count matrices.

1. [Download CellBender](#download-cellbender)
2. [Run remove-background](#run-remove-background)

### Download CellBender
CellBender is most easily downloaded as a docker image as this will automatically download all dependancies. Thankfully, CellBender is available from the Google Container Repository. This can be downloaded via any process capable of pulling docker images; for this I will use apptainer:

***note: this is pulling the latest cellbender docker image (August 2023) and therefore will likely in the future

```bash
module load apptainer
apptainer pull cellbender.sif docker://us.gcr.io/broad-dsde-methods/cellbender:latest 
```

### Run remove-background

Remove background attempts to remove RNA reads stemming from the 'ambient plateau'. this plateau simply refers to the section of the counts-per-droplets VS. Droplet ID ranked by count plot (examples below) where the 'cells' seem to exibit significantly less RNA reads. This isnt always a low amount (see exibit B: high background); however, it should be a signifcant decrease in expression campared to other cells in the plot. this section of the plot represents what is not likely to be a cell, but instead a droplet containing a significant amount of RNA fragments. Again, if the top left portion of this plot represents a probability of being a cell, then the 'plateau' section of this graph is highly unlikely to be a true cell count as its expression is too greatly reduced and therefore likely represents an artifact. 

![ambientplateau](images/ambientplateau.png)


![cellbenderhigh-low](images/cellbenderhigh-low.png)

properly interpreting these graphs is extremely important for maximising the utility of cellbender as the `remove-background` function requires setting an `expected-cells` and `total-droplets included` values. 

Thankfully `cellranger` provides a UMI-curve as an output from within the `web_summary.html` output file. This file provides a bunch of useful information including an estimated cell count, mean reads per cell. For our purposes we simply need `Estimated Number of Cells` and a rough estimate of the middle of the `ambient plateau` these values will represent `expected-cells` and `total-droplets included values` for cellbender. For example, in the sample provided below I would set these values to ~8k and ~15K, repectively. please note that we will be running cellbender in a for loop, so we will need a rough average values for all your samples (round up)!

![web_summary](images/web_sum.png)

to run cellbender from the docker image you will need to run `module load apptainer` then make use of the exec function to execute cellbender. please adjust the `<path/to/cellranger/*.h5>` and `<path/to/cellbender/output.h5>` to the appropriate directories. In addition, please adjust `expected-cells` and `total-droplets-included` as described above. please note: the `--cuda` tag allows cellbender to utilize a GPU which signifcantly reduces the time it takes to run this opperation. if you do not have access to a GPU, please remove this tag. 

***Note: you must use the apptainer --bind tag to access filepaths outside of the docker image (this would include all data for processing)

```bash
module load apptainer 

# Set the root directory where the subdirectories with raw_feature_bc_matrix.h5 files are located
input_directory="<path/to/input/data/directory>"
output_directory="<path/to/desired/output/directory>"

# Iterate through each subdirectory and process the raw_feature_bc_matrix.h5 file
for subdir in "$input_directory"/*; do
  # Check if the raw_feature_bc_matrix.h5 file exists in the current subdirectory
  if [ -f "$subdir/outs/raw_feature_bc_matrix.h5" ]; then
    sample=$(basename "$subdir")
    
    # Process the raw_feature_bc_matrix.h5 file in the current subdirectory
    apptainer exec -C -W ${SLURM_TMPDIR} --nv \
      --bind "$input_directory":"/mnt/input" \
      --bind "$output_directory":"/mnt/output" \
      cellbender.sif cellbender remove-background \
      --input "/mnt/input/$(basename "$subdir")/outs/raw_feature_bc_matrix.h5" \
      --output "/mnt/output/$sample.h5" \
      --cuda \
      --expected-cells 8000 \
      --total-droplets-included 20000 \
      --fpr 0.01 \
      --epochs 150
  fi
done
```

CellBender has a variety of inputs for each sample including '<sample_name>_filtered.h5' which will be utilized in downstream analyses. Before proceeding please investigate the '<sample_name>.pdf' document as this contains much of the quality control outputs. Most important of which is the UMI plot from above with the overlayed prediction made by the CellBinder software. Please variefy that the cells removed from your downstream analysis (show a low probability of being a cell) seems reasonable to you. An example from my dataset is shown below. 

![cellbenderpredict](images/cellbenderpredict.png)

## Scanpy Single-Cell Sequencing

### Loading Single-Cell RNA-seq Count Data: Scanpy

[Scanpy](https://github.com/scverse/scanpy) (as noted by there team) is a scalable toolkit for analyzing single-cell gene expression data. It includes methods for preprocessing, visualization, clustering, pseudotime and trajectory inference, differential expression testing, and simulation of gene regulatory networks. Its Python-based implementation efficiently deals with data sets of more than one million cells. 

As this data is highly interactive, requiring multiple intermediate steps I highly recommend use of a jupyter-notebook (a working example of which can be found at */PythonScripts/main.ipynb)

This may not accurately represent all required dependancies (will update as I go)
```python 
import scanpy as sc
import anndata as ad
from glob import glob
import matplotlib.pyplot as plt
```

Scanpy has weird default values that are often immediately changed by convention. The first is that scanpy by default will only show warning+error messages. We would like to change this as it also offers 'Hints' if you set verbosity to 3.

```python
sc.settings.verbosity = 3
```
Here we are simply iterating through all of the '*/outs/filtered_feature_bc_matrix' directories and importing the single cell dataset into a scanpy object. 

1. I generate a list of directory paths to cellranger output directories
2. I generate a dictionary containing anndata objects with sample_name as  `keys`

*** setting cache=True will significantly speed up reloading these objects a second time and will generate additional `.h5ad` files in the directory where this code is executed. 

```python
dirpaths = glob('path/to/cellbender/h5/files/*_filtered.h5')

adatas = {dirpath.split('/')[10].split('_')[0]: sc.read_h5ad(dirpath, cache=True) for dirpath in dirpaths}
adata = ad.concat(adatas.values(), label='sample')
```

### Pre-Processing Single-Cell RNA-seq Count Data: Scanpy

making an integrated anndata object containing all samples

```python
adata = ad.concat(adatas, label='sample')
adata.obs_names_make_unique
adata.var_names_make_unique
```

A common first pre-processing step is to evaluate top gene expression in your samples, to ensure they align with what you would expect within your samples. It is often difficult to see these distributions because they are heavily right skewed so I like to include the `log=True` tag to show log transformed data. please note the x-axis is log transformed (as displayed below).

```python
sc.pl.highest_expr_genes(adata, n_top=20, show=True, log=True)
```

![preprocess-topgene](images/preprocess-top-gene.png)

Additionally, we should consider the distribution of genes per cell, and the number of cells that express a specific gene. This will be done for each sample within our dataset individually. 

```python

for adata in adatas.values():
    sc.pp.filter_cells(adata, min_genes=0)
    gene_counts_per_cell = adata.obs['n_genes']
    sc.pp.filter_genes(adata, min_cells=0)
    cell_counts_per_gene = adata.var['n_cells']
    
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.hist(gene_counts_per_cell, bins=50)
    plt.xlabel('Number of Genes per Cell')
    plt.ylabel('Frequency')
    
    plt.subplot(1, 2, 2)
    plt.hist(cell_counts_per_gene, bins=50)
    plt.xlabel('Number of Cells per Gene')
    plt.ylabel('Frequency')
    
    plt.tight_layout()
    plt.show()

```
This should output a distribution for each sample like so: 

![gene-cell-distribution](images/genecelldis.png)

You will notice that a significant proportion of the cells captured in your sample, have very few genes. What are these? Do these populations represent a 'true' cell within your sample? likely not! These 'cells' are more likely to be captured mitochondria, ribosomes, or hemaglobin which are likely not of interest. As such, we should assess the mitochondrial, ribosomal, and hemoglobin based gene signatures in each of our computed cells, and eliminate all 'cells' that seem to express a significant fraction of these genes. To do so, we can leverage gene naming convention, and label cells with extremely high rates of transcrition of these genes. Mitochondrial genes all start with `'MT'`, Ribosomal genes all start with `'RPS'` or `'RPL'` and Hemoglobin genes all start with `'^HB[^(P)]'`.

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
```
Now that we have isolated these gene sets, we can compute quality control metrics using [calculate_qc_metrics](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.calculate_qc_metrics.html) like so: 

```python
sc.pp.calculate_qc_metrics(adata,
                           qc_vars=["mt", "ribo", "hb"],
                           inplace=True,
                           log1p=True)

```
Lets now generate violin plots for these QC metrics to assess their prominence in our dataset. I am only displaying the results for `'pct_counts_mt'`; however, please note this can be done for all QC metrics by replacing `'pct_counts_mt'` with `'pct_counts_ribo'` or `'pct_counts_hb'`

```python
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

```

![violin-QC-mt](images/violinQCmt.png)

Additionally, we can produce a scatterplot assessing the disribution of mt, ribo and hb genes in our cells. The same is true here - you can replace `'pct_counts_mt'` with `'pct_counts_ribo'` or `'pct_counts_hb'` 

```python
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
```
![scatter-QC-mt](images/scatterQCmt.png)

Now we have grounds to remove cells that highly express any of these genes. Again, focusing on the scatterplot, we can to eliminate all cells with a high expression of mt, ribo, or hb genes. However, some cells that highly express these genes are likely to represent real biology! As such, scanpy highly recommends that you start with a permissive filtering stratergy. As such, based on my QC plot, I plan to eliminate all cells with less than 100 genes, and all genes that are expressed in less than 5 cells. I can do so by utilizing `'filter_cells'` and `'filter_genes'` like so. 

```python
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=5)
```

### Doublet Prediction: scrublet

Another common preprocessing step for single-cell sequencing analysis is removing doublets. Doublets are simply multiple cells that were co-encapsulated in a single bead when preparing the sample for sequencing. These 'cells' therefore must be removed from our analysis. To do so we will utilize the python package [scrublet](https://github.com/swolock/scrublet). Unfortunately, scrublet is fairly computational expensive, as such I must utilize a compute cluster from 'compute canada'. To do so, I must generate a shell script that generates a python enviroment and then runs scrublet on our samples. Thankfully, the scanpy package has a deployable version of scrublet for our purposes. Therefore, we only need to import scanpy for this process. 

```bash
#!/bin/bash
#SBATCH --account=account-name
#SBATCH --job-name=scrublet
#SBATCH --qos=privileged
#SBATCH --nodes=1                
#SBATCH --tasks-per-node=1         
#SBATCH --mem 50G
#SBATCH --time 10:00:00
#SBATCH --output=path/to/output/file/scrublet.%J.out
#SBATCH --error=path/to/err/file/scrublet.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#

module load StdEnv/2023

#generate python enviroment
module load python/3.10
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip

#load scanpy 
pip install --no-index scanpy

# Run your Python script
python path/to/scrublet/file/slurm_scrublet.py
```
To run scrublet through the scanpy package, we can simply call [sc.external.pp.scrublet](https://scanpy.readthedocs.io/en/1.9.x/generated/scanpy.external.pp.scrublet.html) like so: 

```python
import scanpy as sc
import glob

outpath = '/path/to/scrublet/output/dir'

dirpaths = glob('path/to/h5/files/*.h5')

for dirpath in dirpaths:
    id = dirpath.split('/')[10].split('_')[0]
    adata = sc.read_10x_h5(dirpath) 
    sc.external.pp.scrublet(adata, verbose=True)
    sc.write(f"{outpath}/{id}_scrub.h5", adata)

```
This will generate a new H5 file with all doublets removed from your sample. You will now need to re-import these files into your working enviroment. This can be done as before: 

```python
dirpaths = glob('path/to/scrublet/h5/files/*_scrub.h5')

adatas = {dirpath.split('/')[10].split('_')[0]: sc.read_h5ad(dirpath, cache=True) for dirpath in dirpaths}
adata = ad.concat(adatas.values(), label='sample')
```

### Normalization

Now that we have a reasonable clean single-cell dataset, we are ready for data normalization. For this, we will normalize the count depth of our gene expression and perform a log transformation. first we must save a copy of our unnormalized gene expression by:

```python
adata.layers["counts"] = adata.X.copy()
```
now we can normalize our counts by the median total counts per cell. This is done through the [sc.pp.normalize_total](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html) function in scanpy. Additionally we can log tranform our data using the [sc.pp.log1p](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.log1p.html) function. 

```python
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

### Dimensional Reduction

To reduce the diminsionality of our dataset, a logical first step is to assess the utility of our features, a process called feature selection. For example, a substatial fraction of the genes that will be expressed in your dataset will be uniformly distributed in all cells. These genes, commonly called 'house keeping genes', do not offer any utility in deliniating cell-type, computing trajectories, etc. As such, we are first going to reduce the diminsionality of our dataset by removing uniformly distributed genes. To do so, we will select the most highly variable genes by using the [sc.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html) function. This function computes the dispertion of each genes within your single cell dataset, and returns the top `n` genes withinn your dataset through the `n_top_genes` parameter. We can then plot the selected genes (based on the set cutoff) using sc.pl.highly_variable_genes function.

```python
sc.pp.highly_variable_genes(adata,
                            n_top_genes=2500,
                            batch_key="sample")
sc.pl.highly_variable_genes(adata)
```
![highly-variable-genes](images/highlyvariablegenes.png)

#### Principle Component Analysis

In order to further reduce the diminsionality of our dataset, we should leverage PCA. PCA computes the main axes of variation within our dataset. If you are unfamiliar with PCA analysis, I strongly encourage you to further investigate this proceedure. Here are some great resources:

1. [Greenacre et, al.](https://www.nature.com/articles/s43586-022-00184-w#:~:text=Principal%20component%20analysis%20is%20a,variance%20of%20all%20the%20variables.)
2. [StatQuest](https://youtu.be/FgakZw6K1QQ?si=iyUOR2jcbUEzC-ah)
3. [IBM explaination of PCA](https://www.ibm.com/topics/principal-component-analysis)

To perform PCA on our dataset, we will use the [sc.tl.pca](https://scanpy.readthedocs.io/en/1.9.x/generated/scanpy.tl.pca.html) function like so:

```python
sc.tl.pca(adata)
```

With our principle components computed the next step is the assess the relative influence each PC has on the total variance in our dataset. This can be done through use of the [sc.pl.pca_variance_ratio](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.pca_variance_ratio.html) function in scanpy.

```python
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
```
![explain-variance-ratio](images/explainedvarianceratio.png)

It is often stated that the most optimal cut-off for PCA is the 'elbow' of the plot. This is the point at which the addition of an additional PC adds very little varience to our dataset. In my case, this would likely be around PC 18 - 24. However, as the addition of PC 25 - 50 do not heavily negatively impact results (stated by scanpy) there really is very little harm in including all 50 PCs. 

As a quick sanity check we should ensure that our quality control metrics (things like the distribution of MT genes) do not explain the varience in our top PCs. As such, we can plot our MT fraction of Pairs-Plots of the top 8 PCs. 

```python
sc.pl.pca(
    adata,
    color=['sample', 'sample', 'sample', 'sample',
           'pct_counts_mt', 'pct_counts_mt', 'pct_counts_mt', 'pct_counts_mt'],
    dimensions=[(0, 1), (2, 3), (4, 5), (6, 7),
               (0, 1), (2, 3), (4, 5), (6, 7)],
    ncols=4,
    size=4,
)
```

![PCA-mt](images/PCA-mt.png)

As shown above, there a no dedicated seperate clusters that are fully explained by MT fraction. As such, our QC seems to have worked well! However, in plot 1 (comparing PC1 vs PC2) there seems to be a large degree of seperation between samples (i.e. the varience within our dataset can be partially explained by a 'batch effect'). We will explore this further below!

To assess the degree of 'batch effects' between samples, we will quickly compute the [Nearest Neighbours](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html) within our PCA representation matrix, and plot the output in 2-dimensions using [UMAP](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html). Both of these functions will be outlined in greater detail after we correct for batch effects (Harmony).

```python
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata,
           color="sample",
           size=2)
```

![UMAP-batcheffect](images/umap-batch.png)

The above Uniform Manifold Approximation and Projection plot is labeling our samples (colours) and clustering the cells by there transcriptional expression. As you may note in the above image, our samples seem to be clustering independantly! This is a clear sign of 'batch effects'. To correct for this technical error, we will utilize [Harmony](https://github.com/lilab-bcb/harmony-pytorch) to correct for these effects. 

### Harmony Dataset Integration

Harmony is an algorithm developed by [Ilia Korsunsky et, al.](https://www.nature.com/articles/s41592-019-0619-0) that projects cells into a shared embedding in which cells group by cell type rather than dataset-specific conditions. As such, Harmony is an excellent algorithm for removing the technical errors posed by batch effects. This algorithm can also be used to integrate datasets from multiple sources in increase the power of your analysis. 

Unfortunately, As with Scrublet, Harmony is very computationally expensive. As such, as before, I must use compute canada for such analyses. 

to run this script on SLURM you need to send your shell script using sbatch like so

```bash
sbatch path/to/harmony/shell/script/run_harmony.sh
```

SLURM shell script.
```bash
#!/bin/bash
#SBATCH --account=account-name
#SBATCH --job-name=harmony
#SBATCH --qos=privileged
#SBATCH --nodes=1                # number of Nodes
#SBATCH --tasks-per-node=1        # number of MPI processes per node
#SBATCH --mem 10G
#SBATCH --time 00:30:00
#SBATCH --output=path/to/output/file/harmony.%J.out
#SBATCH --error=path/to/error/file/harmony.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#

module load StdEnv/2023

module load python/3.11
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip

pip install --no-index harmony_pytorch
pip install --no-index scanpy

python path/to/harmony/python/file/slurm_harmony.py
```
SLURM python script 

```python
import scanpy as sc
from harmony import harmonize


adata = sc.read_h5ad('/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/Harmony/Lai/Input/input_harmony.h5ad')
Z = harmonize(adata.obsm['X_pca'], adata.obs, batch_key = 'sample')
adata.obsm['X_harmony'] = Z

sc.write("/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/Harmony/Lai/output/output_harmony.h5ad", adata)
```
