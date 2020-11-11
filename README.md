# Seq2FunR
R package for high-performance shotgun metagenomics data processing



### Note - Seq2FunR is still under development - we cannot guarantee full functionality ### 

## Description 

**Seq2FunR** contains the R functions and libraries to conduct end-to-end analysis of shotgun microbiome sequencing data for taxon classification. It starts with raw reads and conduct quality check of these reads (removing low quality reads, too short reads, trimming low quality score bases, correcting sequencing error for overlapped region of  overlapped paired end reads and joining overlapped paired end reads), followed by translated search against a protein database to identify their homology protein sequences. A taxon abundance table is generated by mapping protein sequences to taxon information. Here, we use a subset of clinical IBD samples as a show case. Seq2FunR can be used as a part of **MicrobiomeAnalystR**. 

## Getting Started 

### Step 1. Install package dependencies 

To use Seq2FunR , first install all package dependencies. Ensure that you are able to download packages from bioconductor. To install package dependencies, use the pacman R package (for those with >R 3.5.1). Note that some of these packages may require additional library dependencies that need to be installed prior to their own successful installation. 

```R
install.packages("pacman")

library(pacman)

pacman::p_load(crayon, magrittr, stringr);
```
### Step 2. Install the package 

Seq2FunR is freely available from GitHub. Note that the Rpackage is currently under construction. The package documentation, including the vignettes for each module and user manual is available within the downloaded R package file. If all package dependencies were installed, you will be able to install the Seq2FunR. There are three options, A) using the R package devtools, B) cloning the github, C) manually downloading the .tar.gz file. Note that the Seq2FunR github will have the most up-to-date version of the package. 

#### Option A) Install the package directly from github using the *devtools* package. Open R and enter:

Due to issues with Latex, some users may find that they are only able to install Seq2FunR without any documentation (i.e. vignettes). 

```R
# Step 1: Install devtools
install.packages("devtools")
library(devtools)

# Step 2: Install MicrobiomeAnalystR WITHOUT documentation
devtools::install_github("xia-lab/Seq2FunR")


```

## Case Studies
still in developing, need to make it as a built in data;

```R
setwd("/working/dir/having/your/reads");
seq2fun(suffixNameR1 = "_R1.fastq.gz",
        suffixNameR2 = "_R2.fastq.gz",
        genemap = "database/gene_taxa_matched.txt",
        tfmi = "database/twins_IGC.geneset.out.renew.matched.clean.taxon.fmi");
```

#For reads quality

```R
to get the graph of quality scores;

getReadsQuality(inputFile = "outputSeq2fun/SRR6468595_report.json")


to get graph of insert size distribution;

getInsertSizeDistribution(inputFile = "outputSeq2fun/SRR6468595_report.json")


to get reads duplication levels;

getDuplicationLevel(inputFile = "outputSeq2fun/SRR6468595_report.json")

```






