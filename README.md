# ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis (Granja JM*, Corces MR*, et al. 2020)

## **Link** : QQQ

## Please cite : Granja JM et al., An integrative and scalable software package for single-cell chromatin accessibility analysis. bioRxiv (2020) <br/>

![](Figure1.png)

# Brief Descriptions of Analysis Scripts

## scATAC Analyses

**scATAC_01** - 

**scATAC_02** - 

**scATAC_03** - 

# scATAC-seq Files

### 10x Version1 vs NextGem (Fragment Files) - https://www.10xgenomics.com/solutions/single-cell-atac/

### Satpathy*, Granja* et al 2019 (Fragment Files) - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785&holding=F1000&otool=stanford

### Granja*, Klemm*, McGinnis* et al 2019 (Fragment Files) - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139369

### Buenrostro et al 2018 (Bam Files) - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96772

### Cusanovich et al 2018 (Bam Files) - http://atlas.gs.washington.edu/mouse-atac/data/

# Benchmarking Rmarkdowns

### PBMC 10k -

### PBMC 20k -

### PBMC 30k -

### PBMC 70k -

### BMMC 30k -

### Mouse Atlas 70k -

**Note 1.** We included 1 replicate Rmarkdown for each of the steps.

**Note 2.** We merged steps as described in the supplemental of the MS for clarity.

# Additional Data Download Links

### These links may be moved if we can find a better host for better download speed

## Notes

**.rds** file is an R binarized object to read into R use readRDS(filename)

**SummarizedExperiment** is a class in R see : <br/>https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html

**deviations** (TF chromVAR) is a class in R see : <br/>https://bioconductor.org/packages/release/bioc/html/chromVAR.html

## Large Hematopoiesis

**ArchR Project Containing Arrow Files** : <br/>https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scATAC-Healthy-Hematopoiesis-191120.rds

**scATAC-seq Hematopoeisis cell x peak Summarized Experiment** : <br/>https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scATAC-Healthy-Hematopoiesis-191120.rds

**scATAC-seq Hematopoeisis cell x gene Summarized Experiment** : <br/>https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scATAC-Healthy-Hematopoiesis-191120.rds

