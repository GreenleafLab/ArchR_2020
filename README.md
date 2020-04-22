# ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis (Granja JM*, Corces MR*, et al. 2020)

## **Link** : QQQ

## Please cite : Granja JM et al., An integrative and scalable software package for single-cell chromatin accessibility analysis. bioRxiv (2020) <br/>

![](Figure1.png)

# Brief Descriptions of Analysis Scripts

## scATAC Analyses (See Code Folder)

**ArchR_Make_Gene_Models.R** - Code for making individual gene score models with ArchR.

**ArchR_Test_Gene_Score_Model_Aggregates.R** - Code for testing gene score models with aggregates of cells.

**ArchR_PreProcess_Laptop_SimulatedPBMC_1M.R** - Code for benchmarking analysis of 1M Simulated PBMC.

**ArchR_PreProcess_Laptop_LargeHeme_220k.R** - Code for benchmarking analysis of Large Hematopoeisis.

**ArchR_PostAnalysis_Large_Heme_Analysis.R** - Code for advanced post analysis of Large Hematopoeisis.

# scATAC-seq Files

**10x Version1 vs NextGem (Fragment Files)** - https://www.10xgenomics.com/solutions/single-cell-atac/

**Satpathy*, Granja* et al 2019 (Fragment Files)** - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785&holding=F1000&otool=stanford

**Granja*, Klemm*, McGinnis* et al 2019 (Fragment Files)** - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139369

**Buenrostro et al 2018 (Bam Files)** - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96772

**Cusanovich et al 2018 (Bam Files)** - http://atlas.gs.washington.edu/mouse-atac/data/

# Benchmarking Results (ArchR, SnapATAC, Signac)

**PBMC 10k** - https://github.com/GreenleafLab/ArchR_2020/tree/master/ArchR_Benchmarking/PBMC-10k

**PBMC 20k** - https://github.com/GreenleafLab/ArchR_2020/tree/master/ArchR_Benchmarking/PBMC-20k

**PBMC 30k** - https://github.com/GreenleafLab/ArchR_2020/tree/master/ArchR_Benchmarking/PBMC-30k

**PBMC 70k** - https://github.com/GreenleafLab/ArchR_2020/tree/master/ArchR_Benchmarking/PBMC-70k

**BMMC 30k** - https://github.com/GreenleafLab/ArchR_2020/tree/master/ArchR_Benchmarking/BMMC-30k

**Mouse Atlas 70k** - https://github.com/GreenleafLab/ArchR_2020/tree/master/ArchR_Benchmarking/MouseAtlas-70k/ArchR

**Note 1.** We included 1 replicate Rmarkdown for the large computational setup since we have results for all softwares for each of the steps.

**Note 2.** We merged steps as described in the supplemental of the MS for clarity.

# Additional Data Download Links

### These links may be moved if we can find a better host for better download speed

## Notes

**.rds** file is an R binarized object to read into R use readRDS(filename)

**SummarizedExperiment** is a class in R see : <br/>https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html

**deviations** (TF chromVAR) is a class in R see : <br/>https://bioconductor.org/packages/release/bioc/html/chromVAR.html

## Large Hematopoiesis

**ArchR Project Containing Arrow Files** : <br/>https://drive.google.com/file/d/1L8hOTLDxqdzgqY4hw_uXr5O2m8qd9Bdn/view?usp=sharing

**scATAC-seq Hematopoeisis cell x peak Summarized Experiment** : <br/>https://drive.google.com/file/d/1zWFJes4z2uMOkgm4kEd7zffjlgBnlSxD/view?usp=sharing

**scATAC-seq Hematopoeisis cell x motif Summarized Experiment** : <br/>https://drive.google.com/file/d/1RUOcJbxric6PF3Ilzx4622uS2Lbi8r31/view?usp=sharing

