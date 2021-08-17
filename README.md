# ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis (Granja JM*, Corces MR*, et al. 2020)

## **Link** : https://www.nature.com/articles/s41588-021-00790-6

## Please cite : Granja JM et al., ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nature Genetics (2021) <br/>

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

**Satpathy+, Granja+ et al 2019 (Fragment Files)** - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785&holding=F1000&otool=stanford

**Granja+, Klemm+, McGinnis+ et al 2019 (Fragment Files)** - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139369

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

**ArchR Project Containing Arrow Files** : <br/>https://www.dropbox.com/s/sijf2votfej629t/Save-Large-Heme-ArchRProject.tar.gz?dl=0
md5sum - 48d2b201b26b31d432314cc7dc8eb4a2

**scATAC-seq Hematopoeisis cell x peak Summarized Experiment** : <br/>https://www.dropbox.com/s/8ghfwumepftu3ll/SE-scATAC-Large-Heme-Peaks.rds?dl=0
md5sum - 7c8abb4edf30b8aef038ddfc18f65af0

**scATAC-seq Hematopoeisis cell x motif Summarized Experiment** : <br/>https://www.dropbox.com/s/owcathhdbcr23jq/SE-scATAC-Large-Heme-Motifs.rds?dl=0
md5sum - 23b7f6b7cdaba5e92ce695b309ec993c

