#Preprocessing of Large Hematopoeisis 220k Cells
#04/22/20
#Adapted from Granja*, Corces*, et al. 
#ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis (2020)
#Created by Jeffrey Granja

#Analysis was performed with version 0.2.1 so some function names and defaults may have changed in current release version
#We document most of these analysis on our tutorial data with better documentation and updated functions in ArchR's
#user manual (https://www.archrproject.com/bookdown/index.html). Please refer here to perform these analyses with the most
#up to date version for best functionality.

library(ArchR)

dir.create("/Volumes/JG_SSD_2/LargeHeme/") #Working Directory
setwd("/Volumes/JG_SSD_2/LargeHeme/")

inputFiles <- system("find /Volumes/JG_SSD_2/Analysis/hg19/ -name *.gz*", intern = TRUE)
names(inputFiles) <- gsub(".fragments.tsv.gz|.tsv.gz", "", basename(inputFiles))

addArchRGenome("hg19")
addArchRThreads(16, force = TRUE)

tstart <- Sys.time()
ArrowFiles <- createArrowFiles(
	inputFiles = inputFiles,
	minFrags = 1000,
	maxFrags = 50000,
	filterFrags = 1000,
	filterTSS = 4
)
saveTime <- Sys.time() - tstart
cat(utils:::capture.output(saveTime), file = "Save-Time-CreateArrows.txt")

tstart2 <- Sys.time()
proj <- ArchRProject(ArrowFiles = ArrowFiles, copyArrows = FALSE)
saveTime2 <- Sys.time() - tstart2
cat(utils:::capture.output(saveTime2), file = "Save-Time-CreateArchRProject.txt")

tstart3 <- Sys.time()
proj <- addIterativeLSI(
	ArchRProj = proj, 
	varFeatures = 25000,
	sampleCellsFinal = 25000,
	saveIterations = FALSE,
	force = TRUE
)
saveTime3 <- Sys.time() - tstart3
cat(utils:::capture.output(saveTime3), file = "Save-Time-ReducedDims.txt")

tstart4 <- Sys.time()
proj <- addClusters(proj, force = TRUE)
saveTime4 <- Sys.time() - tstart4
cat(utils:::capture.output(saveTime4), file = "Save-Time-Clustering.txt")

tstart5 <- Sys.time()
proj <- addUMAP(proj, name = "UMAP", force = TRUE)
saveTime5 <- Sys.time() - tstart5
cat(utils:::capture.output(saveTime5), file = "Save-Time-UMAP.txt")

save.image("Save-Analysis-LargeHeme.Rdata", compress = FALSE)




