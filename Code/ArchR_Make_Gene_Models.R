#Creating Gene Score Models for ArchR
#04/22/20
#Adapted from Granja*, Corces*, et al. 
#ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis (2020)
#Created by Jeffrey Granja

#Analysis was performed with version 0.2.1 so some function names and defaults may have changed in current release version
#We document most of these analysis on our tutorial data with better documentation and updated functions in ArchR's
#user manual (https://www.archrproject.com/bookdown/index.html). Please refer here to perform these analyses with the most
#up to date version for best functionality.

library(ArchR)
#addArchRGenome("hg19")
addArchRGenome("hg38")

genes <- getGenes()
mcols(genes)$name <- mcols(genes)$symbol

#############################################
# Model Promoter Regions
#############################################
dropStrand <- function(gr){
	strand(gr) <- "*"
	gr
}

list(
	modelName = "Promoter_1K",
	FN = "addFeatureMatrix",
	features = dropStrand(resize(resize(genes, 1, "start"), 1000, "center")),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-Promoter-1.rds")

list(
	modelName = "Promoter_2K",
	FN = "addFeatureMatrix",
	features = dropStrand(resize(resize(genes, 1, "start"), 2000, "center")),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-Promoter-2.rds")

list(
	modelName = "Promoter_5K",
	FN = "addFeatureMatrix",
	features = dropStrand(resize(resize(genes, 1, "start"), 5000, "center")),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-Promoter-5.rds")

list(
	modelName = "Promoter_10K",
	FN = "addFeatureMatrix",
	features = dropStrand(resize(resize(genes, 1, "start"), 10000, "center")),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-Promoter-10.rds")

list(
	modelName = "Promoter_25K",
	FN = "addFeatureMatrix",
	features = dropStrand(resize(resize(genes, 1, "start"), 25000, "center")),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-Promoter-25.rds")

list(
	modelName = "Promoter_50K",
	FN = "addFeatureMatrix",
	features = dropStrand(resize(resize(genes, 1, "start"), 50000, "center")),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-Promoter-50.rds")

list(
	modelName = "Promoter_100K",
	FN = "addFeatureMatrix",
	features = dropStrand(resize(resize(genes, 1, "start"), 100000, "center")),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-Promoter-50.rds")

#############################################
# Model Promoter Regions
#############################################
dropStrand <- function(gr){
	strand(gr) <- "*"
	gr
}

list(
	modelName = "GeneBody_0_0",
	FN = "addFeatureMatrix",
	features = dropStrand(extendGR(genes, 0, 0)),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-GeneBody-0_0.rds")

list(
	modelName = "GeneBody_1000_0",
	FN = "addFeatureMatrix",
	features = dropStrand(extendGR(genes, 1000, 0)),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-GeneBody-1000_0.rds")

list(
	modelName = "GeneBody_2000_0",
	FN = "addFeatureMatrix",
	features = dropStrand(extendGR(genes, 2000, 0)),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-GeneBody-2000_0.rds")

list(
	modelName = "GeneBody_5000_0",
	FN = "addFeatureMatrix",
	features = dropStrand(extendGR(genes, 5000, 0)),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-GeneBody-5000_0.rds")

list(
	modelName = "GeneBody_10000_0",
	FN = "addFeatureMatrix",
	features = dropStrand(extendGR(genes, 10000, 0)),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-GeneBody-10000_0.rds")

list(
	modelName = "GeneBody_1000_1000",
	FN = "addFeatureMatrix",
	features = dropStrand(extendGR(genes, 1000, 1000)),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-GeneBody-1000_1000.rds")

list(
	modelName = "GeneBody_2000_2000",
	FN = "addFeatureMatrix",
	features = dropStrand(extendGR(genes, 2000, 2000)),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-GeneBody-2000_2000.rds")

list(
	modelName = "GeneBody_5000_5000",
	FN = "addFeatureMatrix",
	features = dropStrand(extendGR(genes, 5000, 5000)),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-GeneBody-5000_5000.rds")

list(
	modelName = "GeneBody_10000_10000",
	FN = "addFeatureMatrix",
	features = dropStrand(extendGR(genes, 10000, 10000)),
	matrixName = "FeatureMatrix",
	force = TRUE
) %>% saveRDS("Models/Model-GeneBody-10000_10000.rds")


#############################################
# Custom Gene Models From TSS
#############################################

i <- 1

list(
	modelName = paste0("GeneModel-TSS-Constant-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "1",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 5000),
	extendDownstream = c(1000, 5000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-Constant-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-TSS-Constant-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "1",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 10000),
	extendDownstream = c(1000, 10000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-Constant-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-TSS-Constant-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "1",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 25000),
	extendDownstream = c(1000, 25000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-Constant-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-TSS-Constant-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "1",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-Constant-",i,".rds"))


#############################################
# Custom Gene Models From TSS 
#############################################

i <- 1

list(
	modelName = paste0("GeneModel-TSS-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-TSS-Exponential-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-TSS-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-TSS-Exponential-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-TSS-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-TSS-Exponential-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-TSS-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/100000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-TSS-Exponential-",i,".rds"))



#############################################
# Custom Gene Models From TSS No Boundary
#############################################

i <- 1

list(
	modelName = paste0("GeneModel-TSS-NoBoundary-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-TSS-Exponential-NoBoundary-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-TSS-NoBoundary-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-TSS-Exponential-NoBoundary-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-TSS-NoBoundary-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-TSS-Exponential-NoBoundary-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-TSS-NoBoundary-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/100000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-TSS-Exponential-NoBoundary-",i,".rds"))

#############################################
# Custom Gene Models From GB 
#############################################

i <- 1

list(
	modelName = paste0("GeneModel-GB-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/100000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-",i,".rds"))


#############################################
# Custom Gene Models From TSS No Boundary
#############################################

i <- 1

list(
	modelName = paste0("GeneModel-GB-NoBoundary-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-NoBoundary-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-NoBoundary-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-NoBoundary-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-NoBoundary-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-NoBoundary-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-NoBoundary-Exponential-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/100000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-NoBoundary-",i,".rds"))


#############################################
# Custom Gene Models From GB 
#############################################

i <- 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 1000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 2000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 5000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 1000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 2000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 5000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 1000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 2000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))

i <- i + 1

list(
	modelName = paste0("GeneModel-GB-Exponential-Extend-", i),
	FN = "addGeneScoreMatrix",
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 5000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5 #New Param
) %>% saveRDS(paste0("Models/GeneModel-GB-Exponential-Extend-",i,".rds"))






























































































































































































































