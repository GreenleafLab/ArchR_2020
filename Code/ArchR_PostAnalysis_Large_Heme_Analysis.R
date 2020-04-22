#Analysis of Large Hematopoietic data with ArchR
#04/22/20
#Adapted from Granja*, Corces*, et al. 
#ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis (2020)
#Created by Jeffrey Granja

#Analysis was performed with version 0.2.1 so some function names and defaults may have changed in current release version
#We document most of these analysis on our tutorial data with better documentation and updated functions in ArchR's
#user manual (https://www.archrproject.com/bookdown/index.html). Please refer here to perform these analyses with the most
#up to date version for best functionality.

library(ArchR) 
addArchRThreads(20)
addArchRGenome("hg19")

set.seed(1)

#See HTML files for Preprocessing to this step on github
proj <- readRDS("All_Immune_Filter3-ArchR-3-ReduceDimsNoDoub2-25000-Save.rds")
proj@projectMetadata$outputDirectory <- "output"

############################################################
# Cluster / Marker Gene Analysis
############################################################

#Add Clusters with default resolution 0.8
proj <- addClusters(proj, force = TRUE) #23 minutes

#Embedding was already computed see HTML file we now can overlay cluster IDs
p <- plotEmbedding(proj, name = "Clusters")
plotPDF(p, name = "PlotClusters")

#We want to better Resolve The Sub-Clustering AT CD34+ Cells After Inspection
projCD34 <- proj[proj$Clusters %in% c(paste0("C", 16:21))]
projCD34 <- addClusters(projCD34, name = "Cluster34", force = TRUE, resolution = 0.4)

#Create New Clusters
proj$NewClust <- proj$Clusters
proj$NewClust[match(projCD34$cellNames, proj$cellNames)] <- paste0("CD34_", projCD34$Cluster34)

p <- plotEmbedding(proj, name = "NewClust")
plotPDF(p, name = "PlotClusters")

#Plot Some Marker Genes To See if We Can ID The Clusters
proj <- addImputeWeights(proj, nRep = 5)

#Markers
markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", #B-Cell Trajectory
    "CD14", #Monocytes
    "FCGR3A", "GZMB", "SDC1", "GATA2", "CD1C", "CD4",
    "SELL", "GATA2",
    "CD3D", "CD8A", "TBX21", "IL7R", #TCells,
    "B3GAT1", "CD69", "NCAM1", "CD38", "IL2RA", "KLRB1"
  )

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP", 
  bins = 100,
  continuousSet = "horizonExtra", 
  quantileCut = c(0, 0.999),
  log2norm = TRUE
)
plotPDF(p, name = "Plot-Markers")

#Plot Embedding with Re-Assignment of Clusters for CD34+
p <- plotEmbedding(proj, name = "NewClust", labelAsFactors = FALSE)
plotPDF(p, name = "PlotClusters")

#Last Round Merge Extraneous Clusters and Assign Biological Labels
renameList <- list(
  "C1_HSC" = c("CD34_C4"),
  "C2_CMP.LMPP" = c("CD34_C5", "CD34_C6", "CD34_C7"),
  "C3_Early.Ery" = c("CD34_C2"),
  "C4_Late.Ery" = c("CD34_C1", "CD34_C3"),
  "C5_Early.Baso" = paste0("C",15),
  "C6_Late.Baso" = paste0("C",1),
  "C7_CLP.1" = c("CD34_C11"),
  "C8_CLP.2" = paste0("C",13),
  "C9_PreB" = paste0("C",14),
  "C10_B" = paste0("C",c(11,12)),
  "C11_Plasma" = paste0("C",c(10)),
  "C12_pDC" = c("CD34_C12"),
  "C13_GMP" = c("CD34_C8", "CD34_C9", "CD34_C10"),
  "C14_Mono" = paste0("C",c(9, 33, 8, 3, 7, 5, 2, 6)),
  "C15_cDC" = paste0("C",c(4)),
  "C16_CD8N" = paste0("C",c(23)),
  "C17_CD4N" = paste0("C",c(24, 26, 25, 27)),
  "C18_CD4M" = paste0("C",c(29, 28, 32, 22, 31)),
  "C19_CD8_EM" = paste0("C",c(30, 35, 36, 38)),
  "C20_CD8_CM" = paste0("C",c(34)),
  "C21_NK" = paste0("C",c(37, 39))
)

namesV <- lapply(seq_along(renameList), function(x){
  data.frame(val = renameList[[x]], name = names(renameList)[x])
}) %>% Reduce("rbind", .)

saveRDS(namesV, "Remap.rds")

proj$ClustFinal <- mapLabels(proj$NewClust, newLabels = paste0(namesV[,2]), oldLabels = paste0(namesV[,1]))
p <- plotEmbedding(proj, name = "ClustFinal", labelAsFactors = TRUE)
plotPDF(p, name = "PlotClusters")

#Create A nice color Palette for Heme
palHeme <- c(
"FAC0FF
002BFF
C1B62B
7E2AD8
E26DDF
A9A9A9
FF6700
602E01
FF0000
FFB55C
82FC56
A37857
1FBF82
A30D0D
FCE017
FF7D7D
A983F2
036D03
AA059F
43B9F9
262C6B
") %>% {stringr::str_split(., pattern="\n", simplify=TRUE)[1,]}
palHeme <- palHeme[palHeme!=""]
palHeme <- paste0("#", palHeme)
names(palHeme) <- gtools::mixedsort(unique(proj$ClustFinal))

#Plot Embedding with New Palette and Biological Labels
p <- plotEmbedding(proj, name = "ClustFinal", labelAsFactors = TRUE, pal = palHeme)
plotPDF(p, name = "PlotClusters", width = 20, height= 20)

############################################################
# Make A Peak Set and Counts
############################################################
proj <- addGroupCoverages(proj, groupBy = "ClustFinal")
proj <- addReproduciblePeakSet(proj, groupBy = "ClustFinal")
proj <- addPeakMatrix(proj)

############################################################
# Add Motif Deviations Matrix
############################################################
proj <- addMotifAnnotations(proj, version = 1)
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(proj, force = TRUE, threads = 10)

############################################################
# Identify Motif Regulators
############################################################
corTFGS <- correlateMatrices(proj, useMatrix1 = "MotifMatrix", useMatrix2 = "GeneScoreMatrix")
corTFRNA <- correlateMatrices(proj, useMatrix1 = "MotifMatrix", useMatrix2 = "GeneIntegrationMatrix")

#Group Motif Deviations
seGroupMotif <- exportGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "ClustFinal")

#Get Deviation Z-Scores
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
seZ <- seZ[rowData(seZ)$name %in% corTFGS[,1], ]

#Max Delta
maxD <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs
corTFGS$maxDelta <- maxD[match(corTFGS[,1], rowData(seZ)$name)]


#Filter By Absolute Correlation and deduplicate
corTFGS <- corTFGS[order(abs(corTFGS$cor), decreasing = TRUE), ]
corTFGS <- corTFGS[which(!duplicated(gsub("\\-.*","",corTFGS[,2]))), ]

#Identify Sig Relation Ships
corTFGS2 <- corTFGS[which(corTFGS$cor > 0.5 & corTFGS$padj < 0.001), ]
corTFGS$sig <- "No"
corTFGS$sig[which(corTFGS[,1] %in% corTFGS2[,1])] <- "Yes"
corTFGS$sig[which(corTFGS$cor < 0.5)] <- "No"
corTFGS$sig[corTFGS$maxDelta < quantile(corTFGS$maxDelta, 0.5)] <- "No"

pdf("Plot-Cor-TF-GS-2.pdf", width = 4, height = 4, useDingbats = FALSE)
ggplot(data.frame(corTFGS), aes(cor, maxDelta, color = sig)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0) + 
  scale_color_manual(values = c("No"="darkgrey", "Yes"="firebrick3")) +
  xlab("Correlation To GS") +
  ylab("TF Motif Delta") +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(corTFGS$maxDelta)*1.05))
dev.off()

readr::write_tsv(data.frame(corTFGS), "Save-Correlations-TF-To-GS-2.tsv")

############################################################
# Marker Peak Identification
############################################################
seMarkerGenes <- markerFeatures(proj, groupBy = "ClustFinal", maxCells = 1000)
seMarkerPeaks <- markerFeatures(proj, groupBy = "ClustFinal", maxCells = 1000, useMatrix = "PeakMatrix")

p <- markerHeatmap(seMarkerGenes, transpose = TRUE, limits = c(-1.5, 1.5),  cutOff = "FDR <= 0.01 & Log2FC >= 1.5",
  clusterCols = FALSE, pal = paletteContinuous("horizonExtra"))
plotPDF(p, name = "MarkerHeamtapGenes", width = 12, height = 8)

p <- markerHeatmap(seMarkerPeaks, transpose = TRUE, limits = c(-1.5, 1.5), clusterCols = FALSE , cutOff = "FDR <= 0.01 & Log2FC >= 1")
plotPDF(p, name = "MarkerHeamtapPeaks", width = 12, height = 8, addDOC = FALSE)

annoEnrich <- peakAnnoEnrichment(seMarkerPeaks, proj, peakAnnotation = "Motif", cutOff = "FDR <= 0.01 & Log2FC >= 1")
df <- data.frame(readr::read_tsv("Save-Correlations-TF-To-GS.tsv"))
df2 <- df[df$sig=="Yes",]

p <- enrichHeatmap(annoEnrich[df2[,1], ], transpose = TRUE, clusterCols = FALSE, n = 4, rastr = FALSE)
plotPDF(p, name = "PeakAnnoEnrichHeatmap.pdf", width = 6, height = 6)


proj <- addArchRAnnotations(proj, collection="ATAC")
annoEnrich <- peakAnnoEnrichment(seMarkerPeaks, proj, "ATAC")
p <- enrichHeatmap(annoEnrich, transpose = TRUE, clusterCols = FALSE, n = 6, rastr = FALSE)
plotPDF(p, name = "PeakAnnoEnrichHeatmap-Peaks", width = 6, height = 6)

############################################################
# Motif Footprints of Select Motifs
############################################################
p <- enrichHeatmap(annoEnrich[df2[,1], ], transpose = TRUE, clusterCols = FALSE, n = 4, returnMatrix = TRUE)

motifPositions <- getPositions(proj)

seFoot_None <- plotFootprints(
  ArchRProj = proj, 
  positions = motifPositions[stringr::str_split(colnames(p),pattern=" ", simplify=TRUE)[,1]], 
  groupBy = "ClustFinal",
  normMethod = "none",
  plotName = "Footprints-No-Normalization"
)

seFoot_None2 <- plotFootprints(
  ArchRProj = proj, 
  inputSE = seFoot_None,
  positions = motifPositions[stringr::str_split(colnames(p),pattern=" ", simplify=TRUE)[,1]], 
  groupBy = "ClustFinal",
  normMethod = "none",
  plotName = "Footprints-No-Normalization",
  smoothWindow = NULL
)

seFoot_Divide <- plotFootprints(
  ArchRProj = proj, 
  inputSE = seFoot_None,
  positions = motifPositions[stringr::str_split(colnames(p),pattern=" ", simplify=TRUE)[,1]], 
  groupBy = "ClustFinal",
  normMethod = "divide",
  pal = palHeme,
  plotName = "Footprints-Divide-Bias",
  smoothWindow = NULL
)

seFoot_Subtract <- plotFootprints(
  ArchRProj = proj, 
  inputSE = seFoot_None,
  positions = motifPositions[stringr::str_split(colnames(p),pattern=" ", simplify=TRUE)[,1]], 
  groupBy = "ClustFinal",
  normMethod = "subtract",
  pal = palHeme,
  plotName = "Footprints-Subtract-Bias",
  smoothWindow = NULL
)

############################################################
# Plot UMAP w/ JDB Labels
############################################################
idx2 <- grep("JDB_Cell", proj$Sample)
proj$Buenrostro <- "NA"
proj$Buenrostro[idx2] <- stringr::str_split(gsub("JDB_Cell#","",proj$cellNames[idx2]),pattern="\\.",simplify=TRUE)[,1]
palB <- paletteDiscrete(stringr::str_split(gsub("JDB_Cell#","",proj$cellNames[idx2]),pattern="\\.",simplify=TRUE)[,1])
palB[["NA"]] <- "lightgrey"
plotDF <- getEmbedding(proj)
plotDF <- plotDF[proj$cellNames, ]
plotDF$color <- proj$Buenrostro
idx <- c(which(plotDF[,3]=="NA"), which(plotDF[,3]!="NA"))
p <- ggPoint(
  x = plotDF[idx,1], 
  y = plotDF[idx,2], 
  color = paste0(plotDF[idx,3]), 
  pal = palB, 
  rastr = TRUE, 
  xlabel = "UMAP 1", 
  ylabel = "UMAP 2", 
  size = 0.2
)
plotPDF(p, name = "Plot-Buenrostro", width = 8, height = 8)

########################################################################################################################
########################################################################################################################
########################################################################################################################

############################################################
# RNA-Integration
############################################################

seRNA <- readRDS("/scratch/users/jgranja/OtherAnalysis/GeneScores/Heme/scRNA-Healthy-Hematopoiesis-191120.rds")

groupList <- SimpleList(
  list( #Not T/NK-Cells
    ATAC = grep(paste0(paste0("C",1:15,"_"),collapse="|"), unique(proj$ClustFinal), value = TRUE),
    RNA = names(table(colData(seRNA)$BioClassification))[1:18]
  ),
  list( #T/NK-Cells
    ATAC = grep(paste0(paste0("C",16:21,"_"),collapse="|"), unique(proj$ClustFinal), value = TRUE),
    RNA = names(table(colData(seRNA)$BioClassification))[19:25]
  )
)

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  seRNA = seRNA,
  groupRNA = "BioClassification",
  groupATAC = "ClustFinal",
  groupList = groupList,
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  embeddingATAC = getEmbedding(proj, embedding = "UMAP"),
  embeddingRNA = colData(seRNA)[, c("UMAP1", "UMAP2")],
  matrixName = "GeneIntegrationMatrix",
  force = TRUE,
  addToArrow = TRUE
)

#Plot Results for Cells Above Prediction Score of 0.4
proj@embeddings[[2]] <- proj@embeddings[[1]]
proj@embeddings[[2]][[1]] <- proj@embeddings[[2]][[1]][proj$cellNames[proj$predictedScore > 0.4], ]
names(proj@embeddings) <- c("UMAP", "UMAP2")
p <- plotEmbedding(proj, name = "predictedGroup", labelAsFactors = TRUE, embedding = "UMAP2")
plotPDF(p, name = "PlotPredictedGroup", width = 10, height= 10)

############################################################
# Plot UMAP With GeneScores and RNA for Marker Genes
############################################################

#Markers
markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", #B-Cell Trajectory
    "CD14", #Monocytes
    "FCGR3A", "GZMB", "SDC1", "GATA2", "CD1C", "CD4",
    "SELL", "GATA2",
    "CD3D", "CD8A", "TBX21", "IL7R", #TCells,
    "B3GAT1", "CD69", "NCAM1", "CD38", "IL2RA", "KLRB1",
    "AVP", "DERL3", "MPO", "AZU1", "MME", "BLVRB", "CCL5",
    "CCR7"
  )

p1 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP", 
  continuousSet = "horizonExtra", 
  quantileCut = c(0, 0.999),
  log2Norm = TRUE
)

plotPDF(p1, name = "Plot-Markers-GS")

p2 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneIntegrationMatrix", 
  name = markerGenes, 
  embedding = "UMAP", 
  continuousSet = "blueYellow", 
  quantileCut = c(0, 0.999),
  log2Norm = TRUE
)

plotPDF(p2, name = "Plot-Markers-RNA")

############################################################
# Identify Motif Regulators
############################################################

#Get Deviation Z-Scores
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
seZ <- seZ[rowData(seZ)$name %in% corTFRNA[,1], ]

#Max Delta
maxD <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs
corTFRNA$maxDelta <- maxD[match(corTFRNA[,1], rowData(seZ)$name)]

#Filter By Absolute Correlation and deduplicate
corTFRNA <- corTFRNA[order(abs(corTFRNA$cor), decreasing = TRUE), ]
corTFRNA <- corTFRNA[which(!duplicated(gsub("\\-.*","",corTFRNA[,2]))), ]
corTFRNA <- corTFRNA[]

#Identify Sig Relation Ships
corTFRNA2 <- corTFRNA[which(corTFRNA$cor > 0.5 & corTFRNA$padj < 0.001), ]
corTFRNA$sig <- "No"
corTFRNA$sig[which(corTFRNA[,1] %in% corTFRNA2[,1])] <- "Yes"
corTFRNA$sig[which(corTFRNA$cor < 0.5)] <- "No"
corTFRNA$sig[corTFRNA$maxDelta < quantile(corTFRNA$maxDelta, 0.5)] <- "No"

pdf("Plot-Cor-TF-RNA.pdf", width = 4, height = 4, useDingbats = FALSE)
ggplot(data.frame(corTFRNA), aes(cor, maxDelta, color = sig)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0) + 
  scale_color_manual(values = c("No"="darkgrey", "Yes"="firebrick3")) +
  xlab("Correlation To RNA") +
  ylab("TF Motif Delta") +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(corTFRNA$maxDelta)*1.05))
dev.off()

readr::write_tsv(data.frame(corTFRNA), "Save-Correlations-TF-To-RNA.tsv")

df <- data.frame(sig = intersect(corTFRNA[corTFRNA$sig=="Yes",2], corTFGS[corTFGS$sig=="Yes",2]))

readr::write_tsv(data.frame(df), "Save-Correlations-TF-To-RNA-GS-Overlap.tsv")

############################################################
# Plot Motif Regulators With GeneScores, RNA and chromVAR
############################################################

#Markers
markerGenes  <- c(
  readr::read_tsv("Save-Correlations-TF-To-RNA-GS-Overlap.tsv")[,1][[1]]
)

p1 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP", 
  continuousSet = "horizonExtra", 
  quantileCut = c(0, 0.999),
  log2Norm = TRUE
)

plotPDF(p1, name = "Plot-Cor-TF-GS")

p2 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneIntegrationMatrix", 
  name = markerGenes, 
  embedding = "UMAP", 
  continuousSet = "blueYellow", 
  quantileCut = c(0, 0.999),
  log2Norm = TRUE
)

plotPDF(p2, name = "Plot-Cor-TF-RNA")

corTFGS <- data.frame(readr::read_tsv("Save-Correlations-TF-To-GS.tsv"))
df <- data.frame(readr::read_tsv("Save-Correlations-TF-To-RNA-GS-Overlap.tsv"))
markerMotifs <- c(
  corTFGS[match(paste0(df$sig), corTFGS[,2]),1]
)

p3 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "MotifMatrix", 
  name = paste0("z:",markerMotifs), 
  embedding = "UMAP", 
  continuousSet = "solarExtra", 
  quantileCut = c(0.01, 0.99),
  log2Norm = FALSE
)

plotPDF(p3, name = "Plot-Cor-TFs")

############################################################
# Identify Peak To Gene Linkages and Co-Accessibility
############################################################

#Peak-To-Gene Linkage
proj <- addPeak2GeneLinks(proj, knnIteration = 1000)
p <- peak2GeneHeatmap(proj, nPlot = 100000, groupBy = "ClustFinal", palGroup = palHeme)
plotPDF(p, name = "Peak2GeneHeatmap", width = 14, height = 8)

#Add CoAccessbility
proj <- addCoAccessibility(proj, knnIteration = 1000)

############################################################
# Cellular Trajectory Analysis for Lymphoid
############################################################

trajectory <- c("C1_HSC", "C2_CMP.LMPP", "C7_CLP.1", "C8_CLP.2", "C9_PreB", "C10_B")
proj <- addTrajectory(
  proj, 
  name = "Lymphoid", 
  trajectory = trajectory,
  groupBy = "ClustFinal",
  force = TRUE,
  preFilterQuantile = 0.8,
  postFilterQuantile = 0.9,
  useAll = FALSE
)
ArchRPalettes$fireworks2 <- ArchRPalettes$fireworks[-1]
p1 <- plotTrajectory(proj, trajectory = "Lymphoid", name = "Lymphoid", smoothWindow = 11, continuousSet = "greenBlue")
p2 <- plotTrajectory(proj, trajectory = "Lymphoid", name = "Lymphoid", smoothWindow = 11, continuousSet = "greyMagma")
p3 <- plotTrajectory(proj, trajectory = "Lymphoid", name = "Lymphoid", smoothWindow = 11, continuousSet = "purpleOrange")
p4 <- plotTrajectory(proj, trajectory = "Lymphoid", name = "Lymphoid", smoothWindow = 11, continuousSet = "beach")
p5 <- plotTrajectory(proj, trajectory = "Lymphoid", name = "Lymphoid", smoothWindow = 11, continuousSet = "sambaNight")
plotPDF(p3, name = "Plot-Lym-Trajectory")

#Correlated TF Heatmaps
lymMotifs  <- getTrajectory(ArchRProj = proj, name = "Lymphoid", useMatrix = "MotifMatrix", log2Norm = FALSE, smoothWindow = 11)
lymGS  <- getTrajectory(ArchRProj = proj, name = "Lymphoid", useMatrix = "GeneScoreMatrix", log2Norm = TRUE, smoothWindow = 11)
lymRNA  <- getTrajectory(ArchRProj = proj, name = "Lymphoid", useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE, smoothWindow = 11)
lymATAC  <- getTrajectory(ArchRProj = proj, name = "Lymphoid", useMatrix = "PeakMatrix", log2Norm = TRUE, smoothWindow = 11)

############################################################
#Correlate TFs-RNA-GS across Lymphoid
############################################################
corMR <- correlateTrajectories(lymMotifs[grep("z:",rownames(lymMotifs)),], lymRNA, useRanges = FALSE)
corMG <- correlateTrajectories(lymMotifs[grep("z:",rownames(lymMotifs)),], lymGS, useRanges = FALSE)

corDF <- corMR[[1]][,1:3]
corDF$TF <- rownames(corMR$seTrajectory1)[corMR[[1]][,1]]
corDF$Gene <- rownames(corMR$seTrajectory2)[corMR[[1]][,2]]
id1 <- match(corDF$TF, rownames(corMG$seTrajectory1))
id2 <- match(corDF$Gene, rownames(corMG$seTrajectory2))
corDF$CorrelationGS <- corMG[[2]][match(paste0(id1,"#",id2), apply(corMG[[2]][,1:2],1,function(x)paste0(x,collapse="#"))),3]
corDF2 <- corDF[which(corDF$CorrelationGS > 0.5),]
corDF2$CorAvg <- rowMeans(as.matrix(corDF2[,c("Correlation", "CorrelationGS")]))
corDF2 <- corDF2[order(corDF2$CorAvg, decreasing = TRUE), ]
corDF3 <- corDF2[!duplicated(corDF2[,"Gene"]), ]

rowOrder <- trajectoryHeatmap(lymRNA[corDF3$Gene,], maxFeatures = NULL, varCutOff = NULL, pal = paletteContinuous(set = "blueYellow"), returnMat = TRUE)
idxOrder <- match(rownames(rowOrder), corDF3$Gene)

corDF3[idxOrder,]

rowOrder <- trajectoryHeatmap(lymMotifs[corDF3$TF,], maxFeatures = NULL, varCutOff = NULL, pal = paletteContinuous(set = "blueYellow"), returnMat = TRUE)
idxOrder <- match(rownames(rowOrder), corDF3$TF)

tfHeatmap_Motifs <- trajectoryHeatmap(lymMotifs[corDF3$TF,], 
  rowOrder = idxOrder, 
  limits = c(-1.5, 1.5), 
  pal = paletteContinuous(set = "solarExtra")
)

tfHeatmap_RNA <- trajectoryHeatmap(
  lymRNA[corDF3$Gene,], 
  rowOrder = idxOrder, 
  limits = c(-1.5, 1.5), 
  pal = paletteContinuous(set = "blueYellow")
)

tfHeatmap_GS <- trajectoryHeatmap(
  lymGS[corDF3$Gene,], 
  rowOrder = idxOrder, 
  limits = c(-1.5, 1.5), 
  pal = paletteContinuous(set = "horizonExtra")
)

p <- tfHeatmap_Motifs + tfHeatmap_RNA + tfHeatmap_GS
plotPDF(p, name = "Heatmaps-Lym-TF-MRG", width = 12, height = 6, addDOC = FALSE)

############################################################
#Correlate ATAC + RNA Peak-To-Gene Links
############################################################
corAR <- correlateTrajectories(lymATAC, lymRNA, useRanges = TRUE, varCutOff1 = 0.9, varCutOff2 = 0.9)

corPeaks <- corAR[[1]]

a <- lymATAC[corPeaks[,1],]
b <- lymRNA[corPeaks[,2],]
c <- a
assay(c) <- assay(c) * assay(b)
rownames(c) <- paste0(rownames(a), "#", rownames(b))
rowOrderP <- trajectoryHeatmap(c,  pal = paletteContinuous(set = "blueYellow"), returnMat = TRUE, varCutOff = 0)
idxOrderP <- match(rownames(rowOrderP), rownames(corAR$seTrajectory1)[corPeaks[,1]])

peakHeatmap_ATAC <- trajectoryHeatmap(corAR$seTrajectory1[corPeaks[,1],], 
  varCutOff = NULL,
  maxFeatures = NULL,
  rowOrder = idxOrderP, 
  limits = c(-2, 2), 
  pal = paletteContinuous(set = "solarExtra")
)

peakHeatmap_RNA <- trajectoryHeatmap(corAR$seTrajectory2[corPeaks[,2],], 
  varCutOff = NULL,
  maxFeatures = NULL,
  rowOrder = idxOrderP, 
  labelTop = 2000,
  limits = c(-2, 2), 
  pal = paletteContinuous(set = "blueYellow")
)
p <- peakHeatmap_ATAC + peakHeatmap_RNA
plotPDF(p, name = "Heatmaps-Lym-AR10", width = 12, height = 8, addDOC = FALSE)

ARPeaksDF <- corPeaks[idxOrderP, ]
ARPeaksDF$Peak <- rownames(corAR$seTrajectory1)[ARPeaksDF[,1]]
ARPeaksDF$Gene <- rownames(corAR$seTrajectory2)[ARPeaksDF[,2]]
saveRDS(data.frame(ARPeaksDF), "Save-Lymphoid-ATAC-RNA-Linked-P2G.rds")
readr::write_tsv(data.frame(ARPeaksDF), "Save-Lymphoid-ATAC-RNA-Linked-P2G.tsv")

############################################################
# Footprint across Cellular Trajectory
############################################################
proj$LymphoidTrajGroup2 <- paste0("T_", plyr::round_any(ifelse(proj$Lymphoid > 0, proj$Lymphoid - 0.001, proj$Lymphoid), 20, floor))
proj <- addGroupCoverages(proj, groupBy = "LymphoidTrajGroup2")

motifPositions <- getPositions(proj)
motifsFoot <- corDF3$TF
motifsFoot <- gsub("z:","", motifsFoot)

seFoot_None <- plotFootprints(
  ArchRProj = proj, 
  positions = motifPositions[motifsFoot], 
  groupBy = "LymphoidTrajGroup2",
  normMethod = "none",
  plotName = "Lymphoid-Footprints-No-Normalization"
)

palLym <- paletteContinuous("fireworks2", n = 5)
#palLym <- c("#581845", "#900C3F", "#C70039", "#FF5744", "#FFC30F")
names(palLym) <- c("T_0", "T_20", "T_40", "T_60", "T_80")

seFoot_Subtract <- plotFootprints(
  ArchRProj = proj, 
  inputSE = seFoot_None,
  positions = motifPositions[motifsFoot], 
  groupBy = "LymphoidTrajGroup2",
  normMethod = "subtract",
  useGroups = names(palLym),
  pal = palLym,
  plotName = "Lymphoid-Footprints-Subtract-Bias",
  smoothWindow = NULL
)

############################################################
# Plot RNA across Lymphoid Trajectory
############################################################

markerGenes <- c("HMGA1", "BLK")
p1 <- plotTrajectory(proj, trajectory = "Lymphoid", 
  colorBy = "GeneIntegrationMatrix", name = markerGenes[1], smoothWindow = 11, 
  continuousSet = "fireworks2", log2Norm = TRUE)
p2 <- plotTrajectory(proj, trajectory = "Lymphoid", 
  colorBy = "GeneIntegrationMatrix", name = markerGenes[2], smoothWindow = 11, 
  continuousSet = "fireworks2", log2Norm = TRUE)

p11 <- p1[[2]] + coord_equal(ratio = 0.5, xlim = c(0,100), ylim = c(0,6), expand = FALSE)
p21 <- p2[[2]] + coord_equal(ratio = 0.5, xlim = c(0,100), ylim = c(0,3.75), expand = FALSE)

plotPDF(p11, p21, name = "Plot-Lym-Trajectory-Markers")

############################################################
# Project Bulk ATAC
############################################################
peakFile <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74912/suppl/GSE74912_ATACseq_All_Counts.txt.gz"
hemeDF <- data.frame(data.table::fread(peakFile))
rr <- GRanges(hemeDF$Chr, IRanges(start = hemeDF$Start, end = hemeDF$End))
reheader <- readRDS(gzcon(url("https://chang-public-data.s3-us-west-1.amazonaws.com/2016_NatGen_ATAC-AML/NatGen_Heme_Reheader.rds")))
colnames(hemeDF) <- reheader
hemeMat <- hemeDF[, grep("Donor", colnames(hemeDF),value=TRUE)]
seHeme <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(hemeMat)), rowRanges = rr)
projectHeme <- projectBulkATAC(proj, seATAC = seHeme, n = 250)
projectHeme[[1]][,3] <- stringr::str_split(projectHeme[[1]]$Type,pattern="\\_",simplify=TRUE)[,2]
plotProj <- rbind(projectHeme[[2]], projectHeme[[1]])
pal <- paletteDiscrete(unique(as.vector(plotProj[,3])))
pal["scATAC"] <- "lightgrey"
pdf("Plot-Bulk-Heme-Overlay-Projection.pdf")
p <- ggPoint(plotProj[,1], plotProj[,2], as.vector(plotProj[,3]), rastr = TRUE, pal = pal)
.fixPlotSize(p)
dev.off()





