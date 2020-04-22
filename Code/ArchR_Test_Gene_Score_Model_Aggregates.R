#Comparison of Gene Score Models using Aggregates of Cells
#04/22/20
#Adapted from Granja*, Corces*, et al. 
#ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis (2020)
#Created by Jeffrey Granja

#Analysis was performed with version 0.2.1 so some function names and defaults may have changed in current release version
#We document most of these analysis on our tutorial data with better documentation and updated functions in ArchR's
#user manual (https://www.archrproject.com/bookdown/index.html). Please refer here to perform these analyses with the most
#up to date version for best functionality.

library(Seurat)
library(ArchR)
addArchRThreads(16)

######################################################
# Read RNA
######################################################

matRNA <- readRDS("Save-KNN-Groups-scRNA-Matrix.rds")

######################################################
# Read scRNA
######################################################

if(getwd() == "/scratch/users/jgranja/OtherAnalysis/GeneScores/Snakemake/PBMC"){
  seRNA <- readRDS("scRNA-pbmc_10k_v3.rds")
  seRNA <- seRNA[intersect(rownames(seRNA), rownames(matRNA)), ]
  seRNA <- NormalizeData(object = seRNA, verbose = TRUE)
  seRNA$Group <- seRNA$celltype
  RNA <- seRNA
}

if(getwd() == "/scratch/users/jgranja/OtherAnalysis/GeneScores/Snakemake/Heme"){
  seRNA <- readRDS("scRNA-Healthy-Hematopoiesis-191120.rds")
  seRNA <- seRNA[intersect(rownames(seRNA), rownames(matRNA)), ]
  seuratRNA <- CreateSeuratObject(counts = assay(seRNA))
  seuratRNA$Group <- colData(seRNA)[, "BioClassification", drop = TRUE]
  seuratRNA <- NormalizeData(object = seuratRNA, verbose = TRUE)
  RNA <- seuratRNA
  Idents(object=RNA) <- RNA$Group
}

######################################################
# Read Models
######################################################

path <- "output/ArchR/Results/Models"
files <- gtools::mixedsort(list.files(path))

#ATAC Models
matList <- lapply(seq_along(files), function(x){
  message(x)
  as.matrix(readRDS(file.path(path, files[x]))[rownames(matRNA), ])
}) %>% SimpleList
names(matList) <- gsub(".rds","",files)

if(getwd() == "/scratch/users/jgranja/OtherAnalysis/GeneScores/Snakemake/PBMC"){
  modelsCoA <- list.files("../../PBMC/KNN-Group-Matrices", pattern = "distCoA", full.names = TRUE)
}

if(getwd() == "/scratch/users/jgranja/OtherAnalysis/GeneScores/Snakemake/Heme"){
  modelsCoA <- list.files("../../Heme/KNN-Group-Matrices", pattern = "distCoA", full.names = TRUE)
}

matListCoA <- lapply(seq_along(modelsCoA), function(x){
  message(x)
  matx <- readRDS(modelsCoA[x])
  as.matrix(ArchR:::.safeSubset(as(matx, "dgCMatrix"), subsetRows = paste0(rownames(matRNA))))
}) %>% SimpleList
names(matListCoA) <- paste0("CoA-", seq_along(matListCoA))

matList <- c(matList, matListCoA)

matList$SnapATAC <- readRDS("Save-KNN-Groups-SnapATAC-Matrix.rds")
matList$Signac <- readRDS("Save-KNN-Groups-Signac-Matrix.rds")

matList <- matList[rev(seq_along(matList))]
saveRDS(matList, "Save-MatList.rds", compress = FALSE)


######################################################
# Save Image
######################################################
#save.image("Test-GeneScore-Models-March1.Rdata")

######################################################
# Load Image
######################################################
load("Test-GeneScore-Models-March1.Rdata")
matList <- readRDS("Save-MatList.rds")

######################################################
# Test Models in Variable Peaks
######################################################

RNA <- FindVariableFeatures(RNA, nfeatures = 2000)
varGenesAll <- VariableFeatures(RNA)[VariableFeatures(RNA) %in% rownames(matRNA)]

varGenes <- head(varGenesAll, 2000)

corVarGenes_GeneLvl <- lapply(seq_along(matList), function(x){
    
    #message(x)
    mx <- matList[[x]]
    mx <- t(t(mx) / colSums(mx)) * 10^4
    mx <- as.matrix(ArchR:::.safeSubset(as.matrix(data.frame(mx)), paste0(varGenes)))
    
    cor <- ArchR:::rowCorCpp(
      X = log2(mx + 1), 
      Y = matRNA[varGenes,], 
      idxX = seq_along(varGenes),
      idxY = seq_along(varGenes)
    )
    cor[is.na(cor)] <- 0
    cor[cor < 0] <- 0
    
    df <- data.frame(
      model = names(matList)[x], 
      cor = cor,
      gene = varGenes,
      geneRank = seq_along(varGenes)
    )

    df

}) %>% Reduce("rbind", .)

corVarGenes_SampleLvl <- lapply(seq_along(matList), function(x){

  mx <- matList[[x]]
  mx <- t(t(mx) / colSums(mx)) * 10^4
  mx <- as.matrix(ArchR:::.safeSubset(as.matrix(data.frame(mx)), paste0(varGenes)))
  mx <- log2(mx + 1)

  mz <- matRNA[varGenes,]

  corVals <- lapply(seq_len(ncol(matRNA)), function(z){
    cor(mx[, z], mz[, z])
  }) %>% unlist

  corVals[is.na(corVals)] <- 0
  corVals[corVals < 0] <- 0

  df <- data.frame(
    model = names(matList)[x], 
    cor = corVals,
    idx = seq_along(corVals)
  )
  
  df

}) %>% Reduce("rbind", .)

######################################################
# Test Models in Differential Peaks
######################################################

if(!file.exists("scRNA-Markers.rds")){
  library(future)
  plan(strategy = "multicore", workers = 10)
  RNAMarkers <- FindAllMarkers(RNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  saveRDS(RNAMarkers, "scRNA-Markers.rds")
}else{
  RNAMarkers <- readRDS("scRNA-Markers.rds")
}

library(dplyr)
i <- 1
varGenes <- {RNAMarkers %>% dplyr::group_by(cluster) %>% top_n(n = i, wt = avg_logFC)}$gene
while(length(varGenes) < 1000){
  i <- i + 1
  varGenes <- {RNAMarkers %>% dplyr::group_by(cluster) %>% top_n(n = i, wt = avg_logFC)}$gene
  varGenes <- unique(varGenes)
}

corDiffGenes_GeneLvl <- lapply(seq_along(matList), function(x){
    
    #message(x)
    mx <- matList[[x]]
    mx <- t(t(mx) / colSums(mx)) * 10^4
    mx <- as.matrix(ArchR:::.safeSubset(as.matrix(data.frame(mx)), paste0(varGenes)))
    
    cor <- ArchR:::rowCorCpp(
      X = log2(mx + 1), 
      Y = matRNA[varGenes,], 
      idxX = seq_along(varGenes),
      idxY = seq_along(varGenes)
    )
    cor[is.na(cor)] <- 0
    cor[cor < 0] <- 0
    
    df <- data.frame(
      model = names(matList)[x], 
      cor = cor,
      gene = varGenes,
      geneRank = seq_along(varGenes)
    )

    df

}) %>% Reduce("rbind", .)

corDiffGenes_SampleLvl <- lapply(seq_along(matList), function(x){

  mx <- matList[[x]]
  mx <- t(t(mx) / colSums(mx)) * 10^4
  mx <- as.matrix(ArchR:::.safeSubset(as.matrix(data.frame(mx)), paste0(varGenes)))
  mx <- log2(mx + 1)

  mz <- matRNA[varGenes,]

  corVals <- lapply(seq_len(ncol(matRNA)), function(z){
    cor(mx[, z], mz[, z])
  }) %>% unlist

  corVals[is.na(corVals)] <- 0
  corVals[corVals < 0] <- 0

  df <- data.frame(
    model = names(matList)[x], 
    cor = corVals,
    idx = seq_along(corVals)
  )
  
  df

}) %>% Reduce("rbind", .)

######################################################
# Plot Results
######################################################
modelRename <- data.frame(readr::read_tsv("ModelRename.tsv"))

palModel <- paletteDiscrete(values = gtools::mixedsort(unique(paste0(plotDF$model3))))

plotDF <- corDiffGenes_GeneLvl

plotDF$model2 <- mapLabels(
  labels = paste0(plotDF$model), 
  newLabels = paste0(modelRename[,2]), 
  oldLabels = paste0(modelRename[,1])
)
plotDF$model3 <- stringr::str_split(plotDF$model2, pattern = "-", simplify = TRUE)[,1]
dm <- stats::aggregate(plotDF$cor ~ plotDF$model2, FUN = median)
groupOrder <- paste0(dm[,1])[order(dm[,2], decreasing = TRUE)]
plotDF$model2 <- factor(paste0(plotDF$model2), groupOrder)
corMax <- max(dm[,2])

pdf(paste0(basename(getwd()),"-new-cor-diff-gene.pdf"), width = 10, height = 8)
ggplot(plotDF, aes(x = model2, y = cor, fill = model3)) +
  geom_violin(alpha = 1,  color = "black") +
  geom_boxplot(outlier.size = 0, outlier.stroke = 0, fill = NA, color = "black") +
  theme_ArchR(xText90 = TRUE) +
  scale_fill_manual(values = palModel) +
  ylab("Correlation") +
  geom_hline(yintercept = corMax, lty = "dashed", size = 0.5)
dev.off()


plotDF <- corDiffGenes_SampleLvl

plotDF$model2 <- mapLabels(
  labels = paste0(plotDF$model), 
  newLabels = paste0(modelRename[,2]), 
  oldLabels = paste0(modelRename[,1])
)
plotDF$model3 <- stringr::str_split(plotDF$model2, pattern = "-", simplify = TRUE)[,1]
dm <- stats::aggregate(plotDF$cor ~ plotDF$model2, FUN = median)
groupOrder <- paste0(dm[,1])[order(dm[,2], decreasing = TRUE)]
plotDF$model2 <- factor(paste0(plotDF$model2), groupOrder)
corMax <- max(dm[,2])

pdf(paste0(basename(getwd()),"-new-cor-diff-sample.pdf"), width = 10, height = 8)
ggplot(plotDF, aes(x = model2, y = cor, fill = model3)) +
  geom_violin(alpha = 1,  color = "black") +
  geom_boxplot(outlier.size = 0, outlier.stroke = 0, fill = NA, color = "black") +
  theme_ArchR(xText90 = TRUE) +
  scale_fill_manual(values = palModel) +
  ylab("Correlation") +
  geom_hline(yintercept = corMax, lty = "dashed", size = 0.5)
dev.off()



plotDF <- corVarGenes_GeneLvl

plotDF$model2 <- mapLabels(
  labels = paste0(plotDF$model), 
  newLabels = paste0(modelRename[,2]), 
  oldLabels = paste0(modelRename[,1])
)
plotDF$model3 <- stringr::str_split(plotDF$model2, pattern = "-", simplify = TRUE)[,1]
dm <- stats::aggregate(plotDF$cor ~ plotDF$model2, FUN = median)
groupOrder <- paste0(dm[,1])[order(dm[,2], decreasing = TRUE)]
plotDF$model2 <- factor(paste0(plotDF$model2), groupOrder)
corMax <- max(dm[,2])

pdf(paste0(basename(getwd()),"-new-cor-var-gene.pdf"), width = 10, height = 8)
ggplot(plotDF, aes(x = model2, y = cor, fill = model3)) +
  geom_violin(alpha = 1,  color = "black") +
  geom_boxplot(outlier.size = 0, outlier.stroke = 0, fill = NA, color = "black") +
  theme_ArchR(xText90 = TRUE) +
  scale_fill_manual(values = palModel) +
  ylab("Correlation") +
  geom_hline(yintercept = corMax, lty = "dashed", size = 0.5)
dev.off()


plotDF <- corVarGenes_SampleLvl

plotDF$model2 <- mapLabels(
  labels = paste0(plotDF$model), 
  newLabels = paste0(modelRename[,2]), 
  oldLabels = paste0(modelRename[,1])
)
plotDF$model3 <- stringr::str_split(plotDF$model2, pattern = "-", simplify = TRUE)[,1]
dm <- stats::aggregate(plotDF$cor ~ plotDF$model2, FUN = median)
groupOrder <- paste0(dm[,1])[order(dm[,2], decreasing = TRUE)]
plotDF$model2 <- factor(paste0(plotDF$model2), groupOrder)
corMax <- max(dm[,2])

pdf(paste0(basename(getwd()),"-new-cor-var-sample.pdf"), width = 10, height = 8)
ggplot(plotDF, aes(x = model2, y = cor, fill = model3)) +
  geom_violin(alpha = 1,  color = "black") +
  geom_boxplot(outlier.size = 0, outlier.stroke = 0, fill = NA, color = "black") +
  theme_ArchR(xText90 = TRUE) +
  scale_fill_manual(values = palModel) +
  ylab("Correlation") +
  geom_hline(yintercept = corMax, lty = "dashed", size = 0.5)
dev.off()


corDiffGenes_GeneLvl$model2 <- mapLabels(
  labels = paste0(corDiffGenes_GeneLvl$model), 
  newLabels = paste0(modelRename[,2]), 
  oldLabels = paste0(modelRename[,1])
)
diff_geneLvl <- stats::aggregate(corDiffGenes_GeneLvl$cor ~ corDiffGenes_GeneLvl$model2, FUN = median)
diff_geneLvl <- diff_geneLvl[order(diff_geneLvl[,2],decreasing=TRUE),]
diff_geneLvl$rank <- seq_len(nrow(diff_geneLvl))

corDiffGenes_SampleLvl$model2 <- mapLabels(
  labels = paste0(corDiffGenes_SampleLvl$model), 
  newLabels = paste0(modelRename[,2]), 
  oldLabels = paste0(modelRename[,1])
)
diff_sampleLvl <- stats::aggregate(corDiffGenes_SampleLvl$cor ~ corDiffGenes_SampleLvl$model2, FUN = median)
diff_sampleLvl <- diff_sampleLvl[order(diff_sampleLvl[,2],decreasing=TRUE),]
diff_sampleLvl$rank <- seq_len(nrow(diff_sampleLvl))

corVarGenes_GeneLvl$model2 <- mapLabels(
  labels = paste0(corVarGenes_GeneLvl$model), 
  newLabels = paste0(modelRename[,2]), 
  oldLabels = paste0(modelRename[,1])
)
var_geneLvl <- stats::aggregate(corVarGenes_GeneLvl$cor ~ corVarGenes_GeneLvl$model2, FUN = median)
var_geneLvl <- var_geneLvl[order(var_geneLvl[,2],decreasing=TRUE),]
var_geneLvl$rank <- seq_len(nrow(var_geneLvl))

corVarGenes_SampleLvl$model2 <- mapLabels(
  labels = paste0(corVarGenes_SampleLvl$model), 
  newLabels = paste0(modelRename[,2]), 
  oldLabels = paste0(modelRename[,1])
)
var_sampleLvl <- stats::aggregate(corVarGenes_SampleLvl$cor ~ corVarGenes_SampleLvl$model2, FUN = median)
var_sampleLvl <- diff_sampleLvl[order(diff_sampleLvl[,2],decreasing=TRUE),]
var_sampleLvl$rank <- seq_len(nrow(diff_sampleLvl))

rownames(var_sampleLvl) <- var_sampleLvl[,1]
rownames(var_geneLvl) <- var_geneLvl[,1]
rownames(diff_sampleLvl) <- diff_sampleLvl[,1]
rownames(diff_geneLvl) <- diff_geneLvl[,1]


df <- data.frame(
  row.names = var_sampleLvl[,1],
  model = var_sampleLvl[,1],
  Rank_Var_Samples = as.integer(var_sampleLvl[paste0(var_sampleLvl[,1]),"rank"]),
  Rank_Var_Genes = as.integer(var_geneLvl[paste0(var_sampleLvl[,1]),"rank"]),
  Rank_Diff_Samples = as.integer(diff_sampleLvl[paste0(var_sampleLvl[,1]),"rank"]),
  Rank_Diff_Genes = as.integer(diff_geneLvl[paste0(var_sampleLvl[,1]),"rank"])
)

df <- df[order(rowMeans(df[,2:ncol(df)]), decreasing = FALSE), ]
df2 <- df[,2:ncol(df)]
df2 <- apply(df2, 2, paste0)
df3 <- data.frame(row.names = rownames(df), model = stringr::str_split(rownames(df), pattern = "-", simplify=TRUE)[,1])

library(pheatmap)
pdf(paste0(basename(getwd()),"-New-Heatmap-Rank.pdf"), width = 10, height = 10)
pheatmap(df[,2:ncol(df)], 
  annotation_row = df3,
  annotation_colors = list(model = palModel),
  color = rev(paletteContinuous(set = "sambaNight")), 
  border_color = "black",
  number_color = "black",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  display_numbers = df2
)
dev.off()

