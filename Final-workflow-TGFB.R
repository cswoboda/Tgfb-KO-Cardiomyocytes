library(Seurat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(reticulate)
#reticulate::py_discover_config(required_module = "phate")
#reticulate::import("phate")
#library(phateR)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(stringr)
library(Rcpp)
library(SeuratObject)
library(Seurat)
library(Matrix)
library(dplyr)
library(cowplot)
library(ggplot2)
#library(harmony)
library(zellkonverter)
library(RcppCNPy)
##Initial data processing insert cellbender matrices here
##These can be obtained by request or at the count matrices of the objects posted to the GEO
options(future.globals.maxSize = 1000000000000000)
x <- ReadCB_h5("Control2_cellbender_corrected_FPR_0.1_filtered.h5")
control2 <- CreateSeuratObject(x, project = "Control2")
x <- ReadCB_h5("Control3_cellbender_corrected_filtered.h5")
control3 <- CreateSeuratObject(x, project = "Control3")
x <- ReadCB_h5("Mutant2_cellbender_corrected_FPR_0.1_filtered.h5")
mutant2 <- CreateSeuratObject(x, project = "Mutant2")
x <- ReadCB_h5("Mutant3_cellbender_corrected_FPR_0.1_filtered.h5")
mutant3 <- CreateSeuratObject(x, project = "Mutant3")

list <- c(control2, control3, mutant2, mutant3)

for (i in 1:length(list)) {
  list[[i]] <- PercentageFeatureSet(list[[i]], pattern = "^mt-", col.name = "percent.mt")
  z <- quantile(list[[i]]@meta.data[["nFeature_RNA"]], probs = c(0.05, 0.1, 0.15, 0.9, 0.98))
  y <- quantile(list[[i]]@meta.data[["nCount_RNA"]], probs = c(0.05, 0.1, 0.15, 0.9, 0.98))  
  #VlnPlot(list[[i]], features = c("nFeature_RNA"), pt.size = 0)
  #VlnPlot(list[[i]], features = c("nCount_RNA"), pt.size = 0)
  list[[i]] <-  subset(list[[i]], subset = nFeature_RNA <= z[[5]])
  list[[i]] <-  subset(list[[i]], subset = nCount_RNA <= y[[5]])
  list[[i]] <- subset(list[[i]], percent.mt < 20)
}

##Load in in doublet files.
names <- c("Control2", "Control3", "Mutant2", "Mutant3")
for (i in 1:length(list)) {
  fmat <- npyLoad(paste(paste("/Volumes/Millay-Sequencing/Atlas/solo-output", names[i], sep = "/"), "preds.npy", sep = "/"))
  character_fmat <- as.character(fmat)
  character_fmat <- replace(character_fmat, character_fmat == "0", "Singlet")
  character_fmat <- replace(character_fmat, character_fmat == "4.94065645841247e-324", "Doublet")
  list[[i]]$doublet <- character_fmat 
}

##Run processing for individual objects
for(i in 1:length(list)) {
  list[[i]] <- SCTransform(list[[i]])
  list[[i]] <- RunPCA(list[[i]])
  list[[i]] <- FindNeighbors(list[[i]], dims = 1:20)
  list[[i]] <- FindClusters(list[[i]], resolution = 1)
  list[[i]] <- RunUMAP(list[[i]], dims = 1:20)
}


obj <- list[[1]]
tmplist <- list[-1]
merged <- merge(obj, tmplist)

merged <- SCTransform(merged)
merged <- RunPCA(merged, npcs = 60, verbose = F)
merged<- IntegrateLayers(
  object = merged,
  method = CCAIntegration,
  normalization.method = "SCT",
  verbose = F,
  orig.reduction = "pca",
  new.reduction = "integrated.rpcasct60npc"
)


merged <- FindNeighbors(merged, dims = 1:60, reduction = "integrated.rpcasct60npc")
merged <- FindClusters(merged, resolution = 2, cluster.name = "cca-clusters-res2")
merged <- RunUMAP(merged, dims = 1:60, reduction = "integrated.rpcasct60npc", reduction.name = "umap60res2")
saveRDS(merged, "tmp-integrated-object-merged.RDS")

merged <- readRDS("tmp-integrated-object-merged.RDS")
merged <- FindNeighbors(merged, dims = 1:40, reduction = "integrated.rpcasct60npc")
merged <- FindClusters(merged, resolution = 2, cluster.name = "cca-clusters-res2")
merged <- RunUMAP(merged, dims = 1:40, reduction = "integrated.rpcasct60npc", reduction.name = "umap40res2")

merged <- FindNeighbors(merged, dims = 1:40, reduction = "integrated.rpcasct60npc")
merged <- FindClusters(merged, resolution = 1, cluster.name = "cca-clusters-res1")
merged <- RunUMAP(merged, dims = 1:40, reduction = "integrated.rpcasct60npc", reduction.name = "umap40res1")

merged <- FindNeighbors(merged, dims = 1:60, reduction = "integrated.rpcasct60npc")
merged <- FindClusters(merged, resolution = 0.5, cluster.name = "cca-clusters-res05")
merged <- RunUMAP(merged, dims = 1:60, reduction = "integrated.rpcasct60npc", reduction.name = "umap60res1")
##generate umaps
DimPlot(merged, reduction = "umap60res1")
DimPlot(merged, group.by = "doublet")

table(Idents(merged))
##Remove low quality cluster/doublet enriched cluster
merged <- subset(merged, idents = c("15"), invert = TRUE)

Idents(merged) <- "doublet"

##remove doublets, recluster
singlet_obj <- subset(merged, idents = c("Singlet"))


singlet_obj <- FindNeighbors(singlet_obj, dims = 1:60, reduction = "integrated.rpcasct60npc")
singlet_obj <- FindClusters(singlet_obj, resolution = 2, cluster.name = "cca-clusters-res2")
singlet_obj <- RunUMAP(singlet_obj, dims = 1:60, reduction = "integrated.rpcasct60npc", reduction.name = "umap60res2")

singlet_obj <- FindNeighbors(singlet_obj, dims = 1:40, reduction = "integrated.rpcasct60npc")
singlet_obj <- FindClusters(singlet_obj, resolution = 2, cluster.name = "cca-clusters-res2")
singlet_obj <- RunUMAP(singlet_obj, dims = 1:40, reduction = "integrated.rpcasct60npc", reduction.name = "umap40res2")

singlet_obj <- FindNeighbors(singlet_obj, dims = 1:40, reduction = "integrated.rpcasct60npc")
singlet_obj <- FindClusters(singlet_obj, resolution = 1, cluster.name = "cca-clusters-res1")
singlet_obj <- RunUMAP(singlet_obj, dims = 1:40, reduction = "integrated.rpcasct60npc", reduction.name = "umap40res1")

singlet_obj <- FindNeighbors(singlet_obj, dims = 1:60, reduction = "integrated.rpcasct60npc")
singlet_obj <- FindClusters(singlet_obj, resolution = 0.5, cluster.name = "cca-clusters-res05")
singlet_obj <- RunUMAP(singlet_obj, dims = 1:60, reduction = "integrated.rpcasct60npc", reduction.name = "umap60res1")
DimPlot(singlet_obj, reduction = "umap60res1", label =  TRUE)

DefaultAssay(singlet_obj) <- "RNA"
singlet_obj <- NormalizeData(singlet_obj)
singlet_obj <- FindVariableFeatures(singlet_obj)
singlet_obj <- ScaleData(singlet_obj)

singlet_obj <- JoinLayers(singlet_obj, assay = "RNA")
##genrate arkers
markers <- FindAllMarkers(singlet_obj, assay = "RNA", only.pos = TRUE)
write.csv(markers, "all-markers-integrated.csv")

##Feature plots for marker genes
FeaturePlot(singlet_obj, features = c("Dcn"))
ggsave("DCN-Feature-Plot.jpeg", dpi = 1000, width = 10, height = 10)
FeaturePlot(singlet_obj, features = c("Ttn"))
ggsave("TTN-Feature-Plot.jpeg", dpi = 1000, width = 10, height = 10)
FeaturePlot(singlet_obj, features = c("Pecam1"))
ggsave("Pecam1-Feature-Plot.jpeg", dpi = 1000, width = 10, height = 10)
FeaturePlot(singlet_obj, features = c("Mpz"))
ggsave("Mpz-Feature-Plot.jpeg", dpi = 1000, width = 10, height = 10)
FeaturePlot(singlet_obj, features = c("Myh11"))
ggsave("Myh11-Feature-Plot.jpeg", dpi = 1000, width = 10, height = 10)
FeaturePlot(singlet_obj, features = c("Mkx"))
ggsave("Mkx-Feature-Plot.jpeg", dpi = 1000, width = 10, height = 10)
FeaturePlot(singlet_obj, features = c("Ptprc"))
ggsave("Ptprc-Feature-Plot.jpeg", dpi = 1000, width = 10, height = 10)
singlet_obj$rawcluster <- Idents(singlet_obj)
Idents(singlet_obj) <- "rawcluster"


DimPlot(singlet_obj, label = TRUE, group.by = "rawcluster") + NoLegend()
Idents(singlet_obj) <- "rawcluster"
singlet_obj <- RenameIdents(singlet_obj, "0" = "Cardiomyocytes", "1" = "Cardiomyocytes", "4" = "Cardiomyocytes", "5" = "Cardiomyocytes", "9" = "Cardiomyocytes", "11" = "Cardiomyocytes", "7" = "Cardiomyocytes", "18" = "Cardiomyocytes", "15" = "Cardiomyocytes", "8" = "FAPs", "2" = "FAPs", "12" = "FAPs")
##Labels for Celltype
singlet_obj$celltypesimple <- Idents(singlet_obj)
Idents(singlet_obj) <- "rawcluster"
##Labels for subclusters
singlet_obj <- RenameIdents(singlet_obj, "0" = "CM1", "1" = "CM2", "4" = "CM3", "5" = "CM4", "9" = "CM5", "11" = "CM6", "7" = "CM8", "18" = "CM9", "15" = "CM10", "8" = "FAP1", "2" = "FAP2", "12" = "FAP3")
singlet_obj$subcelltype <- Idents(singlet_obj)

##Create subset objects for plotting genes
Idents(singlet_obj) <- "celltypesimple"
cm_subset <- subset(singlet_obj, idents = c("Cardiomyocytes"))
Idents(singlet_obj) <- "celltypesimple"
faps_subset <- subset(singlet_obj, idents = c("FAPs"))

#DefaultAssay(cm_subset) <- "SCT"
cm_subset$sub.orig <- paste(cm_subset$subcelltype, cm_subset$condition, sep = "_")
DefaultAssay(cm_subset) <- "RNA"
x <- FindMarkers(cm_subset, ident.1 = c("CM1_Control"), ident.2 = c("CM1_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "CM1_UP_CONTROL.csv")
x <- FindMarkers(cm_subset, ident.1 = c("CM2_Control"), ident.2 = c("CM2_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "CM2_UP_CONTROL.csv")
x <- FindMarkers(cm_subset, ident.1 = c("CM3_Control"), ident.2 = c("CM3_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "CM3_UP_CONTROL.csv")
x <- FindMarkers(cm_subset, ident.1 = c("CM4_Control"), ident.2 = c("CM4_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "CM4_UP_CONTROL.csv")
x <- FindMarkers(cm_subset, ident.1 = c("CM5_Control"), ident.2 = c("CM5_Mutant"), group.by = "sub.orig")
write.csv(x, "CM5_UP_CONTROL.csv")



x <- FindMarkers(cm_subset, ident.1 = c("CM6_Control"), ident.2 = c("CM6_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "CM6_UP_CONTROL.csv")
x <- FindMarkers(cm_subset, ident.1 = c("CM8_Control"), ident.2 = c("CM8_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "CM8_UP_CONTROL.csv")
x <- FindMarkers(cm_subset, ident.1 = c("CM9_Control"), ident.2 = c("CM9_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "CM9_UP_CONTROL.csv")

##Cardiomyocyte gene expression
VlnPlot(cm_subset, features = c("Myh7"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-myh7-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)
VlnPlot(cm_subset, features = c("Myh6"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-myh6-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)
VlnPlot(cm_subset, features = c("Tnni1"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-tnni1-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)
VlnPlot(cm_subset, features = c("Tnni3"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-tnni3-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)
VlnPlot(cm_subset, features = c("Myl7"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-myl7-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)
VlnPlot(cm_subset, features = c("Myl2"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-myl2-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)
VlnPlot(cm_subset, features = c("Myl4"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-myl4-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)
VlnPlot(cm_subset, features = c("Myl3"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-myl3-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)
##Ligand and Receptor gene expression:

VlnPlot(cm_subset, features = c("Lamc1"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-lamc1-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)
VlnPlot(cm_subset, features = c("Lamb1"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-lamb1-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)
VlnPlot(cm_subset, features = c("Lamb2"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-lamb2-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)
VlnPlot(cm_subset, features = c("Itgb1"), group.by = "celltypesimple", split.by = "condition", pt.size = 0) + NoLegend()
ggsave("vln-Itgb1-cm-grouped-only.pdf", dpi = 1000, width = 10, height = 10)


DefaultAssay(faps_subset) <- "SCT"
faps_subset$sub.orig <- paste(faps_subset$subcelltype, faps_subset$condition, sep = "_")
DefaultAssay(faps_subset) <- "RNA"
x <- FindMarkers(faps_subset, ident.1 = c("FAP1_Control"), ident.2 = c("FAP1_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "FAP1_UP_CONTROL.csv")
x <- FindMarkers(faps_subset, ident.1 = c("FAP2_Control"), ident.2 = c("FAP2_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "FAP2_UP_CONTROL.csv")
x <- FindMarkers(faps_subset, ident.1 = c("FAP3_Control"), ident.2 = c("FAP3_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "FAP3_UP_CONTROL.csv")

x <- FindMarkers(faps_subset, ident.2 = c("FAP1_Control"), ident.1 = c("FAP1_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "FAP1_DOWN_CONTROL.csv")
x <- FindMarkers(faps_subset, ident.2 = c("FAP2_Control"), ident.1 = c("FAP2_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "FAP2_DOWN_CONTROL.csv")
x <- FindMarkers(faps_subset, ident.2 = c("FAP3_Control"), ident.1 = c("FAP3_Mutant"), only.pos = TRUE, group.by = "sub.orig")
write.csv(x, "FAP3_DOWN_CONTROL.csv")

cm_subset_markers <- FindAllMarkers(cm_subset)
DefaultAssay(cm_subset) <- "RNA"



idents <- singlet_obj$orig.ident
x <- which(idents %in% c("Mutant2", "Mutant3"))
idents[x] <- "Mutant"
idents[-x] <- "Control"

singlet_obj$condition <- idents

singlet_obj$celltype_condition <- paste(singlet_obj$celltypesimple, singlet_obj$condition, sep = "_")
y <- FindMarkers(singlet_obj, group.by = "celltype_condition", ident.1 = "Cardiomyocytes_Mutant", ident.2 = "Cardiomyocytes_Control", only.pos = TRUE)
write.csv(y, "UP-in-mutant-cm.csv")
y <- FindMarkers(singlet_obj, group.by = "celltype_condition", ident.2 = "Cardiomyocytes_Mutant", ident.1 = "Cardiomyocytes_Control", only.pos = TRUE)
write.csv(y, "UP-in-control-cm.csv")

y <- FindMarkers(singlet_obj, group.by = "celltype_condition", ident.1 = "FAPs_Mutant", ident.2 = "FAPs_Control", only.pos = TRUE)
write.csv(y, "UP-in-mutant-fibro.csv")
y <- FindMarkers(singlet_obj, group.by = "celltype_condition", ident.2 = "FAPs_Mutant", ident.1 = "FAPs_Control", only.pos = TRUE)
write.csv(y, "UP-in-control-fibro.csv")

##Run metabolic scoring
metabolic <- c("Atp5c1", "Atp5f1", "Ndufab1", "Ndufa12", "Ndufv2", "Cox6a1"   
               "Hspd1", "Higd1a", "Ndufb5", "Acadvl", "Cox8a", "Ech1", "Atp5g2", "Atp5h" "Ndufa8"    "Cycs"      "Ndufb3"    "Ndufa7"   
               "Mtch1", "Ndufb1-ps", "Ndufb2","Cox6b1","Ndufs2","Cox5a","Ndufc2","Ndufv1","Idh3g"     "Acads"     "Cox7c"     "Ndufs7"   
               "Tomm20",  "Atp5j","Bax")

metabolic_list <- list()
metabolic_list[[1]] <- metabolic
singlet_obj <- AddModuleScore(singlet_obj, features = metabolic_list, name = "metabolic-score")

##subset for cellchat
Idents(singlet_obj) <- "celltypesimple"
tmp <- subset(singlet_obj, idents = c("Cardiomyocytes", "FAPs"))
Idents(tmp) <- "condition"

mutant <- subset(tmp, idents = "Mutant")
control <- subset(tmp, idents = "Control")



##generate cellchat objects
library(CellChat)
cellChat_mutant <- createCellChat(object = mutant, group.by = "celltypesimple", assay = "RNA")
cellChat_ctrl <- createCellChat(object = control, group.by = "celltypesimple", assay = "RNA")
CellChatDB.use <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
#CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
##loccalize to ECM interactions
showDatabaseCategory(CellChatDB.use)
CellChatDB.use <- subsetDB(CellChatDB.use, search = c("ECM-Receptor"), key = "annotation") # use Secreted Signaling
cellChat_mutant@DB <- CellChatDB.use
cellChat_ctrl@DB <- CellChatDB.use
cellChat_mutant@meta[["subcelltype"]] <- droplevels(cellChat_mutant@meta[["subcelltype"]])
cellChat_ctrl@meta[["subcelltype"]] <- droplevels(cellChat_ctrl@meta[["subcelltype"]])
showDatabaseCategory(CellChatDB.use)
CellChatDB.use <- subsetDB(CellChatDB.use, search = c("ECM-Receptor"), key = "annotation") # use Secreted Signaling
cellChat_mutant@DB <- CellChatDB.use
cellChat_ctrl@DB <- CellChatDB.use
##Drop levels for any identities not chared between the two
cellChat_mutant@meta[["celltypesimple"]] <- droplevels(cellChat_mutant@meta[["celltypesimple"]])
cellChat_ctrl@meta[["celltypesimple"]] <- droplevels(cellChat_ctrl@meta[["celltypesimple"]])
cellChat_ctrl@idents <- droplevels(cellChat_ctrl@idents)
cellChat_mutant@idents <- droplevels(cellChat_mutant@idents)
cellChat_ctrl@idents <- droplevels(cellChat_ctrl@idents)
cellChat_mutant@idents <- droplevels(cellChat_mutant@idents)

#cellChat_mutant$ v= droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))


cellChat_mutant <- subsetData(cellChat_mutant) # This step is necessary even if using the whole database
cellChat_mutant <- identifyOverExpressedGenes(cellChat_mutant)
cellChat_mutant <- identifyOverExpressedInteractions(cellChat_mutant)
cellChat_ctrl <- subsetData(cellChat_ctrl) # This step is necessary even if using the whole database
cellChat_ctrl <- identifyOverExpressedGenes(cellChat_ctrl)
cellChat_ctrl <- identifyOverExpressedInteractions(cellChat_ctrl)
cellChat_ctrl <- computeCommunProb(cellChat_ctrl, type = "triMean")
cellChat_mutant <- computeCommunProb(cellChat_mutant, type = "triMean")
cellChat_ctrl <- filterCommunication(cellChat_ctrl, min.cells = 10)
cellChat_mutant <- filterCommunication(cellChat_mutant, min.cells = 10)
cellChat_ctrl <- computeCommunProbPathway(cellChat_ctrl)
cellChat_mutant <- computeCommunProbPathway(cellChat_mutant)
cellChat_mutant <- aggregateNet(cellChat_mutant)
cellChat_ctrl <- aggregateNet(cellChat_ctrl)
cellChat_mutant <- netAnalysis_computeCentrality(cellChat_mutant, slot.name = "netP")
cellChat_ctrl <- netAnalysis_computeCentrality(cellChat_ctrl, slot.name = "netP")
object.list <- list(CTRL = cellChat_ctrl, MT = cellChat_mutant)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, targets.use = c(1), sources.use = c(2),  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in MT", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, targets.use = c(1:8), sources.use = c(10:13),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in MT", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "CTRL"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "CTRL",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "MT",ligand.logFC = -0.05, receptor.logFC = NULL)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairtmp, sources.use = c(1:8), targets.use = c(10:13), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:8), targets.use = c(10:13), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2


pathways.show <- c("LAMININ") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

vector <- pairLR.use.up$interaction_name
first_five_chars <- substr(vector, 1, 5)
x <- which(first_five_chars %in% "LAMC1")
pairtmp <- as.data.frame(pairLR.use.up[x, ])
pairtmp$interaction_name <- pairtmp$`pairLR.use.up[x, ]`

pairtmp$`pairLR.use.up[x, ]` <- NULL
#pairLR.use.up = net.up[, "interaction_name", drop = F]


pathways.show <- c("AGRN") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("CTRL", "MT")) # set factor level
plotGeneExpression(cellchat, signaling = "AGRN", split.by = "datasets", colors.ggplot = T, type = "violin")


cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("CTRL", "MT")) # set factor level
plotGeneExpression(cellchat, signaling = "COLLAGEN", split.by = "datasets", colors.ggplot = T, type = "violin")


cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("CTRL", "MT")) # set factor level
plotGeneExpression(cellchat, signaling = "LAMININ", split.by = "datasets", colors.ggplot = T, type = "violin")






##n the colorbar, red represents increased signaling in the second dataset compared to the first one.
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = c("FAP1", "FAP2", "FAP3"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = c("FAP1", "FAP2", "FAP3"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))



num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

df.net <- subsetCommunication(cellChat_ctrl, sources.use = c(1,2,3,4,5,6,7,8,9), targets.use = c(10,11,12))
write.csv(df.net, "cardiomyocytes to faps cellchat control.csv")
df.net <- subsetCommunication(cellChat_ctrl, targets.use = c(1,2,3,4,5,6,7,8,9), sources.use = c(10,11,12))
write.csv(df.net, "faps to cardiomyocytes cellchat control.csv")
df.net <- subsetCommunication(cellChat_mutant, sources.use = c(1,2,3,4,5,6,7,8,9), targets.use = c(10,11,12))
write.csv(df.net, "cardiomyocytes to faps cellchat mutant.csv")
df.net <- subsetCommunication(cellChat_mutant, targets.use = c(1,2,3,4,5,6,7,8,9), sources.use = c(10,11,12))
write.csv(df.net,"faps to cardiomyocytes cellchat mutant.csv")

pathways.show <- c("VEGF") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
pathways.show <- c("TGFb") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
pathways.show <- c("HSPG") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
pathways.show <- c("BMP") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("CTRL", "MT")) # set factor level
plotGeneExpression(cellchat, signaling = "TGFb", split.by = "datasets", colors.ggplot = T, type = "violin")
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("CTRL", "MT")) # set factor level
plotGeneExpression(cellchat, signaling = "BMP", split.by = "datasets", colors.ggplot = T, type = "violin")
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("CTRL", "MT")) # set factor level
plotGeneExpression(cellchat, signaling = "HSPG", split.by = "datasets", colors.ggplot = T, type = "violin")
pos.dataset = "MT"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "MT",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "CTRL",ligand.logFC = -0.05, receptor.logFC = NULL)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 1, targets.use = c(3), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 1, targets.use = c(3), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2


pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 3, targets.use = c(1), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 3, targets.use = c(1), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

