### Vcl paper Transcriptomic code

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### R script
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

store.path <- c("//spatial_analysis/seurat/merge_result/test_YFP/")
store.path <- c("//spatial_analysis/seurat/merge_result/cell-cell/edit/final/")

plot_name <- c("merge")

saveplot <- function(plot_name, plot, suffix, store.path, width, height) {
    filename <-paste(plot_name,"_",suffix,".pdf",sep = "")
	ggsave(filename, plot, path = store.path, width = width, height = height,dpi=800)
}

options(future.globals.maxSize = 100000 * 1024^2) ###Set CPU


### Pre-process the each sample
### 4 sample A&B control C&D Vcl cKO
### In the paper, the sample represent A1~C1; B1~C2, C1~M1,D1~M2
data.dirA <- c("//YFP_data/A1")
data.dirB <- c("//YFP_data/B1")
data.dirC <- c("//YFP_data/C1")
data.dirD <- c("//YFP_data/D1")

seurat_obj_A <- Load10X_Spatial(data.dir = data.dirA, slice = "A1")
seurat_obj_B <- Load10X_Spatial(data.dir = data.dirB, slice = "B1")
seurat_obj_C <- Load10X_Spatial(data.dir = data.dirC, slice = "C1")
seurat_obj_D <- Load10X_Spatial(data.dir = data.dirD, slice = "D1")

seurat_obj_A$orig.ident <- "A1"
seurat_obj_B$orig.ident <- "B1"
seurat_obj_C$orig.ident <- "C1"
seurat_obj_D$orig.ident <- "D1"


### Quality Control
### Figure S6A
ABCD.merge <- merge(seurat_obj_A, y = c(seurat_obj_B,seurat_obj_C,seurat_obj_D))

### nConunt&nFeature in Spatail vision
SpatialSample <- SpatialFeaturePlot(ABCD.merge, features = c("nCount_Spatial", "nFeature_Spatial"))
saveplot(plot_name,SpatialSample,"nCount&nFeature_SpatialPlot",store.path, 18, 12)

### nConunt&nFeature in Vlnplot
VlnSample <- VlnPlot(ABCD.merge, group.by = "orig.ident", features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2, cols = c("#4DAA48","#4DAA48","#FF0000","#FF0000")) + NoLegend()
saveplot(plot_name,VlnSample,"nCount&nFeature_VlnPlot",store.path, 9,6)


### SCTransform Normalization and Clustering for each sample
seurat_obj_A <- SCTransform(seurat_obj_A, assay = "Spatial", verbose = FALSE)
seurat_obj_A <- seurat_obj_A[, seurat_obj_A$nCount_Spatial > 10000 & seurat_obj_A$nFeature_Spatial > 4000]
seurat_obj_A  <- RunPCA(seurat_obj_A , assay = "SCT", verbose = F)
seurat_obj_A  <- FindNeighbors(seurat_obj_A , reduction = "pca", dims = 1:6, verbose = F)
seurat_obj_A  <- FindClusters(seurat_obj_A , verbose = F)
seurat_obj_A  <- RunUMAP(seurat_obj_A , dims = 1:6, verbose = F)

seurat_obj_B <- SCTransform(seurat_obj_B, assay = "Spatial", verbose = FALSE)
seurat_obj_B <- seurat_obj_B[, seurat_obj_B$nCount_Spatial > 10000 & seurat_obj_B$nFeature_Spatial > 4000]
seurat_obj_B  <- RunPCA(seurat_obj_B , assay = "SCT", verbose = F)
seurat_obj_B  <- FindNeighbors(seurat_obj_B , reduction = "pca", dims = 1:6, verbose = F)
seurat_obj_B  <- FindClusters(seurat_obj_B , verbose = F)
seurat_obj_B  <- RunUMAP(seurat_obj_B , dims = 1:6, verbose = F)

seurat_obj_C <- SCTransform(seurat_obj_C, assay = "Spatial", verbose = FALSE)
seurat_obj_C <- seurat_obj_C[, seurat_obj_C$nCount_Spatial > 10000 & seurat_obj_C$nFeature_Spatial > 4000]
seurat_obj_C  <- RunPCA(seurat_obj_C , assay = "SCT", verbose = F)
seurat_obj_C  <- FindNeighbors(seurat_obj_C , reduction = "pca", dims = 1:6, verbose = F)
seurat_obj_C  <- FindClusters(seurat_obj_C , verbose = F)
seurat_obj_C  <- RunUMAP(seurat_obj_C , dims = 1:6, verbose = F)

seurat_obj_D <- SCTransform(seurat_obj_D, assay = "Spatial", verbose = FALSE)
seurat_obj_D <- seurat_obj_D[, seurat_obj_D$nCount_Spatial > 10000 & seurat_obj_D$nFeature_Spatial > 4000]
seurat_obj_D  <- RunPCA(seurat_obj_D , assay = "SCT", verbose = F)
seurat_obj_D  <- FindNeighbors(seurat_obj_D , reduction = "pca", dims = 1:6, verbose = F)
seurat_obj_D  <- FindClusters(seurat_obj_D , verbose = F)
seurat_obj_D  <- RunUMAP(seurat_obj_D , dims = 1:6, verbose = F)


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### Anchor_based Integration of 4 spatial sample s
# Run SCTransform integration workflow: https://satijalab.org/seurat/archive/v3.0/integration.html

ABCD.list <- list(seurat_obj_A,seurat_obj_B,seurat_obj_C,seurat_obj_D)
genes.common <- Reduce(intersect, list(rownames(seurat_obj_A), rownames(seurat_obj_B), rownames(seurat_obj_C), rownames(seurat_obj_D)))
list.features <- SelectIntegrationFeatures(object.list = ABCD.list,
                                           nfeatures = 2000,
                                           assay = c("SCT", "SCT", "SCT", "SCT"))
ABCD.list <- PrepSCTIntegration(object.list = ABCD.list,
                               anchor.features = list.features,
                               assay = "SCT",
                               verbose = F)
ABCD.anchors <- FindIntegrationAnchors(object.list = ABCD.list,
                                      normalization.method = "SCT",
                                      anchor.features = list.features,
                                      verbose = F)                       
ABCD <- IntegrateData(anchorset = ABCD.anchors,
                     features.to.integrate = genes.common,
                     normalization.method = "SCT", 
                     verbose = F)

saveRDS(ABCD,"//spatial_analysis/seurat/merge_result/cell-cell/edit/final/ABCD.rds")
ABCD <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/edit/final/ABCD.rds")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### Rerun dimensionality reduction and clustering on integrated object.

DefaultAssay(ABCD) <- "integrated"
ABCD <- RunPCA(ABCD, verbose = F)
ABCD <- FindNeighbors(ABCD, reduction = "pca", dims = 1:30)
ABCD <- FindClusters(ABCD, resolution = 12, verbose = F)
ABCD <- RunUMAP(ABCD, reduction = "pca", dims = 1:30)


ABCD.subset <- subset(ABCD, idents = c(3,56,96,22,43,28,72,53,30,45)) #RES = 12

SpatialSample <- SpatialDimPlot(ABCD, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 4, ncol = 2)
saveplot(plot_name,SpatialSample,"SpatialDimPlot_ABCD.subset_original",store.path,24,24)

p1 <- DimPlot(ABCD.subset, reduction = "umap", label =TRUE, label.size = 4, label.box = TRUE)
saveplot(plot_name,p1,"UMAP_ABCD.subset_orignal",store.path,24,24)

### remove the spots which are distinct to the heart region on H&E image
ABCD.subset <- subset(ABCD.subset, a1_imagecol>210 ,invert = TRUE)
ABCD.subset <- subset(ABCD.subset, b1_imagecol>210 ,invert = TRUE)
ABCD.subset <- subset(ABCD.subset, c1_imagecol>280&c1_imagerow>200 ,invert = TRUE)
ABCD.subset <- subset(ABCD.subset, c1_imagecol>230&c1_imagerow<200 ,invert = TRUE)
ABCD.subset <- subset(ABCD.subset, d1_imagecol>280&d1_imagerow>200 ,invert = TRUE)

### remove the spots which are distinct to main region in UMAP
ABCD.subset <- ABCD.subset[,ABCD.subset@reductions$umap@cell.embeddings[,2] > -3]

### Draw the SpatialDimplot
Idents(ABCD.subset) <- "integrated_snn_res.12"
SpatialSample <- SpatialDimPlot(ABCD.subset, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 4, ncol = 2, cols =c(
"3" = "#FF6C67","30" = "#E48C00","53" = "#A3A500" ,"56" = "#00B900" ,"96" = "#00C377","22" = "#00C3C6" ,"28" = "#00B3FD","43" = "#9690FF" ,"45" = "#F961FA" ,"72" = "#FF52BF" ))
saveplot(plot_name,SpatialSample,"SpatialDimPlot_ABCD_original",store.path,24,24)

Sample <- subset(ABCD.subset,idents = c(3,30,53,56,96))
SpatialSample <- SpatialDimPlot(Sample,crop = FALSE, label = FALSE, pt.size.factor = 1.2, ncol = 4)
saveplot(plot_name,SpatialSample,"SpatialDimPlot_OFT",store.path,24,24)


saveRDS(ABCD.subset,"//spatial_analysis/seurat/merge_result/cell-cell/edit/final/ABCD.subset1.rds")


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### OFT and non_OFT annotation
ABCD.subset@meta.data$OFT <- "Null"

ABCD.subset@meta.data[ABCD.subset@meta.data$integrated_snn_res.12=="3",]$OFT <- "OFT_region"
ABCD.subset@meta.data[ABCD.subset@meta.data$integrated_snn_res.12=="56",]$OFT <- "OFT_region"
ABCD.subset@meta.data[ABCD.subset@meta.data$integrated_snn_res.12=="96",]$OFT <- "OFT_region"
ABCD.subset@meta.data[ABCD.subset@meta.data$integrated_snn_res.12=="53",]$OFT <- "OFT_region"
ABCD.subset@meta.data[ABCD.subset@meta.data$integrated_snn_res.12=="30",]$OFT <- "OFT_region"

ABCD.subset@meta.data[ABCD.subset@meta.data$integrated_snn_res.12=="22",]$OFT <- "non_OFT_region"
ABCD.subset@meta.data[ABCD.subset@meta.data$integrated_snn_res.12=="43",]$OFT <- "non_OFT_region"
ABCD.subset@meta.data[ABCD.subset@meta.data$integrated_snn_res.12=="28",]$OFT <- "non_OFT_region"
ABCD.subset@meta.data[ABCD.subset@meta.data$integrated_snn_res.12=="72",]$OFT <- "non_OFT_region"
ABCD.subset@meta.data[ABCD.subset@meta.data$integrated_snn_res.12=="45",]$OFT <- "non_OFT_region"


### Draw the SpatialDimplot
Idents(ABCD.subset) <- "OFT"
SpatialSample <- SpatialDimPlot(ABCD.subset, crop = FALSE, label = FALSE, pt.size.factor = 1, label.size = 4, ncol = 2, cols=c("OFT_region"="#A9D18E","non_OFT_region"="#FF9999"))
saveplot(plot_name,SpatialSample,"SpatialDimPlot_OFT&non_OFT",store.path,10,10)


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
Idents(ABCD.subset) <- "integrated_snn_res.12"
levels(ABCD.subset)<-c(3,30,53,56,96,22,28,43,45,72)
DefaultAssay(ABCD.subset) <- "SCT"
p1 <- VlnPlot(ABCD.subset,features = c("YFP"))
saveplot(plot_name,p1,"VlnPlot_YFP_SCT",store.path,5,5)

DefaultAssay(ABCD.subset) <- "Spatial"
p1 <- VlnPlot(ABCD.subset,features = c("YFP"))
saveplot(plot_name,p1,"VlnPlot_YFP_Spatial",store.path,5,5)

DefaultAssay(ABCD.subset) <- "integrated"
p1 <- VlnPlot(ABCD.subset,features = c("YFP"))
saveplot(plot_name,p1,"VlnPlot_YFP_integrated",store.path,5,5)


Feature_MES <- list(c("Lum", "Pdgfra", "Sox9"))
Feature_VSMC <- list(c("Myh11", "Cxcl12"))
Feature_Neurons <- list(c("Phox2b", "Th"))
Feature_Schwann <- list(c("Sox10", "Fabp7"))

Idents(ABCD.subset) <- "integrated_snn_res.12"
levels(ABCD.subset)<-c(3,30,53,56,96,22,28,43,45,72)
DefaultAssay(ABCD.subset) <- "SCT"

ABCD.subset <- AddModuleScore(object = ABCD.subset, features = Feature_MES, name = "Feature_MES")
ABCD.subset <- AddModuleScore(object = ABCD.subset, features = Feature_VSMC, name = "Feature_VSMC")
ABCD.subset <- AddModuleScore(object = ABCD.subset, features = Feature_Neurons, name = "Feature_Neurons")
ABCD.subset <- AddModuleScore(object = ABCD.subset, features = Feature_Schwann, name = "Feature_Schwann")

Sample <- VlnPlot(ABCD.subset,features = c("Feature_MES1","Feature_VSMC1","Feature_Neurons1","Feature_Schwann1"),ncol = 1)
saveplot(plot_name,Sample,"AddModuleScore_VlnPlot_Feature_SCT",store.path,5,16)


Sample <- VlnPlot(ABCD.subset,features = c("Feature_MES1"),ncol = 1,cols =c(
"3" = "#FF6C67","30" = "#E48C00","53" = "#A3A500" ,"56" = "#00B900" ,"96" = "#00C377","22" = "#00C3C6" ,"28" = "#00B3FD","43" = "#9690FF" ,"45" = "#F961FA" ,"72" = "#FF52BF" ) ) + ylab("Score") + labs(title ="AddModuleScore of MES markers")
saveplot(plot_name,Sample,"AddModuleScore_VlnPlot_Feature_SCT_MES",store.path,6,5)

Sample <- VlnPlot(ABCD.subset,features = c("Feature_VSMC1"),ncol = 1,cols =c(
"3" = "#FF6C67","30" = "#E48C00","53" = "#A3A500" ,"56" = "#00B900" ,"96" = "#00C377","22" = "#00C3C6" ,"28" = "#00B3FD","43" = "#9690FF" ,"45" = "#F961FA" ,"72" = "#FF52BF" ) ) + ylab("Score") + labs(title ="AddModuleScore of VSMC markers")
saveplot(plot_name,Sample,"AddModuleScore_VlnPlot_Feature_SCT_VSMC",store.path,6,5)


Sample <- VlnPlot(ABCD.subset,features = c("Feature_Neurons1"),ncol = 1,cols =c(
"3" = "#FF6C67","30" = "#E48C00","53" = "#A3A500" ,"56" = "#00B900" ,"96" = "#00C377","22" = "#00C3C6" ,"28" = "#00B3FD","43" = "#9690FF" ,"45" = "#F961FA" ,"72" = "#FF52BF" ) ) + ylab("Score") + labs(title ="AddModuleScore of Neurons markers")
saveplot(plot_name,Sample,"AddModuleScore_VlnPlot_Feature_SCT_Neurons",store.path,6,5)


Sample <- VlnPlot(ABCD.subset,features = c("Feature_Schwann1"),ncol = 1,cols =c(
"3" = "#FF6C67","30" = "#E48C00","53" = "#A3A500" ,"56" = "#00B900" ,"96" = "#00C377","22" = "#00C3C6" ,"28" = "#00B3FD","43" = "#9690FF" ,"45" = "#F961FA" ,"72" = "#FF52BF" ) ) + ylab("Score") + labs(title ="AddModuleScore of Schwann markers")
saveplot(plot_name,Sample,"AddModuleScore_VlnPlot_Feature_SCT_Schwann",store.path,6,5)


# OFT_MES <- list(c("Lum", "Pdgfra", "Sox9"))
# OFT_VSMC <- list(c("Myh11", "Cxcl12"))
# OFT_Neurons <- list(c("Phox2a", "Th"))
# OFT_Schwann <- list(c("Sox10", "Fabp7"))

# Idents(ABCD.subset) <- "integrated_snn_res.12"
# levels(ABCD.subset)<-c(3,30,53,56,96,22,28,43,45,72)
# DefaultAssay(ABCD.subset) <- "integrated"

# ABCD.subset <- AddModuleScore(object = ABCD.subset, features = OFT_MES, name = "OFT_MES")
# ABCD.subset <- AddModuleScore(object = ABCD.subset, features = OFT_VSMC, name = "OFT_VSMC")
# ABCD.subset <- AddModuleScore(object = ABCD.subset, features = OFT_Neurons, name = "OFT_Neurons")
# ABCD.subset <- AddModuleScore(object = ABCD.subset, features = OFT_Schwann, name = "OFT_Schwann")

# Sample <- VlnPlot(ABCD.subset,features = c("OFT_MES1","OFT_VSMC1","OFT_Neurons1","OFT_Schwann1"),ncol = 1)
# saveplot(plot_name,Sample,"AddModuleScore_VlnPlot_OFT_integrated",store.path,5,16)


saveRDS(ABCD.subset,"//spatial_analysis/seurat/merge_result/cell-cell/edit/final/ABCD.subset2.rds")


### OTF region subset
Idents(ABCD.subset) <- "integrated_snn_res.12"
ABCD_OFT <- subset(ABCD.subset,idents=c("3","30","53","56","96"))

saveRDS(ABCD_OFT,"//spatial_analysis/seurat/merge_result/cell-cell/edit/final/ABCD.OFT1.rds")

ABCD_OFT <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/edit/final/ABCD.OFT1.rds")
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### deconvlution
### scRNA-seq pre_process

CNCC_EMBO <- load("//spatial_analysis/CNCC_dataset/all.cncc.combined.EMBO.mapped.Rdata")

cncc.combined.sub.list <- SplitObject(cncc.combined.sub, split.by = "stage")

seuset_2 <- cncc.combined.sub.list$E13.5

Idents(seuset_2) <- "orig.ident"
seuset_2_ct <- subset(seuset_2,idents= c("ct2"))
seuset_2_mt <- subset(seuset_2,idents= c("vcl"))


Idents(seuset_2_ct) <- "EMBO_anno"
Idents(seuset_2_mt) <- "EMBO_anno"

seuset_2_ct <- SCTransform(seuset_2_ct, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

seuset_2_mt <- SCTransform(seuset_2_mt, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30) 

saveRDS(seuset_2_ct,"//spatial_analysis/seurat/merge_result/cell-cell/test/seuset_2_ct.rds")
saveRDS(seuset_2_mt,"//spatial_analysis/seurat/merge_result/cell-cell/test/seuset_2_mt.rds")

seuset_2_ct <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/test/seuset_2_ct.rds")
seuset_2_mt <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/test/seuset_2_mt.rds")

### subset yfp+ cell
cells_pos <- ABCD_OFT@meta.data[ABCD_OFT@assays$SCT@data["YFP",]>0.36,]
cells_neg <- ABCD_OFT@meta.data[ABCD_OFT@assays$SCT@data["YFP",]<=0.36,]

ABCD_OFT@meta.data$YFP_express<-"CC"

ABCD_OFT@meta.data[rownames(cells_pos),]$YFP_express<- "YFP+"

ABCD_OFT@meta.data[rownames(cells_neg),]$YFP_express<- "YFP-"

Idents(ABCD_OFT) <- "YFP_express"
ABCD_OFT.YFP <- subset(ABCD_OFT,idents=c("YFP+"))

Idents(ABCD_OFT.YFP) <- "integrated_snn_res.12"
SpatialSample <- SpatialDimPlot(ABCD_OFT.YFP, crop = FALSE, label = FALSE, pt.size.factor = 1.2, label.size = 4, ncol = 2, cols =c(
"3" = "#FF6C67","30" = "#E48C00","53" = "#A3A500" ,"56" = "#00B900" ,"96" = "#00C377"))
saveplot(plot_name,SpatialSample,"SpatialDimPlot_OFT_YFP+",store.path,24,24)



Idents(ABCD_OFT.YFP) <- "orig.ident"
A1_YFP <- subset(ABCD_OFT.YFP,idents=c("A1"))
B1_YFP <- subset(ABCD_OFT.YFP,idents=c("B1"))
C1_YFP <- subset(ABCD_OFT.YFP,idents=c("C1"))
D1_YFP <- subset(ABCD_OFT.YFP,idents=c("D1"))


anchors <- FindTransferAnchors(reference = seuset_2_ct, query = A1_YFP, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = seuset_2_ct$EMBO_anno, prediction.assay = TRUE,
    weight.reduction = A1[["pca"]], dims = 1:30,k.weight = 30)
A1_YFP[["predictions"]] <- predictions.assay

anchors <- FindTransferAnchors(reference = seuset_2_ct, query = B1_YFP, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = seuset_2_ct$EMBO_anno, prediction.assay = TRUE,
    weight.reduction = B1[["pca"]], dims = 1:30,k.weight = 30)
B1_YFP[["predictions"]] <- predictions.assay

anchors <- FindTransferAnchors(reference = seuset_2_mt, query = C1_YFP, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = seuset_2_mt$EMBO_anno, prediction.assay = TRUE,
    weight.reduction = C1[["pca"]], dims = 1:30,k.weight = 30)
C1_YFP[["predictions"]] <- predictions.assay

anchors <- FindTransferAnchors(reference = seuset_2_mt, query = D1_YFP, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = seuset_2_mt$EMBO_anno, prediction.assay = TRUE,
    weight.reduction = D1[["pca"]], dims = 1:30,k.weight = 30)
D1_YFP[["predictions"]] <- predictions.assay


DefaultAssay(A1_YFP) <- "predictions"
DefaultAssay(B1_YFP) <- "predictions"
DefaultAssay(C1_YFP) <- "predictions"
DefaultAssay(D1_YFP) <- "predictions"

A1_YFP@images$B1 <-NULL
A1_YFP@images$C1 <-NULL
A1_YFP@images$D1 <-NULL

B1_YFP@images$A1 <-NULL
B1_YFP@images$C1 <-NULL
B1_YFP@images$D1 <-NULL

C1_YFP@images$A1 <-NULL
C1_YFP@images$B1 <-NULL
C1_YFP@images$D1 <-NULL

D1_YFP@images$B1 <-NULL
D1_YFP@images$C1 <-NULL
D1_YFP@images$A1 <-NULL

### Figure S7A

Sample <- SpatialFeaturePlot(A1_YFP,pt.size.factor = 1.2,features=c("VSMC-1","VSMC-2","VSMC-3","VSMC-4","VSMC-5","VSMC-6","VSMC-7","VSMC-8") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_A1_VSMC",store.path,10,40)
Sample <- SpatialFeaturePlot(A1_YFP,pt.size.factor = 1.2,features=c("MES-1","MES-2","MES-3","MES-4") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_A1_MES",store.path,10,40)
Sample <- SpatialFeaturePlot(A1_YFP,pt.size.factor = 1.2,features=c("SWN","Mural","MLA","Neuron") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_A1_Other",store.path,10,40)

Sample <- SpatialFeaturePlot(B1_YFP,pt.size.factor = 1.2,features=c("VSMC-1","VSMC-2","VSMC-3","VSMC-4","VSMC-5","VSMC-6","VSMC-7","VSMC-8") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_B1_VSMC",store.path,10,40)
Sample <- SpatialFeaturePlot(B1_YFP,pt.size.factor = 1.2,features=c("MES-1","MES-2","MES-3","MES-4") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_B1_MES",store.path,10,40)
Sample <- SpatialFeaturePlot(B1_YFP,pt.size.factor = 1.2,features=c("SWN","Mural","MLA","Neuron") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_B1_Other",store.path,10,40)

Sample <- SpatialFeaturePlot(C1_YFP,pt.size.factor = 1.2,features=c("VSMC-1","VSMC-2","VSMC-3","VSMC-4","VSMC-5","VSMC-6","VSMC-7","VSMC-8") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_C1_VSMC",store.path,10,40)
Sample <- SpatialFeaturePlot(C1_YFP,pt.size.factor = 1.2,features=c("MES-1","MES-2","MES-3","MES-4") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_C1_MES",store.path,10,40)
Sample <- SpatialFeaturePlot(C1_YFP,pt.size.factor = 1.2,features=c("SWN","Mural","MLA","Neuron") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_C1_Other",store.path,10,40)


Sample <- SpatialFeaturePlot(D1_YFP,pt.size.factor = 1.2,features=c("VSMC-1","VSMC-2","VSMC-3","VSMC-4","VSMC-5","VSMC-6","VSMC-7","VSMC-8") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_D1_VSMC",store.path,10,40)
Sample <- SpatialFeaturePlot(D1_YFP,pt.size.factor = 1.2,features=c("MES-1","MES-2","MES-3","MES-4") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_D1_MES",store.path,10,40)
Sample <- SpatialFeaturePlot(D1_YFP,pt.size.factor = 1.2,features=c("SWN","Mural","MLA","Neuron") ,ncol = 1, crop = FALSE,alpha = c(0.1, 1))
saveplot(plot_name,Sample,"Predictions_D1_Other",store.path,10,40)




data<- data.frame(Cell_Type=rep(c("MESCNCC","VSMCCNCC","SHF","Other"),2),
  sample_group = rep(c("Control","VclcKO"),each = 4),
  count = c(0.077633592,0.5010258,0.421052632,0.000287976,0.203794773,0.422249852,0.35483871,0.019116665))

data$Cell_Type <- factor(data$Cell_Type, levels = c("VSMCCNCC", "MESCNCC", "SHF", "Other"))


gg <- ggplot(data, aes(x = sample_group, y = count, fill = Cell_Type))+
  geom_bar(stat = "identity", width = 0.7)+ # 柱状图绘制
  geom_flow(aes(alluvium = Cell_Type), alpha = 0.5) + # 添加柱状图后的条带
  scale_fill_manual(values = c("SHF"="#FE3CDB","VSMCCNCC"="#548235","MESCNCC"="#C5E0B4","Other"="#D9D9D9"),labels = c(expression(VSMC[CNCC]),expression(MES[CNCC]),"SHF","Other"))+
  labs(fill = "Cell Type")+
  theme_bw()+ # 将主题调整为白色背景和浅灰色网格线
  xlab("")+ # 去掉x轴的标题
  ylab("Relative Cell Composition (%)")+ # 设置y轴的标签
  scale_x_discrete(labels = c("Control",italic("Vcl")~"cKO"))+
  theme(panel.grid.major.x = element_blank(),
      axis.title.y = element_text(size=16),
        axis.text.x = element_text(size = rel(1.2)),
        axis.text.y = element_text(size=rel(0.8)),
        legend.text = element_text(size = rel(1),hjust =0))+ # 更改x轴、y轴的字体大小、刻度线等
  scale_y_continuous(labels = scales::percent_format(scale = 100),expand = c(0,0.05))+
  geom_text(aes(label = paste0(round(count*100), "%")), position = position_stack(vjust = 0.5), size = 3)

saveplot(plot_name,gg,"percentage",store.path,4,4)


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### Python Script




#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### Cell-Cell interaction using Cellchat R package
### CNCC-SHF level interaction analysis


### divivde the integrated clustes into YFP positive and YFP negative
cells_pos <- ABCD_OFT@meta.data[ABCD_OFT@assays$SCT@data["YFP",]>0.36,]
cells_neg <- ABCD_OFT@meta.data[ABCD_OFT@assays$SCT@data["YFP",]<=0.36,]


### rename the cluster in (clusterID) + (YFP+/-)
ABCD_OFT@meta.data$cluster_YFP<-"Null"
ABCD_OFT@meta.data[rownames(cells_pos),]$cluster_YFP<- paste0(ABCD_OFT@meta.data[rownames(cells_pos),]$integrated_snn_res.12,"_YFP")
ABCD_OFT@meta.data[rownames(cells_neg),]$cluster_YFP<- paste0(ABCD_OFT@meta.data[rownames(cells_neg),]$integrated_snn_res.12)


### rename the cluster in VSMC&MES&SHF level
ABCD_OFT@meta.data$anno_VSMC_MES_SHF <- "Null"
ABCD_OFT@meta.data[ABCD_OFT@meta.data$cluster_YFP=="3_YFP",]$anno_VSMC_MES_SHF <- "VSMC_CNCC"
ABCD_OFT@meta.data[ABCD_OFT@meta.data$cluster_YFP=="3",]$anno_VSMC_MES_SHF <- "SHF"
ABCD_OFT@meta.data[ABCD_OFT@meta.data$cluster_YFP=="30_YFP",]$anno_VSMC_MES_SHF <- "VSMC_CNCC"
ABCD_OFT@meta.data[ABCD_OFT@meta.data$cluster_YFP=="30",]$anno_VSMC_MES_SHF <- "SHF"
ABCD_OFT@meta.data[ABCD_OFT@meta.data$cluster_YFP=="53_YFP",]$anno_VSMC_MES_SHF <- "VSMC_CNCC"
ABCD_OFT@meta.data[ABCD_OFT@meta.data$cluster_YFP=="53",]$anno_VSMC_MES_SHF <- "SHF"
ABCD_OFT@meta.data[ABCD_OFT@meta.data$cluster_YFP=="56_YFP",]$anno_VSMC_MES_SHF <- "MES_CNCC"
ABCD_OFT@meta.data[ABCD_OFT@meta.data$cluster_YFP=="56",]$anno_VSMC_MES_SHF <- "SHF"
ABCD_OFT@meta.data[ABCD_OFT@meta.data$cluster_YFP=="96_YFP",]$anno_VSMC_MES_SHF <- "MES_CNCC"
ABCD_OFT@meta.data[ABCD_OFT@meta.data$cluster_YFP=="96",]$anno_VSMC_MES_SHF <- "SHF"


### rename the cluster in CNCC&SHF level
ABCD_OFT@meta.data$anno_CNCC_SHF <- ABCD_OFT@meta.data$anno_VSMC_MES_SHF
ABCD_OFT@meta.data[ABCD_OFT@meta.data$anno_VSMC_MES_SHF=="VSMC_CNCC",]$anno_CNCC_SHF <- "CNCC"
ABCD_OFT@meta.data[ABCD_OFT@meta.data$anno_VSMC_MES_SHF=="MES_CNCC",]$anno_CNCC_SHF <- "CNCC"


saveRDS(ABCD_OFT,"//spatial_analysis/seurat/merge_result/cell-cell/edit/final/ABCD.OFT3.rds")
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### Cell-Cell interaction using Cellchat R package
### CNCC-SHF level interaction analysis
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(CellChat)
library(ComplexHeatmap)

saveplot <- function(plot_name, plot, suffix, store.path, width, height) {
    filename <-paste(plot_name,"_",suffix,".pdf",sep = "")
  ggsave(filename, plot, path = store.path, width = width, height = height,dpi=600)
}

store.path <- c("//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_cnccshf/")

plot_name <- c("cellchat")

ABCD.subset <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/edit/ABCD.subset_cut.rds")


Idents(ABCD_OFT) <- "anno_CNCC_SHF"

### Split 4 sample and merge control group and mutant group

ABCD_OFT.list <- SplitObject(ABCD_OFT, split.by = "orig.ident")

DefaultAssay(ABCD_OFT.list$A1) <- "SCT"
DefaultAssay(ABCD_OFT.list$B1) <- "SCT"
DefaultAssay(ABCD_OFT.list$C1) <- "SCT"
DefaultAssay(ABCD_OFT.list$D1) <- "SCT"

subset_ct <- merge(ABCD_OFT.list$A1,ABCD_OFT.list$B1)
subset_mt <- merge(ABCD_OFT.list$C1,ABCD_OFT.list$D1)

#
data.input_ct = GetAssayData(subset_ct, slot = "data", assay = "SCT") # normalized data matrix
data.input_mt = GetAssayData(subset_mt, slot = "data", assay = "SCT") # normalized data matrix

Idents(subset_mt)<- "anno_CNCC_SHF"

meta_ct = data.frame(labels = Idents(subset_ct), row.names = names(Idents(subset_ct))) # manually create a dataframe consisting of the cell labels
meta_mt = data.frame(labels = Idents(subset_mt), row.names = names(Idents(subset_mt))) # manually create a dataframe consisting of the cell labels

unique(meta_ct$labels) # check the cell labels
unique(meta_mt$labels) # check the cell labels

### Get each sample coordinates
spatial.locs_A = GetTissueCoordinates(ABCD_OFT@images$A1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs_B = GetTissueCoordinates(ABCD_OFT@images$B1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs_C = GetTissueCoordinates(ABCD_OFT@images$C1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs_D = GetTissueCoordinates(ABCD_OFT@images$D1, scale = NULL, cols = c("imagerow", "imagecol")) 


### Add large constant to one of Control and Muntant group to make sure the coodinates can't overlap
### Addtionally, large constant avoid the two sample interaction in one group.
spatial.locs_B[1] = spatial.locs_B[1]+5000
spatial.locs_B[2] = spatial.locs_B[2]+5000

spatial.locs_D[1] = spatial.locs_D[1]+5000
spatial.locs_D[2] = spatial.locs_D[2]+5000

spatial.locs_ct <- rbind(spatial.locs_A,spatial.locs_B)
spatial.locs_mt <- rbind(spatial.locs_C,spatial.locs_D)

scale.factors_ct = jsonlite::fromJSON(txt = file.path("/Vcl_spatial_Elly_2022_Jul/cellranger/10XSpatial_Vcl/EVSE220525-VclCN-A1-1_report/outs/spatial", 'scalefactors_json.json'))
scale.factors_mt = jsonlite::fromJSON(txt = file.path("/Vcl_spatial_Elly_2022_Jul/cellranger/10XSpatial_Vcl/EVSE220525-VclMT-C1-1_report/outs/spatial", 'scalefactors_json.json'))

scale.factors_ct = list(spot.diameter = 65, spot = scale.factors_ct$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors_ct$fiducial_diameter_fullres, hires = scale.factors_ct$tissue_hires_scalef, lowres = scale.factors_ct$tissue_lowres_scalef)

scale.factors_mt = list(spot.diameter = 65, spot = scale.factors_mt$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors_mt$fiducial_diameter_fullres, hires = scale.factors_mt$tissue_hires_scalef, lowres = scale.factors_mt$tissue_lowres_scalef)


#### Create spatial Cellchat object
cellchat_ct <- createCellChat(object = data.input_ct, meta = meta_ct, group.by = "labels",datatype = "spatial", coordinates = spatial.locs_ct, scale.factors = scale.factors_ct)
cellchat_mt <- createCellChat(object = data.input_mt, meta = meta_mt, group.by = "labels",datatype = "spatial", coordinates = spatial.locs_mt, scale.factors = scale.factors_mt)

### Using Mouse database
CellChatDB.use <- CellChatDB.mouse # simply use the default CellChatDB
cellchat_ct@DB <- CellChatDB.use
cellchat_mt@DB <- CellChatDB.use


cellchat_ct <- subsetData(cellchat_ct) # This step is necessary even if using the whole database
cellchat_mt <- subsetData(cellchat_mt) # This step is necessary even if using the whole database

### We would like to use the same background gene to test the cell-cell interaction change between Control and Vcl cKO group.
### We transfered the background gene which infered by control group to Vcl cKO group
cellchat_ct <- identifyOverExpressedGenes(cellchat_ct)
cellchat_mt@var.features <- cellchat_ct@var.features

cellchat_ct <- identifyOverExpressedInteractions(cellchat_ct)
cellchat_mt <- identifyOverExpressedInteractions(cellchat_mt)

cellchat_ct <- computeCommunProb(cellchat_ct, type = "truncatedMean", trim = 0.1,distance.use = TRUE, interaction.length = 200, scale.distance = 0.01,k.min = 5)
cellchat_mt <- computeCommunProb(cellchat_mt, type = "truncatedMean", trim = 0.1,distance.use = TRUE, interaction.length = 200, scale.distance = 0.01,k.min = 5)

cellchat_ct <- computeCommunProbPathway(cellchat_ct)
cellchat_mt <- computeCommunProbPathway(cellchat_mt)

cellchat_ct <- aggregateNet(cellchat_ct)
cellchat_mt <- aggregateNet(cellchat_mt)

cellchat_ct <- netAnalysis_computeCentrality(cellchat_ct, slot.name = "netP")
cellchat_mt <- netAnalysis_computeCentrality(cellchat_mt, slot.name = "netP")


saveRDS(cellchat_ct, file = "//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_cnccshf/cellchat_ct_lowthreshold_distance_transfer.rds")
saveRDS(cellchat_mt, file = "//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_cnccshf/cellchat_mt_lowthreshold_distance_transfer.rds")

cellchat_ct <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_cnccshf/cellchat_ct_lowthreshold_distance_transfer.rds")
cellchat_mt <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_cnccshf/cellchat_mt_lowthreshold_distance_transfer.rds")

cellchat_ct@meta$labels = factor(cellchat_ct@meta$labels, levels = c("CNCC","SHF"))
cellchat_mt@meta$labels = factor(cellchat_mt@meta$labels, levels = c("CNCC","SHF"))

cellchat_ct <- updateClusterLabels(cellchat_ct,new.order = c("CNCC","SHF"))
cellchat_mt <- updateClusterLabels(cellchat_mt,new.order = c("CNCC","SHF"))

cellchat_ct <- updateCellChat(cellchat_ct)
cellchat_mt <- updateCellChat(cellchat_mt)

### Supplymental table 1
### Interaction in Control and Vcl cKO group
df.net_ct <- subsetCommunication(cellchat_ct,thresh = 0.05)
df.net_mt <- subsetCommunication(cellchat_mt,thresh = 0.05)
write.csv(df.net_ct, "//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_cnccshf/CCI_CNCC&SHF_Control.csv")
write.csv(df.net_mt, "//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_cnccshf/CCI_CNCC&SHF_VclcKO.csv")



### Figure 7c
### orginal Cellchat figure
object.list <- list(Control = cellchat_ct, Mutant = cellchat_mt)
cellchat <- mergeCellChat(object.list, add.names = c("Control","VclcKO")

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), color.use = c("#4DAA48","#FF0000"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight", color.use = c("#4DAA48","#FF0000"))
pl <- gg1 + gg2
saveplot(plot_name,pl,"compareInteractions",store.path,10,10)


### Replot in paper version
library(ggplot2)
color.use <- c("Control"="#47AA45","Vcl cKO"="#FF0000")
ylabel = "Number of inferred interactions"
color.alpha = 1
size.text = 14
data <- data.frame(category = c("Control", "Vcl cKO"), value = c(20,9))

# 绘制柱状图
gg <- ggplot(data, aes(x = category, y = value))
gg <- gg+  geom_bar(stat = "identity",fill = color.use,colour="black",size = 1.5,width = 0.6)
gg<- gg+ ylab(ylabel) + theme_classic() +  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
    theme(text = element_text(size = size.text), axis.text = element_text(colour="black"))

gg<- gg+ geom_text(aes(label = value), vjust = -0.5, size = 4)
# gg<- gg + scale_y_continuous(limits = c(0, max(data$value)*1.1))
gg<- gg+ scale_y_continuous(expand = c(0,0),limits = c(0, max(data$value)*1.1))
gg<- gg+theme(panel.background = element_rect(fill = "white"))

gg<- gg+theme(axis.title.x=element_blank())
gg<- gg+scale_x_discrete(labels = c("Control",italic("Vcl")~"cKO"))

saveplot(plot_name,gg,"cnccshf_interaction_number",store.path,2,4)





### Figure 7d
### Bubble plot with CNCC&SHF Ligand-Receptor
pdf("//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_cnccshf/clowthreshold_distance_transfer_bubble.pdf",width = 7,height = 5)
netVisual_bubble(cellchat,comparison = c(1,2),sources.use = c(1,2),targets.use = c(1,2) ,angle.x = 45,font.size = 10, color.text = c("#47AA45","#FF0000"))
dev.off()








#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### Cell-Cell interaction in VSMC_CNCC&MES_CNCC&SHF

saveplot <- function(plot_name, plot, suffix, store.path, width, height) {
    filename <-paste(plot_name,"_",suffix,".pdf",sep = "")
  ggsave(filename, plot, path = store.path, width = width, height = height,dpi=600)
}

store.path <- c("//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_vsmcmesshf/")

plot_name <- c("cellchat")

ABCD_OFT <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/edit/ABCD.subset_cut.rds")

cells_pos <- ABCD.subset_cut@meta.data[ABCD.subset_cut@assays$SCT@data["YFP",]>0.36,]
cells_neg <- ABCD.subset_cut@meta.data[ABCD.subset_cut@assays$SCT@data["YFP",]<=0.36,]


ABCD.subset_cut@meta.data$anno_spatial<-"CC"

ABCD.subset_cut@meta.data[rownames(cells_pos),]$anno_spatial<- paste0(ABCD.subset_cut@meta.data[rownames(cells_pos),]$integrated_snn_res.12,"_YFP")

ABCD.subset_cut@meta.data[rownames(cells_neg),]$anno_spatial<- paste0(ABCD.subset_cut@meta.data[rownames(cells_neg),]$integrated_snn_res.12)


ABCD.subset_cut@meta.data$anno_yfp <- "CC"

ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_spatial=="3_YFP",]$anno_yfp <- "VSMC_CNCC"
ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_spatial=="3",]$anno_yfp <- "VSMC_SHF"

ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_spatial=="30_YFP",]$anno_yfp <- "VSMC_CNCC"
ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_spatial=="30",]$anno_yfp <- "VSMC_SHF"

ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_spatial=="53_YFP",]$anno_yfp <- "VSMC_CNCC"
ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_spatial=="53",]$anno_yfp <- "VSMC_SHF"

ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_spatial=="56_YFP",]$anno_yfp <- "MES_CNCC"
ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_spatial=="56",]$anno_yfp <- "MES_SHF"

ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_spatial=="96_YFP",]$anno_yfp <- "MES_CNCC"
ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_spatial=="96",]$anno_yfp <- "MES_SHF"


ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_yfp=="VSMC_SHF",]$anno_yfp <- "SHF"
ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_yfp=="MES_SHF",]$anno_yfp <- "SHF"

ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_yfp=="VSMC_CNCC",]$anno_yfp <- "CNCC"
ABCD.subset_cut@meta.data[ABCD.subset_cut@meta.data$anno_yfp=="MES_CNCC",]$anno_yfp <- "CNCC"

Idents(ABCD.subset_cut) <- "anno_yfp"

ABCD.subset_cut.list <- SplitObject(ABCD.subset_cut, split.by = "orig.ident")

DefaultAssay(ABCD.subset_cut.list$A1) <- "SCT"
DefaultAssay(ABCD.subset_cut.list$B1) <- "SCT"
DefaultAssay(ABCD.subset_cut.list$C1) <- "SCT"
DefaultAssay(ABCD.subset_cut.list$D1) <- "SCT"

subset_ct <- merge(ABCD.subset_cut.list$A1,ABCD.subset_cut.list$B1)
subset_mt <- merge(ABCD.subset_cut.list$C1,ABCD.subset_cut.list$D1)

#
data.input_ct = GetAssayData(subset_ct, slot = "data", assay = "SCT") # normalized data matrix
data.input_mt = GetAssayData(subset_mt, slot = "data", assay = "SCT") # normalized data matrix

Idents(subset_ct)<- "anno_yfp"
Idents(subset_mt)<- "anno_yfp"

meta_ct = data.frame(labels = Idents(subset_ct), row.names = names(Idents(subset_ct))) # manually create a dataframe consisting of the cell labels
meta_mt = data.frame(labels = Idents(subset_mt), row.names = names(Idents(subset_mt))) # manually create a dataframe consisting of the cell labels

unique(meta_ct$labels) # check the cell labels
unique(meta_mt$labels) # check the cell labels


spatial.locs_A = GetTissueCoordinates(ABCD.subset_cut@images$A1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs_B = GetTissueCoordinates(ABCD.subset_cut@images$B1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs_C = GetTissueCoordinates(ABCD.subset_cut@images$C1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs_D = GetTissueCoordinates(ABCD.subset_cut@images$D1, scale = NULL, cols = c("imagerow", "imagecol")) 

spatial.locs_B[1] = spatial.locs_B[1]+5000
spatial.locs_B[2] = spatial.locs_B[2]+5000

spatial.locs_D[1] = spatial.locs_D[1]+5000
spatial.locs_D[2] = spatial.locs_D[2]+5000

spatial.locs_ct <- rbind(spatial.locs_A,spatial.locs_B)
spatial.locs_mt <- rbind(spatial.locs_C,spatial.locs_D)

scale.factors_ct = jsonlite::fromJSON(txt = file.path("/Vcl_spatial_Elly_2022_Jul/cellranger/10XSpatial_Vcl/EVSE220525-VclCN-A1-1_report/outs/spatial", 'scalefactors_json.json'))
scale.factors_mt = jsonlite::fromJSON(txt = file.path("/Vcl_spatial_Elly_2022_Jul/cellranger/10XSpatial_Vcl/EVSE220525-VclMT-C1-1_report/outs/spatial", 'scalefactors_json.json'))

scale.factors_ct = list(spot.diameter = 65, spot = scale.factors_ct$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors_ct$fiducial_diameter_fullres, hires = scale.factors_ct$tissue_hires_scalef, lowres = scale.factors_ct$tissue_lowres_scalef)

scale.factors_mt = list(spot.diameter = 65, spot = scale.factors_mt$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors_mt$fiducial_diameter_fullres, hires = scale.factors_mt$tissue_hires_scalef, lowres = scale.factors_mt$tissue_lowres_scalef)


#####
cellchat_ct <- createCellChat(object = data.input_ct, meta = meta_ct, group.by = "labels",datatype = "spatial", coordinates = spatial.locs_ct, scale.factors = scale.factors_ct)
cellchat_mt <- createCellChat(object = data.input_mt, meta = meta_mt, group.by = "labels",datatype = "spatial", coordinates = spatial.locs_mt, scale.factors = scale.factors_mt)


CellChatDB.use <- CellChatDB.mouse # simply use the default CellChatDB

cellchat_ct@DB <- CellChatDB.use
cellchat_mt@DB <- CellChatDB.use


cellchat_ct <- subsetData(cellchat_ct) # This step is necessary even if using the whole database
cellchat_mt <- subsetData(cellchat_mt) # This step is necessary even if using the whole database

cellchat_ct <- identifyOverExpressedGenes(cellchat_ct)
cellchat_mt <- identifyOverExpressedGenes(cellchat_mt)

cellchat_ct <- identifyOverExpressedInteractions(cellchat_ct)
cellchat_mt <- identifyOverExpressedInteractions(cellchat_mt)

cellchat_ct <- computeCommunProb(cellchat_ct, type = "truncatedMean", trim = 0.1,distance.use = TRUE, interaction.length = 200, scale.distance = 0.01,k.min = 5)
cellchat_mt <- computeCommunProb(cellchat_mt, type = "truncatedMean", trim = 0.1,distance.use = TRUE, interaction.length = 200, scale.distance = 0.01,k.min = 5)

cellchat_ct <- filterCommunication(cellchat_ct, min.cells = 1)
cellchat_mt <- filterCommunication(cellchat_mt, min.cells = 1)

cellchat_ct <- computeCommunProbPathway(cellchat_ct)
cellchat_mt <- computeCommunProbPathway(cellchat_mt)

cellchat_ct <- aggregateNet(cellchat_ct)
cellchat_mt <- aggregateNet(cellchat_mt)

cellchat_ct <- netAnalysis_computeCentrality(cellchat_ct, slot.name = "netP")
cellchat_mt <- netAnalysis_computeCentrality(cellchat_mt, slot.name = "netP")


saveRDS(cellchat_ct, file = "//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_vsmcmesshf/cellchat_ct_spatial.rds")
saveRDS(cellchat_mt, file = "//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_vsmcmesshf/cellchat_mt_spatial.rds")

cellchat_ct <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_vsmcmesshf/cellchat_ct_spatial.rds")
cellchat_mt <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_vsmcmesshf/cellchat_mt_spatial.rds")

cellchat_ct@meta$labels = factor(cellchat_ct@meta$labels, levels = c("VSMC_CNCC","MES_CNCC","SHF"))
cellchat_mt@meta$labels = factor(cellchat_mt@meta$labels, levels = c("VSMC_CNCC","MES_CNCC","SHF"))

cellchat_ct <- updateClusterLabels(cellchat_ct,new.order = c("VSMC_CNCC","MES_CNCC","SHF"))
cellchat_mt <- updateClusterLabels(cellchat_mt,new.order = c("VSMC_CNCC","MES_CNCC","SHF"))

cellchat_ct <- updateCellChat(cellchat_ct)
cellchat_mt <- updateCellChat(cellchat_mt)

df.net_ct <- subsetCommunication(cellchat_ct,thresh = 0.05)
df.net_mt <- subsetCommunication(cellchat_mt,thresh = 0.05)
write.csv(df.net_ct, "//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_vsmcmesshf/net_lr_ct.csv")
write.csv(df.net_mt, "//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_vsmcmesshf/net_lr_mt.csv")



object.list <- list(Control = cellchat_ct, Mutant = cellchat_mt)
cellchat <- mergeCellChat(object.list, add.names = c("Control","VclcKO")


### Figure 7J

pathways.show <- c("FN1") 
weight.max <- getMaxWeight(object.list[1], slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
pdf("//spatial_analysis/seurat/merge_result/cell-cell/edit/cellchat_vsmcmesshf/Chord_FN1.pdf")

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord",color.use = c("#548235","#C5E0B4","#FF00F7"), signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()




#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### Cell-Cell interaction in VSMC_CNCC & SHF which surrounded it
library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(CellChat)
library(ComplexHeatmap)


saveRDS(ABCD_OFT,"//spatial_analysis/seurat/merge_result/cell-cell/edit/final/ABCD.OFT3.rds")


ABCD.subset_cut <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/edit/final/ABCD.OFT3.rds")

saveplot <- function(plot_name, plot, suffix, store.path, width, height) {
    filename <-paste(plot_name,"_",suffix,".pdf",sep = "")
  ggsave(filename, plot, path = store.path, width = width, height = height)
}

store.path <- c("//spatial_analysis/seurat/merge_result/cell-cell/edit/stuility/")
plot_name <- c("")


ABCD_OFT.list <- SplitObject(ABCD_OFT, split.by = "orig.ident")

cluster_ids <- unique(ABCD_OFT.list$A1@meta.data$anno_VSMC_MES_SHF)
i=1
j=3
cells_i <- which(ABCD_OFT.list$A1@meta.data$anno_VSMC_MES_SHF == cluster_ids[i])
cells_j <- which(ABCD_OFT.list$A1@meta.data$anno_VSMC_MES_SHF == cluster_ids[j])
data <- rbind((ABCD_OFT.list$A1@images$A1@coordinates[4:5])[cells_i,],(ABCD_OFT.list$A1@images$A1@coordinates[4:5])[cells_j,])
dist_euclidean <- dist(data,method="euclidean")
dist_intergroup <- as.matrix(dist_euclidean)[1:table(ABCD_OFT.list$A1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]],(table(ABCD_OFT.list$A1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]]+1):(table(ABCD_OFT.list$A1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]]+table(ABCD_OFT.list$A1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[j]])]
num_pairs <- sum(dist_intergroup < 232.5)
barcodes_row <- rownames(dist_intergroup)
barcodes_col <- colnames(dist_intergroup)
barcodes_row_200_A1 <- barcodes_row[apply(dist_intergroup, 1, function(x) any(x < 232.5,na.rm = TRUE))]
barcodes_col_200_A1 <- barcodes_col[apply(dist_intergroup, 2, function(x) any(x < 232.5,na.rm = TRUE))]

cluster_ids <- unique(ABCD_OFT.list$B1@meta.data$anno_VSMC_MES_SHF)
i=2
j=1
cells_i <- which(ABCD_OFT.list$B1@meta.data$anno_VSMC_MES_SHF == cluster_ids[i])
cells_j <- which(ABCD_OFT.list$B1@meta.data$anno_VSMC_MES_SHF == cluster_ids[j])
data <- rbind((ABCD_OFT.list$B1@images$B1@coordinates[4:5])[cells_i,],(ABCD_OFT.list$B1@images$B1@coordinates[4:5])[cells_j,])
dist_euclidean <- dist(data,method="euclidean")
dist_intergroup <- as.matrix(dist_euclidean)[1:table(ABCD_OFT.list$B1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]],(table(ABCD_OFT.list$B1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]]+1):(table(ABCD_OFT.list$B1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]]+table(ABCD_OFT.list$B1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[j]])]
num_pairs <- sum(dist_intergroup < 232.5)
barcodes_row <- rownames(dist_intergroup)
barcodes_col <- colnames(dist_intergroup)
barcodes_row_200_B1 <- barcodes_row[apply(dist_intergroup, 1, function(x) any(x < 232.5,na.rm = TRUE))]
barcodes_col_200_B1 <- barcodes_col[apply(dist_intergroup, 2, function(x) any(x < 232.5,na.rm = TRUE))]


cluster_ids <- unique(ABCD_OFT.list$C1@meta.data$anno_VSMC_MES_SHF)
i=2
j=1
cells_i <- which(ABCD_OFT.list$C1@meta.data$anno_VSMC_MES_SHF == cluster_ids[i])
cells_j <- which(ABCD_OFT.list$C1@meta.data$anno_VSMC_MES_SHF == cluster_ids[j])
data <- rbind((ABCD_OFT.list$C1@images$C1@coordinates[4:5])[cells_i,],(ABCD_OFT.list$C1@images$C1@coordinates[4:5])[cells_j,])
dist_euclidean <- dist(data,method="euclidean")
dist_intergroup <- as.matrix(dist_euclidean)[1:table(ABCD_OFT.list$C1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]],(table(ABCD_OFT.list$C1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]]+1):(table(ABCD_OFT.list$C1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]]+table(ABCD_OFT.list$C1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[j]])]
num_pairs <- sum(dist_intergroup < 232.5)
barcodes_row <- rownames(dist_intergroup)
barcodes_col <- colnames(dist_intergroup)
barcodes_row_200_C1 <- barcodes_row[apply(dist_intergroup, 1, function(x) any(x < 232.5,na.rm = TRUE))]
barcodes_col_200_C1 <- barcodes_col[apply(dist_intergroup, 2, function(x) any(x < 232.5,na.rm = TRUE))]


cluster_ids <- unique(ABCD_OFT.list$D1@meta.data$anno_VSMC_MES_SHF)
i=2
j=1
cells_i <- which(ABCD_OFT.list$D1@meta.data$anno_VSMC_MES_SHF == cluster_ids[i])
cells_j <- which(ABCD_OFT.list$D1@meta.data$anno_VSMC_MES_SHF == cluster_ids[j])
data <- rbind((ABCD_OFT.list$D1@images$D1@coordinates[4:5])[cells_i,],(ABCD_OFT.list$D1@images$D1@coordinates[4:5])[cells_j,])
dist_euclidean <- dist(data,method="euclidean")
dist_intergroup <- as.matrix(dist_euclidean)[1:table(ABCD_OFT.list$D1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]],(table(ABCD_OFT.list$D1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]]+1):(table(ABCD_OFT.list$D1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[i]]+table(ABCD_OFT.list$D1@meta.data$anno_VSMC_MES_SHF)[cluster_ids[j]])]
num_pairs <- sum(dist_intergroup < 232.5)
barcodes_row <- rownames(dist_intergroup)
barcodes_col <- colnames(dist_intergroup)
barcodes_row_200_D1 <- barcodes_row[apply(dist_intergroup, 1, function(x) any(x < 232.5,na.rm = TRUE))]
barcodes_col_200_D1 <- barcodes_col[apply(dist_intergroup, 2, function(x) any(x < 232.5,na.rm = TRUE))]


merged_barcodes <- Reduce(union, list(barcodes_row_200_A1,barcodes_col_200_A1,barcodes_row_200_B1,barcodes_col_200_B1,barcodes_row_200_C1,barcodes_col_200_C1,barcodes_row_200_D1,barcodes_col_200_D1))
merged_barcodes <- as.character(merged_barcodes)

ABCD_OFT_boundary_VSMC_SHF <- subset(ABCD_OFT, cells = merged_barcodes)
saveRDS(ABCD_OFT_boundary_VSMC_SHF,"//spatial_analysis/seurat/merge_result/cell-cell/edit/final/ABCD_OFT_boundary_VSMC_SHF.rds")



Idents(ABCD_OFT_boundary_VSMC_SHF) <- "anno_VSMC_MES_SHF"
ABCD_OFT_boundary_VSMC_SHF.list <- SplitObject(ABCD_OFT_boundary_VSMC_SHF, split.by = "orig.ident")

DefaultAssay(ABCD_OFT_boundary_VSMC_SHF.list$A1) <- "SCT"
DefaultAssay(ABCD_OFT_boundary_VSMC_SHF.list$B1) <- "SCT"
DefaultAssay(ABCD_OFT_boundary_VSMC_SHF.list$C1) <- "SCT"
DefaultAssay(ABCD_OFT_boundary_VSMC_SHF.list$D1) <- "SCT"


ABCD_OFT_boundary_VSMC_SHF@meta.data$compare<-ABCD_OFT_boundary_VSMC_SHF@meta.data$orig.ident

Idents(ABCD_OFT_boundary_VSMC_SHF)<-"orig.ident"

control<-WhichCells(ABCD_OFT_boundary_VSMC_SHF,idents=c("A1","B1"))
mutant<- WhichCells(ABCD_OFT_boundary_VSMC_SHF,idents=c("C1","D1"))

ABCD_OFT_boundary_VSMC_SHF@meta.data[control,]$compare<-"control"
ABCD_OFT_boundary_VSMC_SHF@meta.data[mutant,]$compare<-"mutant"

ABCD_OFT_boundary_VSMC_SHF@meta.data$compare_heart<-paste0(ABCD_OFT_boundary_VSMC_SHF@meta.data$compare,"_",ABCD_OFT_boundary_VSMC_SHF@meta.data$anno_VSMC_MES_SHF)

DefaultAssay(ABCD_OFT_boundary_VSMC_SHF) <- "SCT"
data.input = GetAssayData(ABCD_OFT_boundary_VSMC_SHF, slot = "data", assay = "SCT") # normalized data matrix
Idents(ABCD_OFT_boundary_VSMC_SHF)<- "compare_heart"
meta = data.frame(labels = Idents(ABCD_OFT_boundary_VSMC_SHF), row.names = names(Idents(ABCD_OFT_boundary_VSMC_SHF))) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels


spatial.locs_A = GetTissueCoordinates(ABCD_OFT_boundary_VSMC_SHF@images$A1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs_B = GetTissueCoordinates(ABCD_OFT_boundary_VSMC_SHF@images$B1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs_C = GetTissueCoordinates(ABCD_OFT_boundary_VSMC_SHF@images$C1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs_D = GetTissueCoordinates(ABCD_OFT_boundary_VSMC_SHF@images$D1, scale = NULL, cols = c("imagerow", "imagecol")) 

spatial.locs_B[1] = spatial.locs_B[1]+5000
spatial.locs_B[2] = spatial.locs_B[2]+5000

spatial.locs_C[1] = spatial.locs_C[1]+10000
spatial.locs_C[2] = spatial.locs_C[2]+10000

spatial.locs_D[1] = spatial.locs_D[1]+15000
spatial.locs_D[2] = spatial.locs_D[2]+15000


spatial.locs <- rbind(spatial.locs_A,spatial.locs_B,spatial.locs_C,spatial.locs_D)

scale.factors = jsonlite::fromJSON(txt = file.path("/Vcl_spatial_Elly_2022_Jul/cellranger/10XSpatial_Vcl/EVSE220525-VclCN-A1-1_report/outs/spatial", 'scalefactors_json.json'))

scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef)


#####
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)

CellChatDB.use <- CellChatDB.mouse # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)

write.csv(cellchat@var.features, "//spatial_analysis/seurat/merge_result/cell-cell/edit/play/deg_vsmc.csv")

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,distance.use = TRUE, interaction.length = 200, scale.distance = 0.01,k.min = 5)

cellchat <- filterCommunication(cellchat, min.cells = 1)
#
cellchat <- computeCommunProbPathway(cellchat)


cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


saveRDS(cellchat, file = "//spatial_analysis/seurat/merge_result/cell-cell/edit/play/cellchat_spatial_vsmc.rds")
cellchat <- readRDS("//spatial_analysis/seurat/merge_result/cell-cell/edit/play/cellchat_spatial_vsmc.rds")
cellchat@meta$labels = factor(cellchat@meta$labels, levels = c("control_VSMC_CNCC","control_SHF","mutant_VSMC_CNCC","mutant_SHF"))
cellchat <- updateClusterLabels(cellchat,new.order = c("control_VSMC_CNCC","control_SHF","mutant_VSMC_CNCC","mutant_SHF"))
cellchat <- updateCellChat(cellchat)
df.net <- subsetCommunication(cellchat,thresh = 1)
write.csv(df.net, "//spatial_analysis/seurat/merge_result/cell-cell/edit/play/net_lr_vsmc.csv")


### Figure 7I
##all
i = 1
# 组合来自不同数据集的所有已识别的信号通路 
ht1 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "all", height = 18)
pdf("//spatial_analysis/seurat/merge_result/cell-cell/edit/play/vsmc_all.pdf",width =11,height = 20)
draw(ht1)
dev.off()











