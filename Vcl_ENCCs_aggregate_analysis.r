library(Seurat)
# library(singleCellTools)
library(Matrix)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(monocle)
library(scater)
library(corrplot)
# library(easyGgplot2)
library(devtools)
library(clusterProfiler) # GO & KEGG
library(gtools) # combination & permutations

set.seed(49)
options(stringsAsFactors = F)
options(repr.plot.width=4, repr.plot.height=3) # repr.plot.res = 300
options(mc.cores = 3)

myColors <- brewer.pal(8,"Set2")
myColors

# lineage
myColors5 <- c('#FDB462','#B3DE69','#FCCDE5',"#ff309a",'#D9D9D9')
names(myColors5) <- c("GP","BP","NPearly","NPlate","ENMFB")

Ctrl_YFP_ENCC <- "../../control/ENCC/raw/Ctrl_YFP_ENCC_report/outs/filtered_gene_bc_matrices/mm10"
Vcl_YFP_ENCC <- "raw/Vcl_YFP_ENCC_report/outs/filtered_gene_bc_matrices/mm10"

ctrl_rawdata <- Read10X(Ctrl_YFP_ENCC)
colnames(ctrl_rawdata) <- paste("ctrl", colnames(ctrl_rawdata), sep="_")
print(dim(ctrl_rawdata))
ctrl_rawdata[1:3,1:3]

E135Matrix.anno.st <- data.frame(cellName=colnames(ctrl_rawdata), origin.cluster="", 
                              stage="Control E13.5", platform="10x Genomics")
head(E135Matrix.anno.st)

vcl_rawdata <- Read10X(Vcl_YFP_ENCC)
colnames(vcl_rawdata) <- paste("vcl", colnames(vcl_rawdata), sep="_")
print(dim(vcl_rawdata))
vcl_rawdata[1:3,1:3]

vclMatrix.anno.st <- data.frame(cellName=colnames(vcl_rawdata), origin.cluster="", 
                              stage="Vcl cKO", platform="10x Genomics")
head(vclMatrix.anno.st)

metadata <- rbind(E135Matrix.anno.st, vclMatrix.anno.st)
rownames(metadata) <- metadata$cellName

set.seed(49)
all_vcl_ENCC_rawdata <- cbind(as.matrix(ctrl_rawdata),
                          as.matrix(vcl_rawdata))
print(dim(all_vcl_ENCC_rawdata))

options(repr.plot.width=5.5, repr.plot.height=3)
hist(log10(colSums(as.matrix(all_vcl_ENCC_rawdata))), breaks = 100, xlab = "log10 UMIs per cell", 
     main="all_CNCC_rawdata")

options(repr.plot.width=5.5, repr.plot.height=3)
hist(colSums(as.matrix(all_vcl_ENCC_rawdata)>0), breaks = 100, xlab = "gene number per cell", 
     main="all_CNCC_rawdata")

all_vcl_ENCC_rawdata <- all_vcl_ENCC_rawdata[,rownames(metadata)]

dim(all_vcl_ENCC_rawdata)
dim(metadata)

save(all_vcl_ENCC_rawdata, metadata, file="all.raw.vcl.integration.Rdata")

table(metadata$stage)

vcl.integration <- CreateSeuratObject(counts = all_vcl_ENCC_rawdata, meta.data = metadata)

vcl.list <- SplitObject(vcl.integration, split.by = "stage")

vcl.list <- lapply(X = vcl.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

vcl.anchors <- FindIntegrationAnchors(object.list = vcl.list, dims = 1:20)

vcl.combined <- IntegrateData(anchorset = vcl.anchors, dims = 1:20)

DefaultAssay(vcl.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
vcl.combined <- ScaleData(vcl.combined, verbose = FALSE)
vcl.combined <- RunPCA(vcl.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
vcl.combined <- RunUMAP(vcl.combined, reduction = "pca", dims = 1:20)
vcl.combined <- FindNeighbors(vcl.combined, reduction = "pca", dims = 1:20)

vcl.combined <- RunTSNE(vcl.combined, reduction = "pca", dims = 1:20)

vcl.combined <- FindClusters(vcl.combined, resolution = 0.3)

table(vcl.combined$seurat_clusters)

getwd()

print(load("kif7.final.clusters.Rdata"))

colnames(vcl.combined)[1:5]

tmp.clusters <- as.character(kif7.final.clusters[colnames(vcl.combined),])

tmp.clusters[is.na(tmp.clusters)] <- "others"

table(tmp.clusters)

names(tmp.clusters) <- colnames(vcl.combined)

vcl.combined$origin.cluster <- tmp.clusters

options(repr.plot.width=10, repr.plot.height=5)
p1 <- DimPlot(vcl.combined, reduction = "umap", group.by = "stage")
p2 <- DimPlot(vcl.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

options(repr.plot.width=10, repr.plot.height=5)
p1 <- DimPlot(vcl.combined, reduction = "pca", group.by = "stage")
p2 <- DimPlot(vcl.combined, reduction = "pca", label = TRUE)
plot_grid(p1, p2)

options(repr.plot.width=10, repr.plot.height=5)
p1 <- DimPlot(vcl.combined, reduction = "tsne", group.by = "stage")
p2 <- DimPlot(vcl.combined, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)

# lineage
myColors6 <- c('#FDB462','#B3DE69','#FCCDE5',"#ff309a",'#D9D9D9','white')
names(myColors6) <- c("GP","BP","NPearly","NPlate","ENMFB","others")

options(repr.plot.width=12, repr.plot.height=5)
p1 <- DimPlot(vcl.combined, reduction = "tsne", group.by = "origin.cluster", cols = myColors6)
p2 <- DimPlot(vcl.combined, reduction = "umap", group.by = "origin.cluster", cols = myColors6)
plot_grid(p1, p2)

vcl.combined$lineage.sub <- plyr::mapvalues(vcl.combined$seurat_clusters, 
                                                        from=0:8,
                                                        to=c("NPearly","GP","BP","NPlate","BP","BP","GP","ENMFB","ENMFB"))

vcl.combined$lineage.sub <- factor(vcl.combined$lineage.sub, levels = c("BP", "GP", "NPearly","NPlate","ENMFB"))

options(repr.plot.width=6, repr.plot.height=5)
p <- DimPlot(vcl.combined, reduction = "tsne", group.by = "lineage.sub", cols = myColors5)
p

options(repr.plot.width=12, repr.plot.height=5)
p1 <- DimPlot(vcl.combined, reduction = "umap", group.by = "lineage.sub", cols = myColors5)
p2 <- DimPlot(vcl.combined, reduction = "pca", group.by = "lineage.sub", cols = myColors5)
plot_grid(p1, p2)

options(repr.plot.width=10, repr.plot.height=5)
p <- DimPlot(vcl.combined, reduction = "tsne", group.by = "lineage.sub", cols = myColors5, split.by = "stage")
p

final.genes <- c('Sox10',"Plp1", 'Fabp7','Sparc','Rgcc',
                 'Ret','Gal','Nefm',
                 'Tubb3','Phox2b','Meg3','Stmn2','Stmn3',
                 'Pcsk1n','Rtn1','Cartpt','Prph','Sncg',
                 'Acta2','Mylk','Tagln','Col1a1','Dlk1','Myl9','Tpm2','Tpm1')

options(repr.plot.width=8, repr.plot.height=7)
hp <- DoHeatmap(vcl.combined, features = final.genes, group.by = "lineage.sub", hjust=0.5, 
        # cells = rownames(subset(all_umap, stage %in% c("E13.5","Kif7 cKO","E16.5"))),
        group.bar = T, group.colors = myColors5, 
        disp.min = -1.5, disp.max = 1.5, 
        slot = "scale.data", label = T,
        size = 3, angle = 0, raster = F,
        draw.lines = T, lines.width = 100, group.bar.height = 0.03, combine = T) + 
        # NoLegend() +
        scale_fill_gradientn(colors = PurpleAndYellow(), na.value = "white") + # PurpleAndYellow(), c("blue", "white", "red")
        theme(axis.text.y = element_text(size = 15, face="italic"))
hp

all_umap <- DimPlot(vcl.combined, reduction = "umap", group.by = "lineage.sub", cols = myColors5, split.by = "stage")$data

head(all_umap)

print(load("keyRdata/vcl.encc.integrated.Rdata"))

markers.to.plot <- c("Ret","Phox2b","Elavl4","Tubb3","Stmn2","Cartpt","Prph",
                     "Sox10","Fabp7","Plp1",
                     "Acta2","Tagln","Col1a1","Dlk1")

source("https://github.com/leezx/iterbi/raw/main/notebooks/unpackaged-code/dotplot.R")

options(repr.plot.width=8, repr.plot.height=4)
DotPlot_order(vcl.combined, features = markers.to.plot, dot.scale = 8, group.by = "lineage.sub", assay = "RNA") +
    RotatedAxis()

# dir.create("Figures")

ggsave(filename = "Figures/vcl.encc.dotplot.pdf", width = 8, height = 4)

source("https://github.com/leezx/iterbi/raw/main/notebooks/unpackaged-code/violin.R")

markers.to.plot <- c("Ret","Phox2b","Elavl4","Tubb3","Stmn2","Cartpt","Prph","Sox10","Fabp7","Plp1")

main.color

# options(repr.plot.width=13, repr.plot.height=5)
# p1 <- StackedVlnPlot.rowCluster(vcl.combined, markers.to.plot, group.by = "lineage.sub", cols = rev(myColors5))
# p1

vcl.combined$sample <- plyr::mapvalues(vcl.combined$stage, from = c("Control E13.5"), to = c("Control"))

head(all.meta)

all.meta <- vcl.combined@meta.data

all.meta$group <- paste(all.meta$sample, all.meta$lineage.sub)

exprM <- as.data.frame(data.frame(cbind(as.matrix(t(vcl.combined@assays$RNA@counts[markers.to.plot,]))), 
                             group.sub=as.character(all.meta[,c("group")])))

head(exprM)

exprM_melt <- reshape::melt(exprM, id.vars = c("group.sub"))

head(exprM_melt)

table(exprM_melt$group.sub)

exprM_melt$group.sub <- factor(exprM_melt$group.sub, levels = c("Control BP","Vcl cKO BP",
                                                                "Control NPearly","Vcl cKO NPearly",
                                                                "Control NPlate","Vcl cKO NPlate",
                                                                "Control GP","Vcl cKO GP",
                                                                "Control ENMFB","Vcl cKO ENMFB"
                                                               ))

options(repr.plot.width=6.5, repr.plot.height=5/6*8)
plot <- ggplot(exprM_melt, aes(x=variable, y=log2(1+value))) + 
    #geom_hline(yintercept = c(2,4,6), linetype="dashed", color = "grey80", size=0.3) +
    facet_wrap( ~ group.sub, ncol=1, labeller = label_context, scales = "free_y", strip.position = "left") +
    labs(x = "", y = "log2(UMIs+1)\n") +
    theme_bw() +
    # add background color to show the significance
    #geom_rect(data=subset(exprM_melt, group.sub %in% c("Control BP","Kif7 cKO BP")), xmin=0.4, xmax=1.5+1, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    #geom_rect(data=subset(exprM_melt, group.sub %in% c("Control BP","Kif7 cKO BP")), xmin=0.5+7, xmax=1.5+7, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    #geom_rect(data=subset(exprM_melt, group.sub %in% c("Control NPearly","Kif7 cKO NPearly")), xmin=0.4, xmax=1.5+3, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    #geom_rect(data=subset(exprM_melt, group.sub %in% c("Control NPlate","Kif7 cKO NPlate")), xmin=0.4, xmax=1.5+6, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    #geom_rect(data=subset(exprM_melt, group.sub %in% c("Control GP","Kif7 cKO GP")), xmin=0.5+1, xmax=1.5+1, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    #geom_rect(data=subset(exprM_melt, group.sub %in% c("Control GP","Kif7 cKO GP")), xmin=0.5+7, xmax=2.5+7, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    #
    geom_violin(trim = F, na.rm = T, aes(fill = group.sub),colour = "black", scale="width", size=0.3, bw=0.3) + 
    theme(strip.background = element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.spacing=unit(.4, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
    theme(axis.text.x  = element_text(face="italic", angle=90, size = 18, color = "black", vjust=0.5),
        axis.text.y  = element_text(face="plain", size = 8, color = "black"), # y axis scale
        axis.title = element_text(size = 18)) +
    theme(strip.background = element_rect(fill = NA, color = NA), legend.position = "none") +
    theme(strip.placement = "outside") +
    theme(strip.text.y.left = element_text(angle=0, size = 15, face="bold"),
        # strip.text.y = element_text(angle=90, margin = margin(1,1,1,1, "mm"), size = 15, face="bold"), 
          axis.ticks.x = element_blank()) +
    scale_y_continuous(position="right", limits=c(0, 8), breaks = seq(0, 8, length.out = 3)) #+ #
    #scale_y_continuous(position="right") #+
plot

ggsave(filename = "Figures/vcl.encc.violin.pdf", width = 6.5, height = 5/6*8)



library(dplyr)
countTable <- all_umap %>% dplyr::group_by(stage) %>% dplyr::count(lineage.sub) %>% dplyr::mutate(prop = n/sum(n))

options(repr.plot.width=4, repr.plot.height=7)
library(scales)
ggplot(data=countTable, aes(x=stage, y=n, fill=lineage.sub)) +
  geom_bar(stat="identity", position="fill", alpha=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "",y = "Percentage of Population", title = " ") + 
  theme(axis.text.x  = element_text(face="plain", angle=70, size = 15, color = "black", vjust=0.5),
        axis.text.y  = element_text(face="plain", size = 15, color = "black"),
        axis.title =element_text(size = 15)) +
  scale_y_continuous(labels = percent_format()) +
  # coord_flip() +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values=rev(myColors5[c("BP","GP","NPearly","NPlate","ENMFB")]), guide = guide_legend(reverse=F))

countTable

print(load("keyRdata/vcl.encc.integrated.Rdata"))

table(vcl.combined$lineage.sub)

vcl.combined@active.ident <- vcl.combined$lineage.sub

vcl.encc.markers <- FindAllMarkers(vcl.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

head(vcl.encc.markers)

table(vcl.encc.markers$cluster)

vcl.encc.markers$`pct.1` <- NULL
vcl.encc.markers$`pct.2` <- NULL

library(xlsx)
tmp.file <- "unsupervised.markers.5.clusters.encc.E13.5.xlsx"
write.xlsx(vcl.encc.markers, file=tmp.file, sheetName="markers", row.names=T)



save(vcl.combined, all_umap, vcl.encc.markers, file = "keyRdata/vcl.encc.integrated.Rdata")





# prepare data
all_tsne <- all_umap
all_tsne$group <- all_tsne$stage
all_tsne$cluster <- all_tsne$lineage.sub

seuset <- vcl.combined

# seuset@assays$RNA@counts

table(all_tsne$group)
table(all_tsne$cluster)

DEGs_list_full <- list()
two_group <- c("Control E13.5","Vcl cKO")

for (i in unique(all_tsne$cluster)) {
    tmp_all_tsne <- all_tsne[all_tsne$cluster==i & all_tsne$group %in% two_group,]
    # get raw count data
    tmp_rawdata <- seuset@assays$RNA@counts[,rownames(tmp_all_tsne)]
    # down sampling
    num <- min(min(table(tmp_all_tsne$group)[table(tmp_all_tsne$group)!=0]),500)
    set.seed(49)
    cells <- c(sample(rownames(tmp_all_tsne[tmp_all_tsne$group==two_group[1],]), num),
              sample(rownames(tmp_all_tsne[tmp_all_tsne$group==two_group[2],]), num))
    tmp_all_tsne <- tmp_all_tsne[cells,]
    tmp_rawdata <- tmp_rawdata[,cells]
    row_anno <- data.frame(gene_short_name=rownames(tmp_rawdata), row.names = rownames(tmp_rawdata))
    # break
    # from Monocle
    pd <- new("AnnotatedDataFrame", data = tmp_all_tsne)
    fd <- new("AnnotatedDataFrame", data = row_anno)

    Obj <- newCellDataSet(
        tmp_rawdata, 
        phenoData = pd, 
        featureData = fd,
        expressionFamily = negbinomial.size()
    )
    Obj <- estimateSizeFactors(Obj)
    Obj <- estimateDispersions(Obj)

    res <- differentialGeneTest(Obj, fullModelFormulaStr = "~group")

    pVals <- res[,3]
    names(pVals) <- rownames(res)
    pVals <- p.adjust(pVals, method = "fdr")

    DEGs <- data.frame(pVals, gene=names(pVals))
    # check
    if (table(tmp_all_tsne$group)[two_group[1]] != table(tmp_all_tsne$group)[two_group[2]]) {print("count error!");break;}
    # add 1 to avoid NA and finite value
    DEGs$log2FC <- log2((1+apply(tmp_rawdata[,tmp_all_tsne$group==two_group[2]],1,sum)) / (1+apply(tmp_rawdata[,tmp_all_tsne$group==two_group[1]],1,sum)))[DEGs$gene]
    DEGs$obj.pct <- (rowSums(tmp_rawdata[,tmp_all_tsne$group==two_group[2]]>0) / sum(tmp_all_tsne$group==two_group[2]))[DEGs$gene]
    DEGs$bcg.pct <- (rowSums((tmp_rawdata[,tmp_all_tsne$group==two_group[1]])>0) / sum(tmp_all_tsne$group==two_group[1]))[DEGs$gene]
    # DEGs <- subset(DEGs, pVals<0.05)
    #sig_DEGs <- rbind(subset(DEGs, pVals<0.05 & log2FC>0.5 & obj.pct>0.25), subset(DEGs, pVals<0.05 & log2FC< -0.5 & bcg.pct>0.25))

    DEGs_list_full[[i]] <- DEGs
    #sig_DEGs_list_full[[i]] <- sig_DEGs
    # sig_DEGs[[i]] <- rownames(sig_DEGs)
    # DEGsList <- list()
    #sig_DEGs_up[[i]] <- sig_DEGs[sig_DEGs$log2FC>0,]$gene
    #sig_DEGs_down[[i]] <- sig_DEGs[sig_DEGs$log2FC<0,]$gene
    # break
    }

save(DEGs_list_full, file = "result/vcl_encc_DEGs_list_full.Rdata")

print(load("keyRdata/vcl_encc_DEGs_list_full.Rdata"))

names(DEGs_list_full)

DEG.sig <- lapply(DEGs_list_full, function(x) {
    subset(x, pVals<0.05 & abs(log2FC)>0)[,c("gene","pVals","log2FC")]
})

out.file <- "cluster.specific.DEGs.xlsx"
tmp.list <- DEG.sig
sample.list <- unique(names(tmp.list))
for (i in sample.list) {
    print(i)
    tmp.df <- tmp.list[[i]]
    # tmp.df <- subset(tmp.df, pvalue<0.05)
    print(dim(tmp.df))
    if (i==sample.list[1]) {
        write.xlsx(tmp.df, file=out.file, sheetName=i, row.names=F)
    } else {
        write.xlsx(tmp.df, file=out.file, sheetName=i, append=TRUE, row.names=F)
    }
}



# redo sig
sig_DEGs_up <- list()
sig_DEGs_down <- list()

for (i in names(DEGs_list_full)) {
    sig_DEGs_up[[i]] <- subset(DEGs_list_full[[i]], pVals<0.05 & log2FC>0)$gene
    sig_DEGs_down[[i]] <- subset(DEGs_list_full[[i]], pVals<0.05 & log2FC< 0)$gene
    }

names(sig_DEGs_up) <- paste(names(sig_DEGs_up), " (", lapply(sig_DEGs_up, length), " genes up)", sep="")
names(sig_DEGs_down) <- paste(names(sig_DEGs_down), " (", lapply(sig_DEGs_down, length), " genes down)", sep="")

lapply(sig_DEGs_up, length)
lapply(sig_DEGs_down, length)

result.up <- ora.go.kegg.clusterProfiler(geneList = sig_DEGs_up, organism="mm")

options(repr.plot.width=10, repr.plot.height=6)
go_list <- plot.ora.GO.KEGG.barplot.batch(result.up$go_list, type="GO")

options(repr.plot.width=10, repr.plot.height=6)
kegg_list <- plot.ora.GO.KEGG.barplot.batch(result.up$kegg_list, type="KEGG")

result.down <- ora.go.kegg.clusterProfiler(geneList = sig_DEGs_down, organism="mm")

options(repr.plot.width=10, repr.plot.height=6)
go_list <- plot.ora.GO.KEGG.barplot.batch(result.down$go_list, type="GO")

options(repr.plot.width=10, repr.plot.height=6)
kegg_list <- plot.ora.GO.KEGG.barplot.batch(result.down$kegg_list, type="KEGG")

save(result.up, result.down, file="DEG.GO.KEGG.anno.Rdata")

print(load("keyRdata/DEG.GO.KEGG.anno.Rdata"))

tmp.go.NPearly <- result.down$go_list$`NPearly (188 genes down)`@result[c("GO:0007409","GO:0050770"),] # ,"GO:0007411"

tmp.go.NPearly$type <- "GO"

tmp.kegg.NPearly <- result.down$kegg_list$`NPearly (188 genes down)`@result[c("mmu04010","mmu04360"),]

# result.down$kegg_list$`NPearly (188 genes down)`@result

tmp.kegg.NPearly$type <- "KEGG"



tmp.go.NPlate <- result.down$go_list$`NPlate (248 genes down)`@result[c("GO:0031103","GO:1990138","GO:0048675"),]

tmp.go.NPlate$type <- "GO"

tmp.kegg.NPlate <- result.down$kegg_list$`NPlate (248 genes down)`@result[c("mmu04360","mmu04514","mmu04070"),]

# result.down$kegg_list$`NPlate (248 genes down)`@result

tmp.kegg.NPlate$type <- "KEGG"



tmp.NPearly <- rbind(tmp.go.NPearly, tmp.kegg.NPearly)

tmp.NPlate <- rbind(tmp.go.NPlate, tmp.kegg.NPlate)



tmp.NPearly

source("https://github.com/leezx/iterbi/raw/main/notebooks/unpackaged-code/pathway_enrichment.R")

options(repr.plot.width=5, repr.plot.height=3)
plot.GO.barplot(tmp.NPearly)

plot.GO.barplot2(tmp.NPearly)



options(repr.plot.width=5, repr.plot.height=3)
plot.GO.barplot.pair(tmp.NPearly, color = 1:2)

ggsave(filename = "Figures/NPearly.GO.KEGG.pdf", width = 5, height = 3)

options(repr.plot.width=5, repr.plot.height=3.5)
plot.GO.barplot.pair(tmp.NPlate, color = 3:4)

ggsave(filename = "Figures/NPlate.GO.KEGG.pdf", width = 5, height = 3.5)





result <- result.up
for (i in names(result$go_list)) {
    print(i)
    tmp.df <- result$go_list[[i]]@result[,c("Description","pvalue","geneID","Count")]
    write.csv(tmp.df, file=paste("result/", i, ".csv", sep=""))
    # break
}

result <- result.down
for (i in names(result$go_list)) {
    print(i)
    tmp.df <- result$go_list[[i]]@result[,c("Description","pvalue","geneID","Count")]
    write.csv(tmp.df, file=paste("result/", i, ".csv", sep=""))
    # break
}

options(warn = -1)

result <- result.up
for (i in names(result$kegg_list)) {
    print(i)
    tmp.df <- result$kegg_list[[i]]@result[,c("Description","pvalue","geneID","Count")]
    tmp.df$Genes <- unlist(lapply(tmp.df$geneID, function(x){ ENTREZIDtoSYMBOL(x,organism = "mm", returnVector = F) }[1][1]))
    write.csv(tmp.df, file=paste("result/", i, ".csv", sep=""))
    # break
}

result <- result.down
for (i in names(result$kegg_list)) {
    print(i)
    tmp.df <- result$kegg_list[[i]]@result[,c("Description","pvalue","geneID","Count")]
    tmp.df$Genes <- unlist(lapply(tmp.df$geneID, function(x){ ENTREZIDtoSYMBOL(x,organism = "mm", returnVector = F) }[1][1]))
    write.csv(tmp.df, file=paste("result/", i, ".csv", sep=""))
    # break
}

ENTREZIDtoSYMBOL <- function (ID, organism = "hs", returnVector = T) 
{
    ID <- unlist(lapply(strsplit(ID, "/"), as.integer))
    if (organism == "hs") {
        library(org.Hs.eg.db)
        gene.df <- bitr(ID, fromType = "ENTREZID", toType = c("SYMBOL", 
            "ENSEMBL"), OrgDb = org.Hs.eg.db)
    }
    else if (organism == "mm") {
        library(org.Mm.eg.db)
        gene.df <- bitr(ID, fromType = "ENTREZID", toType = c("SYMBOL", 
            "ENSEMBL"), OrgDb = org.Mm.eg.db)
    }
    else {
        stop("only support hs and mm now!")
    }
    gene.df <- gene.df[!duplicated(gene.df$ENTREZID), ]
    rownames(gene.df) <- gene.df$ENTREZID
    if (returnVector) {
        return(gene.df$SYMBOL)
    }
    else {
        genes <- paste(gene.df$SYMBOL, collapse = "/")
        genes
    }
}

ls()

packageVersion("CellChat")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

head(CellChatDB$interaction, n=10)
dim(CellChatDB$interaction)

# table(CellChatDB$interaction$receptor)

tmp.df <- CellChatDB$interaction

tmp.df$geneL <- unlist(lapply(tmp.df$interaction_name_2, function(x){
    y1 <- strsplit(x, split = " - ")[[1]][1]
    y1 <- trimws(y1)
}) )

tmp.df$geneR1 <- ""
tmp.df$geneR2 <- ""

tmp.df1 <- tmp.df[!grepl("\\+", tmp.df$interaction_name_2), ]
head(tmp.df1)

tmp.df2 <- tmp.df[grepl("\\+", tmp.df$interaction_name_2), ]
head(tmp.df2)

# receptor class 1
tmp.df1$geneR1 <- unlist(lapply(tmp.df1$interaction_name_2, function(x){
    y1 <- strsplit(x, split = " - ")[[1]][2]
    y1 <- trimws(y1)
}) )



# receptor class 2
tmp.df2$geneR1 <- unlist(lapply(tmp.df2$interaction_name_2, function(x){
    y1 <- strsplit(x, split = " - ")[[1]][2]
    y1 <- trimws(y1)
    y.len <- nchar(y1)
    y1 <- substr(y1, start = 2, stop = y.len-1)
    y2 <- strsplit(y1, split = "\\+")[[1]][1]
    # unlist(y2)
}) )

# receptor class 2
tmp.df2$geneR2 <- unlist(lapply(tmp.df2$interaction_name_2, function(x){
    y1 <- strsplit(x, split = " - ")[[1]][2]
    y1 <- trimws(y1)
    y.len <- nchar(y1)
    y1 <- substr(y1, start = 2, stop = y.len-1)
    y2 <- strsplit(y1, split = "\\+")[[1]][2]
    # unlist(y2)
}) )

head(tmp.df2)

tmp.melt.df <- rbind(tmp.df1, tmp.df2)

# tmp.melt.df

tmp.melt.df$L_DE <- tmp.melt.df$geneL %in% c(sig_DEGs_down$`NPearly (188 genes down)`, sig_DEGs_down$`NPlate (248 genes down)`)

tmp.melt.df$R1_DE <- tmp.melt.df$geneR1 %in% c(sig_DEGs_down$`NPearly (188 genes down)`, sig_DEGs_down$`NPlate (248 genes down)`)

tmp.melt.df$R2_DE <- tmp.melt.df$geneR2 %in% c(sig_DEGs_down$`NPearly (188 genes down)`, sig_DEGs_down$`NPlate (248 genes down)`)

tmp.melt.df$DE_count <- rowSums(tmp.melt.df[,c("L_DE","R1_DE","R2_DE")])

cellchat.DE <- tmp.melt.df[tmp.melt.df$DE_count > 0,]

cellchat.DE <- cellchat.DE[,c("pathway_name","annotation","interaction_name_2","L_DE","R1_DE","R2_DE","DE_count")]

write.csv(cellchat.DE, file = "Vcl.ENCC.cellchat.DE.csv")

table(cellchat.DE$pathway_name)

length(table(cellchat.DE$pathway_name))

cellchat.DE <- read.csv("Vcl.ENCC.cellchat.DE.csv", row.names = 1, stringsAsFactors = F)

head(cellchat.DE)

length(unique(cellchat.DE$pathway_name))

print(load("keyRdata/vcl.encc.integrated.Rdata"))

table(vcl.combined$lineage.sub)

vcl.combined.sub <- subset(vcl.combined, lineage.sub %in% c("BP","GP","NPearly","NPlate"))

# Normalized data (e.g., library-size normalization and then log-transformed with a pseudocount of 1
data.input <- normalizeData(vcl.combined.sub@assays$RNA@counts)

dim(data.input)

meta <- vcl.combined.sub@meta.data

head(meta)

meta$labels <- as.character(meta$lineage.sub)

table(meta$labels)

levels(meta$labels)

table(meta$stage)

cell.use <- rownames(subset(meta, stage=="Control E13.5"))
# cell.use <- rownames(subset(meta, stage=="Vcl cKO"))
head(cell.use)

# Prepare input data for CelChat analysis
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# speed up using specific pathways
CellChatDB.use <- subsetDB(CellChatDB, search = cellchat.DE$pathway_name, key = "pathway_name")

# speed up using specific pathways
CellChatDB.use <- subsetDB(CellChatDB, search = rownames(cellchat.DE), key = "interaction_name")

head(CellChatDB.use$interaction)

table(CellChatDB.use$interaction$pathway_name)

length(unique(CellChatDB.use$interaction$pathway_name))

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 3) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# mat <- cellchat@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }

# pathways.show <- c("TGFb") 
# # Hierarchy plot
# # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
# vertex.receiver = seq(1,4) # a numeric vector. 
# netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# # Circle plot
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# print(load("keyRdata/cellchat.ctrl.Rdata"))

pathways.show <- c("PTN") 

pathways.show <- c("CADM") 

pathways.show <- c("NCAM") 

# pathways.show <- c("MAG") 

# Chord diagram
par(mfrow=c(1,1))
p <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", color.use = myColors,
                    show.legend = F)
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# ggsave(filename = "Figures/TGFb3.circle.pdf", width = 5, height = 5, dpi = 800)

# this is a recordedplot OBJ
pander::redrawPlot(p)

pdf(paste("Figures/", pathways.show, "_Control.pdf", sep = ""), width=6, height=6)
pander::redrawPlot(p)
dev.off()

# # Heatmap
# par(mfrow=c(1,1))
# netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# #> Do heatmap based on a single object

# # Chord diagram
# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# names(group.cellType) <- levels(cellchat@idents)
# netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
# #> Plot the aggregated cell-cell communication network at the signaling pathway level
# #> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netAnalysis_contribution(cellchat, signaling = pathways.show)

cellchat@netP$pathways

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
pathways.show.all <- c("CADM","NCAM","PTN")

netVisual(cellchat, signaling = "CADM")

netAnalysis_contribution(cellchat, signaling = "NCAM")

netVisual_individual(cellchat, signaling = "CADM", layout = "circle") # circle

pathways.show <- "PTN" # "CADM","NCAM","PTN"

pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[1,]

pairLR

LR.show <- "PTN_SDC3"

g <- netVisual_individual(cellchat, signaling = pathways.show, layout = "hierarchy", pairLR.use = LR.show,
                     vertex.receiver = seq(1,2)) # circle

pdf(paste("Figures/", pathways.show, "_hierarchy.pdf", sep = ""), width=8, height=4)
pander::redrawPlot(g)
dev.off()

# # check the order of cell identity to set suitable vertex.receiver
# levels(cellchat@idents)
# vertex.receiver = seq(1,4)
# for (i in 1:length(pathways.show.all)) {
#   # Visualize communication network associated with both signaling pathway and individual L-R pairs
#   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
#   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#   ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
# }

cellchat <- computeCommunProbPathway(cellchat)

# Heatmap
# par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "CADM", color.heatmap = "Reds")

netVisual_bubble(cellchat, remove.isolate = FALSE, signaling = c("CADM","PTN","NCAM","MAG"))

ggsave(filename="Figures/dotplot.pdf", width = 5, height = 5, units = 'in', dpi = 300)



plotGeneExpression(cellchat, signaling = "PTN")

# # Compute the network centrality scores
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
# netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

save(cellchat, file = "keyRdata/cellchat.ctrl.Rdata")



print(load("keyRdata/cellchat.ctrl.Rdata"))

pathways.show <- c("PTN") 

pathways.show <- c("NCAM") 

# # Heatmap
# par(mfrow=c(1,1))
# netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# #> Do heatmap based on a single object

# # Chord diagram
# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# names(group.cellType) <- levels(cellchat@idents)
# netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
# #> Plot the aggregated cell-cell communication network at the signaling pathway level
# #> Note: The first link end is drawn out of sector 'Inflam. FIB'.

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)

pairLR

LR.show <- "PTN_NCL"

LR.show <- "NCAM1_L1CAM"

LR.show <- "NCAM1_NCAM1"

table(cellchat@meta$labels)

cellchat@idents <- factor(cellchat@idents, levels = c("BP","GP","NPearly","NPlate"))

g <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, color.use = myColors5[c(2,1,3,4)],
                     layout = "hierarchy", vertex.receiver = seq(1,2)) # circle

pdf(paste("Figures/", LR.show, "_hierarchy.pdf", sep = ""), width=8, height=4)
pander::redrawPlot(g)
dev.off()

# g <- netVisual_individual(cellchat, signaling = pathways.show, layout = "circle", pairLR.use = LR.show,
#                      vertex.receiver = seq(1,2)) # circle

g <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, color.use = myColors5[1:4],
                     layout = "chord", vertex.receiver = seq(1,2)) # circle

pdf(paste("Figures/", LR.show, "_chord.pdf", sep = ""), width=6, height=6)
pander::redrawPlot(g)
dev.off()

print(load("keyRdata/vcl.encc.integrated.Rdata"))

table(vcl.combined$stage)

vcl.combined$sample <- plyr::mapvalues(vcl.combined$stage, from = c("Control E13.5"), to = c("Control"))

tmp.genes <- c("Ptn","Ncl","Ncam1","Sdc3","Sdc1","Ptprz1","L1cam","Fgfr1","Cadm3")

vcl.combined <- ScaleData(vcl.combined, features = tmp.genes, assay = "RNA")

cellanno <- data.frame(group = vcl.combined$sample, 
                       row.names = colnames(vcl.combined), 
                       cluster = vcl.combined$lineage.sub)

head(cellanno)

cellannoExp <- cbind(cellanno, t(as.matrix(vcl.combined@assays$RNA@data[tmp.genes, rownames(cellanno)])))

# cellannoExp <- cbind(cellanno, t(as.matrix(vcl.combined@assays$RNA@scale.data[tmp.genes, rownames(cellanno)])))

head(cellannoExp)

cellannoExpMelt <- reshape2::melt(cellannoExp, id.vars = c("group", "cluster"))

colnames(cellannoExpMelt) <- c('group', 'cluster', 'gene', 'expression')

head(cellannoExpMelt)

cellannoExpMelt <- subset(cellannoExpMelt, cluster %in% c("BP","NPearly","NPlate"))



# all scale to 0-1
cellannoExpMelt2 <- data.frame()
cellannoExpMelt$scale <- 0
for (i in unique(cellannoExpMelt$gene)) {
    print(i)
    tmp.df <- subset(cellannoExpMelt, gene==i)
    tmp.df <- tmp.df %>% filter(expression < quantile(tmp.df$expression, 0.95))
    tmp.df$scale <- (tmp.df$expression - min(tmp.df$expression)) / (max(tmp.df$expression) - min(tmp.df$expression))
    tmp.df[is.na(tmp.df)] <- 0 
    cellannoExpMelt2 <- rbind(cellannoExpMelt2, tmp.df)
    # break
}

dat_text <- data.frame()
for (i in unique(cellannoExpMelt$cluster)) {
    for (j in unique(cellannoExpMelt$gene)) {
        direction <- "UP"
        # print(i)
        # print(j)
        log2fc <- log2(1+mean(subset(cellannoExpMelt, cluster==i & gene==j & group=="Vcl cKO")$expression)) - 
                      log2(1+mean(subset(cellannoExpMelt, cluster==i & gene==j & group=="Control")$expression))
        # log2fc <- (mean(subset(cellannoExpMelt, cluster==i & gene==j & group=="PHOX2B-7Ala")$expression)) - 
        #              (mean(subset(cellannoExpMelt, cluster==i & gene==j & group=="Control")$expression))
        # print(log2fc)
        if (log2fc<0) direction <- "DOWN"
        # if (abs(log2fc) < 0.02) {next}
        dat_text <- rbind(dat_text, c(i,j,log2fc,direction))
    }
}
colnames(dat_text) <- c("cluster", "gene", "log2FC","label")

dat_text$cluster <- factor(dat_text$cluster)
dat_text$gene <- factor(dat_text$gene)
dat_text$group <- "Control"

library(ggpubr)

head(cellannoExpMelt)

# install.packages("rstatix")

?wilcox_test

library(rstatix)
stat.test <- cellannoExpMelt %>% 
                group_by(cluster, gene) %>%
                rstatix::t_test(formula = expression ~ group) %>%
                adjust_pvalue(method = 'bonferroni') %>%
                add_significance()
# stat.test

head(stat.test)

stat.test <- stat.test %>% add_xy_position(x = "cluster")

stat.test <- cellannoExpMelt %>% 
                mutate(variable2 = paste(cluster, gene, sep = "-"))%>% 
                group_by(cluster, gene)%>% 
                t_test(formula = expression ~ group)%>%
                add_xy_position() %>%
                adjust_pvalue(method = 'bonferroni') %>%
                add_significance()

stat.test <- as.data.frame(stat.test)

rownames(stat.test) <- paste(stat.test$cluster, stat.test$gene, sep = "-")

head(stat.test)

stat.test$group <- "Control"



options(repr.plot.width=10, repr.plot.height=6)
ggplot(cellannoExpMelt2, aes(x=group, y=scale, fill=group)) +
    facet_grid(cluster~gene, switch="y", scales = "free") +
    scale_y_continuous(position="right", limits = c(-0.2, 1.2)) + # c(-0.7, 1.5)
    # geom_jitter(size=0.1) +
    geom_violin(trim = F) +
    geom_boxplot(outlier.size = 0, width=0.3) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = "", y = "", title = "") +
    stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.35, label.y = 1.1) +
    theme(strip.background = element_rect(fill = "gray97", color = NA), 
            strip.placement = "outside", 
            strip.text.x = element_text(face="italic", size = 10),
            strip.text.y = element_text(face="bold", size = 11)) +
    theme(axis.ticks.x = element_blank(), axis.ticks = element_line(size = 0.1), 
            axis.text.x  = element_text(face="plain", angle=70, size = 10, color = "black", vjust=0.5),
            axis.text.y  = element_text(face="plain", size = 8, color = "black"),
            axis.title = element_text(size = 12)) +
    scale_fill_manual(values=c('white','skyblue')) +
    geom_text(data = dat_text, x = 1.5, y = 0.95, size = 3.5, mapping = aes(label = label),
              color = ifelse(as.numeric(dat_text$log2FC) < 0,'blue','red')) 



rownames(dat_text) <- paste(dat_text$cluster, dat_text$gene, sep = "-")

dat_text <- cbind(dat_text, stat.test[rownames(dat_text),c("p","p.adj","p.adj.signif")])

head(dat_text)

# tmp.genes <- c("Ptn","Ncl","Ncam1","Sdc3","Sdc1","Ptprz1","L1cam","Fgfr1","Cadm3")

# dat_text

options(repr.plot.width=10, repr.plot.height=6)
ggplot(cellannoExpMelt2, aes(x=group, y=scale, fill=group)) +
    facet_grid(cluster~gene, switch="y", scales = "free") +
    scale_y_continuous(position="right", limits = c(-0.2, 1.2)) + # c(-0.7, 1.5)
    # geom_jitter(size=0.1) +
    geom_violin(trim = F) +
    geom_boxplot(outlier.size = 0, width=0.3) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = "", y = "", title = "") +
    # stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01) +
    theme(strip.background = element_rect(fill = "gray97", color = NA), 
            strip.placement = "outside", 
            strip.text.x = element_text(face="italic", size = 10),
            strip.text.y = element_text(face="bold", size = 11)) +
    theme(axis.ticks.x = element_blank(), axis.ticks = element_line(size = 0.1), 
            axis.text.x  = element_text(face="plain", angle=70, size = 10, color = "black", vjust=0.5),
            axis.text.y  = element_text(face="plain", size = 8, color = "black"),
            axis.title = element_text(size = 12)) +
    scale_fill_manual(values=c('white','skyblue')) +
    geom_text(data = dat_text, x = 1.5, y = 0.95, size = 3.5, mapping = aes(label = label),
              color = ifelse(as.numeric(dat_text$log2FC) < 0,'blue','red')) +
    geom_text(data = dat_text, x = 1.5, y = 1.1, size = 3.5, mapping = aes(label = p.adj.signif),
              color = ifelse(as.numeric(dat_text$log2FC) < 0,'blue','red')) 

ggsave(filename = "Figures/violin.key.genes.pdf", width = 10, height = 6, dpi = 800)



plot.violin.barplot.indep.scale <- function(tmp.df, tmp.dat_text, tmp.gene) {
    #
    tmp.df <- tmp.df %>% filter(scale < quantile(tmp.df$scale, 0.95))
    max.y <- max(tmp.df$scale)
    #
    p <- ggplot(tmp.df, aes(x=group, y=scale, fill=group)) +
            # facet_grid(cluster~gene, switch="y", scales = "free") +
            facet_wrap(.~cluster, ncol = 1) +
            scale_y_continuous(position="right", limits = c(0, max.y*1.1)) +
            # geom_jitter(size=0.1) +
            geom_violin() +
            geom_boxplot(outlier.size = 0, width=0.2) +
            theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            labs(x = "", y = "", title = tmp.gene) +
            stat_compare_means(label = "p.signif", method = "wilcox.test", label.x = 1.4, label.y = max.y*0.85) +
            theme(strip.background = element_rect(fill = "gray97", color = NA), 
                    strip.placement = "outside", strip.text.x = element_text(face="plain", size = 10),
                    strip.text.y = element_text(face="bold", size = 11)) +
            theme(plot.title = element_text(hjust = 0.5, size = 12, face="bold.italic")) +
            theme(axis.ticks.x = element_blank(), axis.ticks = element_line(size = 0.1), 
                    axis.text.x  = element_text(face="plain", angle=70, size = 10, color = "black", vjust=0.5),
                    axis.text.y  = element_text(face="plain", size = 6, color = "black"),
                    axis.title = element_text(size = 12)) +
            theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "pt")) + # top, right, bottom, left
            scale_fill_manual(values=c('white','skyblue')) +
            geom_text(data = tmp.dat_text, x = 1.5, y = max.y*1.05, size = 3, mapping = aes(label = label),
                      color = ifelse(as.numeric(tmp.dat_text$log2FC) < 0,'blue','red')) 
    return(p)
}

p.list <- list()
for (tmp.gene in unique(cellannoExpMelt$gene)) {
    print(tmp.gene)
    tmp.df <- subset(cellannoExpMelt, gene==tmp.gene)
    tmp.dat_text <- subset(dat_text, gene==tmp.gene)
    tmp.p <- plot.violin.barplot.indep.scale(tmp.df, tmp.dat_text, tmp.gene)
    p.list[[tmp.gene]] <- tmp.p
}

options(repr.plot.width=11, repr.plot.height=8)
cowplot::plot_grid(plotlist = p.list, nrow = 1)


