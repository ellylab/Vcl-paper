# Vcl CNCC project
Data and source codes for Vcl CNCC project

# R script and code [under Vcl-CNCC-project folder]
- 2021_Wen_Chen_heart_CNCC_integration.html: code for [public CNCC single-cell RNA-seq data](https://www.embopress.org/doi/full/10.15252/embr.202152389) integration **[Figure 5A-E]**
- Vcl_CNCC_scRNA-seq_analysis.html: code for our single-cell RNA-seq data analysis **[Figure 5F-M]**
- Vcl-ST-Seurat-visium.ipynb and Vcl_Spatial_analysis.R: Spatial transcriptome data analysis **[Figure 7A-J]**

# R processed data
Zenodo link
- https://doi.org/10.5281/zenodo.8202751 [processed ST data]
- https://doi.org/10.5281/zenodo.7629719 [processed scRNA-seq data]
        
        
Below are the scripts showing the labels of cells in two datasets.
```
#######################
library(Seurat)

# seurat object for vcl ENCC
print(load("keyRdata/vcl.encc.integrated.Rdata"))

# annotation subtypes in each lineage
table(vcl.combined$lineage.sub)

# condition, normal vs Vcl cKO ENCC
table(vcl.combined$stage)

#######################
# seurat object for vcl CNCC
print(load("E13.5_CNCC_merged_updated.Rdata"))

# sample, ct2 is control and vcl is Vcl cKO CNCC
table(seuset$orig.ident)

# cell annotation file
# unsupervised clusters, see more information in HKU thesis
table(all_tsne$ident)
```

```
#######################
# seurat object for public CNCC from E10.5 - P7
print(load("all.cncc.combined.EMBO.mapped.Rdata"))

# annotated cell types based on markers from EMBO paper
table(cncc.combined.sub@active.ident)

# stages of CNCC cells
table(cncc.combined.sub$stage)
```
Data source (fastq): https://www.embopress.org/doi/full/10.15252/embr.202152389
        
        
        
        
           
#######################
## spatial data of Vcl (with GFP)
- Vcl-ST-Seurat-visium.ipynb
- data in zenodo database

# Raw fastqs in Sequence Read Archive (SRA)
- PRJNAXXXXX (scRNA-seq fastqs)
- PRJNA1004428 (ST raw fastqs)
        
# Contact
- ellyngan engan@hku.hk
- Mingxuan LIANG u3608243@connect.hku.hk
- Zhixin Li zxlee@hku.hk
