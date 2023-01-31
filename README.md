# Vcl-paper
source codes for Vcl paper

unpublished dataset (Seurat Object) for Vcl ENCC and CNCC project
Zenodo link: https://doi.org/10.5281/zenodo.7470702

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
