#read in data matrix
rnaMat = readRDS("ADCC_RNA_Matrix_cleaned.rds")

#load seurat
library(Seurat)
seu <- CreateSeuratObject(rnaMat, min.cells = 3, min.features = 200)

#filter and preprocess data
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

options(repr.plot.width=12, repr.plot.height=6)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

options(repr.plot.width=9, repr.plot.height=6)
ggplot(seu@meta.data, aes(nCount_RNA, nFeature_RNA)) +
  geom_hex(bins = 100) +
  scale_fill_scico(palette = "devon", direction = -1, end = 0.9) +
  scale_x_log10(breaks = breaks_log(12)) + 
  scale_y_log10(breaks = breaks_log(12)) + annotation_logticks() +
  labs(x = "Total UMI counts", y = "Number of genes detected") +
  theme(panel.grid.minor = element_blank())

ggplot(seu@meta.data, aes(nCount_RNA, percent.mt)) +
  geom_pointdensity() +
  scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
  labs(x = "Total UMI counts", y = "Percentage mitochondrial")


seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 9) 
length(grep("A431R", colnames(seu)))
length(grep("A431S", colnames(seu)))
length(grep("FaDuR", colnames(seu)))
length(grep("FaDuS", colnames(seu)))
length(grep("SKOV3R", colnames(seu)))
length(grep("SKOV3S", colnames(seu)))


seu = NormalizeData(seu) 
cells = colnames(seu)
cells = unlist(lapply(cells, function(x){stringr::str_split(x, stringr::coll("."))[[1]][1]}))
seu = AddMetaData(seu, factor(cells,
                              levels = c("A431S", "A431R", "SKOV3S", "SKOV3R", "FaDuS",
                                         "FaDuR")), col.name = 'letter.idents2')
seu = ScaleData(seu)

#save seurat object
saveRDS(seu, "seuratObject1.rds")

#plot expression of EGFR and HER2 (ERBB2)
VlnPlot(seu, c("EGFR", "ERBB2"), group.by = 'letter.idents2', pt.size = 0)

#differential expression between sensitive and resistant cells
deFaDu = FindMarkers(seu, ident.1 = "FaDuR", ident.2 = "FaDuS", group.by = 'letter.idents')
deA431 = FindMarkers(seu, ident.1 = "A431R", ident.2 = "A431S", group.by = 'letter.idents')
deSKOV3 = FindMarkers(seu, ident.1 = "SKOV3R", ident.2 = "SKOV3S", group.by = 'letter.idents')

fullResultsFaDu = FindMarkers(seu, ident.1 = "FaDuR", ident.2 = "FaDuS", group.by = 'letter.idents', logfc.threshold = 0, min.pct = 0, min.diff.pct = 0)
fullResultsA431 = FindMarkers(seu, ident.1 = "A431R", ident.2 = "A431S", group.by = 'letter.idents', logfc.threshold = 0, min.pct = 0, min.diff.pct = 0)
fullResultsSKOV3 = FindMarkers(seu, ident.1 = "SKOV3R", ident.2 = "SKOV3S", group.by = 'letter.idents', logfc.threshold = 0, min.pct = 0, min.diff.pct = 0)

saveRDS(fullResultsFaDu, "fullDEResultsFaDu.rds")
saveRDS(fullResultsA431, "fullDEResultsA431.rds")
saveRDS(fullResultsSKOV3, "fullDEResultsSKOV3.rds")

#differential expression between resistant cell lines
deA431FaDu = FindMarkers(seu, ident.1 = "A431R", ident.2 = "FaDuR", group.by = 'letter.idents')
deA431SKOV3 = FindMarkers(seu, ident.1 = "A431R", ident.2 = "SKOV3R", group.by = 'letter.idents')
deSKOV3FaDu = FindMarkers(seu, ident.1 = "SKOV3R", ident.2 = "FaDuR", group.by = 'letter.idents')

deList = list(deFaDu,deA431, deSKOV3, deA431FaDu, deA431SKOV3, deSKOV3FaDu)
names(deList) = c("FaDuR v FaDuS", "A431R v A431S", "SKOv3R v SKOV3S", "A431R v FaDuR", "A431R v SKOV3R", "SKOv3R v FaDuR")
openxlsx::write.xlsx(deList, "cellLineDEComps.xlsx")


#volcano plots
EnhancedVolcano::EnhancedVolcano(fullResultsFaDu, lab = rownames(fullResultsFaDu), x = 'avg_logFC', y = 'p_val_adj', FCcutoff = 0.25, title = NULL, subtitle = NULL)
EnhancedVolcano::EnhancedVolcano(fullResultsA431, lab = rownames(fullResultsA431), x = 'avg_logFC', y = 'p_val_adj', FCcutoff = 0.25, , title = NULL, subtitle = NULL)
EnhancedVolcano::EnhancedVolcano(fullResultsSKOV3, lab = rownames(fullResultsSKOV3), x = 'avg_logFC', y = 'p_val_adj', FCcutoff = 0.25, , title = NULL, subtitle = NULL)

pdf("cellLineCompVolcanoPlots.pdf", width = 8, height = 8)
  EnhancedVolcano::EnhancedVolcano(deFaDu, lab = rownames(deFaDu), x = 'avg_logFC', y = 'p_val_adj', FCcutoff = 0.5, , title = "FaDuR v FaDuS", subtitle = NULL)
  EnhancedVolcano::EnhancedVolcano(deA431, lab = rownames(deA431), x = 'avg_logFC', y = 'p_val_adj', FCcutoff = 0.5, , title = "A431R v A431S", subtitle = NULL)
  EnhancedVolcano::EnhancedVolcano(deSKOV3, lab = rownames(deSKOV3), x = 'avg_logFC', y = 'p_val_adj', FCcutoff = 0.5, , title = "SKOV3R v SKOV3S", subtitle = NULL)
  EnhancedVolcano::EnhancedVolcano(deA431FaDu, lab = rownames(deA431FaDu), x = 'avg_logFC', y = 'p_val_adj', FCcutoff = 1, , title = "A431R v FaDuR", subtitle = NULL)
  EnhancedVolcano::EnhancedVolcano(deA431SKOV3, lab = rownames(deA431SKOV3), x = 'avg_logFC', y = 'p_val_adj', FCcutoff = 1, , title = "A431R v SKOV3R", subtitle = NULL)
  EnhancedVolcano::EnhancedVolcano(deSKOV3FaDu, lab = rownames(deSKOV3FaDu), x = 'avg_logFC', y = 'p_val_adj', FCcutoff = 1, , title = "SKOV3R v FaDuR", subtitle = NULL)
dev.off()


#find overlapping up and down expressed genes among all
#up in resistant, greater than 0.5 LFC
FaDuUpRes = rownames(deFaDu)[which(deFaDu$p_val_adj < 0.05 & deFaDu$avg_logFC > 0.5)] 

A431UpRes = rownames(deA431)[which(deA431$p_val_adj < 0.05 & deA431$avg_logFC > 0.5)] 

SKOV3UpRes = rownames(deSKOV3)[which(deSKOV3$p_val_adj < 0.05 & deSKOV3$avg_logFC > 0.5)]

#down in resistant, less than -0.5 LFC
FaDuDnRes = rownames(deFaDu)[which(deFaDu$p_val_adj < 0.05 & deFaDu$avg_logFC < -0.5)] 

A431DnRes = rownames(deA431)[which(deA431$p_val_adj < 0.05 & deA431$avg_logFC < -0.5)] 

SKOV3DnRes = rownames(deSKOV3)[which(deSKOV3$p_val_adj < 0.05 & deSKOV3$avg_logFC < -0.5)]

#up in all 3
upRes3 = intersect(FaDuUpRes,intersect(A431UpRes, SKOV3UpRes)) #none
#down in all 3
dnRes3 = intersect(FaDuDnRes,intersect(A431DnRes, SKOV3DnRes)) #none

#pairs
upResA431SKOV3 = intersect(A431UpRes, SKOV3UpRes)
upResA431FaDu = intersect(A431UpRes, FaDuUpRes)
upResFaDuSKOV3 = intersect(FaDuUpRes, SKOV3UpRes)

dnResA431SKOV3 = intersect(A431DnRes, SKOV3DnRes)
dnResA431FaDu = intersect(A431DnRes, FaDuDnRes)
dnResFaDuSKOV3 = intersect(FaDuDnRes, SKOV3DnRes)


#smaller ES cutoff
FaDuUpRes = rownames(deFaDu)[which(deFaDu$p_val_adj < 0.05 & deFaDu$avg_logFC > 0.25)] 

A431UpRes = rownames(deA431)[which(deA431$p_val_adj < 0.05 & deA431$avg_logFC > 0.25)] 

SKOV3UpRes = rownames(deSKOV3)[which(deSKOV3$p_val_adj < 0.05 & deSKOV3$avg_logFC > 0.25)]

#down in resistant, 
FaDuDnRes = rownames(deFaDu)[which(deFaDu$p_val_adj < 0.05 & deFaDu$avg_logFC < -0.25)] 

A431DnRes = rownames(deA431)[which(deA431$p_val_adj < 0.05 & deA431$avg_logFC < -0.25)] 

SKOV3DnRes = rownames(deSKOV3)[which(deSKOV3$p_val_adj < 0.05 & deSKOV3$avg_logFC < -0.25)]

#up in all 3
upRes3 = intersect(FaDuUpRes,intersect(A431UpRes, SKOV3UpRes)) # 5 out of 6 are in immune system GO term and in innate immune system reactome pathway

seu = AddMetaData(seu, factor(samplesSub, levels = levels(as.factor(samplesSub))[c(2,1,6,5,4,3)]), col.name = 'letter.idents2')

#plot for each
VlnPlot(seu, features = upRes3, group.by = 'letter.idents2', pt.size = 0)

DotPlot(seu, features = upRes3, group.by = 'letter.idents') + RotatedAxis()



#down in all 3
dnRes3 = intersect(FaDuDnRes,intersect(A431DnRes, SKOV3DnRes)) 

VlnPlot(seu, features = dnRes3, group.by = 'letter.idents2', pt.size = 0)


#by pairs
upResA431SKOV3 = intersect(A431UpRes, SKOV3UpRes)
upResA431FaDu = intersect(A431UpRes, FaDuUpRes)
upResFaDuSKOV3 = intersect(FaDuUpRes, SKOV3UpRes)

dnResA431SKOV3 = intersect(A431DnRes, SKOV3DnRes)
dnResA431FaDu = intersect(A431DnRes, FaDuDnRes)
dnResFaDuSKOV3 = intersect(FaDuDnRes, SKOV3DnRes)

#no ES cutoff
FaDuResNC = rownames(deFaDu)[which(deFaDu$p_val_adj < 0.05)] 

A431ResNC = rownames(deA431)[which(deA431$p_val_adj < 0.05)] 

SKOV3ResNC = rownames(deSKOV3)[which(deSKOV3$p_val_adj < 0.05)]


de3NC = intersect(FaDuResNC,intersect(A431ResNC, SKOV3ResNC))

saveRDS(de3NC, "DEGenesAllCellLines_NoESCutoff.rds")

#
deProteins = c("CD44", "EGFR", "F3", "ITGA6")

deFaDu[which(rownames(deFaDu) %in% deProteins),]

deA431[which(rownames(deA431) %in% deProteins),]

deSKOV3[which(rownames(deSKOV3) %in% deProteins),]

VlnPlot(seu, features = deProteins, group.by = 'letter.idents2', pt.size = 0)


#cd39 and cd73
fullResultsSKOV3[which(rownames(fullResultsSKOV3) %in% c("ENTPD1", "NT5E")),]


#gene set overrepresentation
source("geneORA.R")
oraList = vector("list", 12)
for(i in 1:6){
  oraList[[i]] = geneORA(rownames(deList[[i]][which(deList[[i]]$avg_logFC > 0),]), significance_threshold = 0.05, species = "human")
}

for(i in 7:12){
  oraList[[i]] = geneORA(rownames(deList[[i-6]][which(deList[[i-6]]$avg_logFC < 0),]), significance_threshold = 0.05, species = "human")
}

oraList = lapply(oraList, sort)
#convert to data frames
oraList2 = lapply(oraList, function(x){
  col1 = names(x)
  df = cbind(col1, as.numeric(unlist(x, use.names = F)))
  colnames(df) = c("Pathway", "p-value")
  df
})

names(oraList2) = c("FaDuR_up v FaDuS", "A431R_up v A431S", "SKOv3R_up v SKOV3S",
                    "A431R_up v FaDuR", "A431R_up v SKOV3R", "SKOv3R_up v FaDuR",
                    "FaDuR_down v FaDuS", "A431R_down v A431S", "SKOv3R_down v SKOV3S",
                    "A431R_down v FaDuR", "A431R_down v SKOV3R", "SKOv3R_down v FaDuR")

openxlsx::write.xlsx(oraList2, "cellLinePathwayORA.xlsx")

