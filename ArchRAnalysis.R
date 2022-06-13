#read in atac data
library(ArchR)
set.seed(1) #setting seed for reproducibility

addArchRThreads(threads = 1) 

addArchRGenome("hg19")

ArrowFiles = readRDS("arrowFileNames.rds")

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ADCC",
  copyArrows = TRUE 
)


p1 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1

p2 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p2

plotPDF(p1,p2, name = "QC-Sample-Statistics.pdf", ArchRProj = proj, addDOC = FALSE, width = 6, height = 8)


#subset to a slightly higher TSS enrichment
projSub = proj[which(proj$TSSEnrichment > 1.3)]

#latent sematic indexing dimensionality reduction
projSub = addIterativeLSI(ArchRProj = projSub, useMatrix = "GeneScoreMatrix", name = "IterativeLSI",
                          clusterParams = list(resolution = c(0.2), sampleCells = 10000, maxClusters = 30, n.start
                                               = 10), varFeatures = 10000, dimsToUse = 1:30, iterations =3, force=T)

#seurat style clustering
projSub = addClusters(input = projSub, reducedDims = "IterativeLSI", force=T)


#embed umap
projSub = addUMAP(ArchRProj = projSub, reducedDims = "IterativeLSI", minDist = 0.0001, force=T, threads = 1, n_threads=1)
#plot bulk vs CD45pos vs CD45neg on umap
p = plotEmbedding(ArchRProj = projSub, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

#plot cluster on umap
p1 = plotEmbedding(ArchRProj = projSub, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

ggAlignPlots(p, p1, type = "h")

#save to pdf
plotPDF(p,p1, name = "Plot-UMAP-Clusters.pdf",
        ArchRProj = projSub, addDOC = FALSE, width = 5, height = 5)

p2 = plotEmbedding(ArchRProj = projSub, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")


#get marker genes of each sample
markersGS = getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#write marker region/gene list
markerList = getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")


openxlsx::write.xlsx(markerList, "geneMarkerListBySample.xlsx")

#save proj
saveRDS(proj, "fullADCCArchRProj.rds")
saveRDS(projSub, "TSSCutADCCArchRProj.rds")




##differential accessibility between treatments 
g1List = list(c("A431R"), c("FaDuR"), c("SKOV3R"), c("A431R"), c("FaDuR"), c("SKOV3R"))
g2List = list(c("A431S"), c("FaDuS"), c("SKOV3S"), c("FaDuR", "SKOV3R"),
              c("A431R", "SKOV3R"), c("A431R", "FaDuR"))
outputList = vector("list", length(g1List)+1)
pv = vector("list", length(g1List)+1)

for(i in seq_along(g1List)){
  de = getMarkerFeatures(
    ArchRProj = projSub, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Sample",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    useGroups = g1List[[i]],
    bgdGroups = g2List[[i]]
  )
  
  pv[[i]] <- markerPlot(seMarker = de, name = g1List[[i]], cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
  deList = getMarkers(de, cutOff = "FDR <= 0.05")
  deOut = as.data.frame(deList[[1]])
  outputList[[i]] = deOut
}

#get all res vs all sens
projSub$Resistant = projSub$Sample
projSub$Resistant[which(projSub$Resistant %in% c("A431R", "FaDuR", "SKOV3R"))] = "Res"
projSub$Resistant[which(projSub$Resistant %in% c("A431S", "FaDuS", "SKOV3S"))] = "Sens"


de = getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Resistant",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "Res",
  bgdGroups = "Sens"
)

pv[[length(pv)]] <- markerPlot(seMarker = de, name = "Res", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
deList = getMarkers(de, cutOff = "FDR <= 0.05")
deOut = as.data.frame(deList[[1]])
outputList[[length(outputList)]] = deOut

library(openxlsx)
names(outputList) = c("A431R v A431S", "FaDuR v FaDuS", "SKOV3R v SKOV3S",
                      "A431R vs other Res", "FaDuR vs other Res",
                      "SKOV3R v other Res", "All Res v All Sens")
write.xlsx(outputList, file = "differentialAccessibilityADCC.xlsx")
plotPDF(pv, name = "VolcanoPlotsADCCATACdata.pdf", width = 5, height = 5, ArchRProj = projSub, addDOC = FALSE)

#are the 6 genes we found to be consistent in Res v Sens in RNA in the ATAC data?
outputList[[7]][which(outputList[[7]]$name %in% c("ISG15", "S100A11", "ANXA2", "FTH1", "IGFBP6", "PSMB9")),]

#CD39 and CD73 in the SKOV3R vs S
outputList[[3]][which(outputList[[3]]$name %in% c("ENTPD1", "NT5E")),]


#find overlap with DE genes
A431DE = read.xlsx("cellLineDEComps.xlsx", sheet = 2)
FaDuDE = read.xlsx("cellLineDEComps.xlsx", sheet = 1)
SKOV3DE = read.xlsx("cellLineDEComps.xlsx", sheet = 3)


length(intersect(A431DE[,1], outputList[[1]]$name))/dim(A431DE)[1]
length(intersect(FaDuDE[,1], outputList[[2]]$name))/dim(FaDuDE)[1]
length(intersect(SKOV3DE[,1], outputList[[3]]$name))/dim(SKOV3DE)[1]

upRegIntA431 = intersect(A431DE[which(A431DE$avg_logFC > 0),1], outputList[[1]]$name[which(outputList[[1]]$Log2FC > 0)])
upRegIntFaDu = intersect(FaDuDE[which(FaDuDE$avg_logFC > 0),1], outputList[[2]]$name[which(outputList[[2]]$Log2FC > 0)])
upRegIntSKOV3 = intersect(SKOV3DE[which(SKOV3DE$avg_logFC > 0),1], outputList[[3]]$name[which(outputList[[3]]$Log2FC > 0)])

dnRegIntA431 = intersect(A431DE[which(A431DE$avg_logFC < 0),1], outputList[[1]]$name[which(outputList[[1]]$Log2FC < 0)])
dnRegIntFaDu = intersect(FaDuDE[which(FaDuDE$avg_logFC < 0),1], outputList[[2]]$name[which(outputList[[2]]$Log2FC < 0)])
dnRegIntSKOV3 = intersect(SKOV3DE[which(SKOV3DE$avg_logFC < 0),1], outputList[[3]]$name[which(outputList[[3]]$Log2FC < 0)])

intList = list(upRegIntA431, upRegIntFaDu, upRegIntSKOV3, dnRegIntA431, dnRegIntFaDu, dnRegIntSKOV3)
names(intList) = c("upResBothA431", "upResBothFaDu", "upResBothSKOV3", "upSensBothA431", "upSensBothFaDu", "upSensBothSKOV3")
write.xlsx(intList, file = "differentialAccessibilityExpressionOverlapADCC.xlsx")

#do the intersecting gene lists overlap significantly with any GO pathways?
source("geneORA.R")


oraList = vector("list", length(intList))
for(i in seq_along(intList)){
  oraList[[i]] = geneORA(intList[[i]], significance_threshold = 0.05, species = "human")
}