inputFiles = c("/RErb-A431-R-ATAC/outs/fragments.tsv.gz",
               "/RErb-A431-S-ATAC/outs/fragments.tsv.gz",
               "/RErb-FaDu-R-ATAC/outs/fragments.tsv.gz",
               "/RErb-FaDu-S-ATAC/outs/fragments.tsv.gz",
               "/RErb-SKOV3-R-ATAC/outs/fragments.tsv.gz",
               "/RErb-SKOV3-S-ATAC/outs/fragments.tsv.gz")

names(inputFiles) = c("A431R", "A431S", "FaDuR", "FaDuS", "SKOV3R", "SKOV3S")

library(ArchR)
addArchRGenome("hg19")

addArchRThreads(threads = 16) 

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 1, #Dont set this too high because you can always increase later
  filterFrags = 2000,
  maxFrags = 1e8,
  addTileMat = FALSE,
  addGeneScoreMat = TRUE,
  force=TRUE
)

saveRDS(ArrowFiles, "arrowFileNames.rds")