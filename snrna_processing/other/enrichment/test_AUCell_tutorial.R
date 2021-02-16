
# install -----------------------------------------------------------------
# To support paralell execution:
# BiocManager::install(c("AUCell", "doMC", "doRNG"))
# For the main example:
# BiocManager::install(c("mixtools", "GEOquery", "SummarizedExperiment"))
# For the examples in the follow-up section of the tutorial:
# BiocManager::install(c("DT", "plotly", "NMF", "d3heatmap", "shiny", "rbokeh",
#                        "dynamicTreeCut","R2HTML","Rtsne", "zoo"))
packages <- c("AUCell", "doMC", "doRNG",
              "mixtools", "GEOquery", "SummarizedExperiment",
              "DT", "plotly", "NMF", "d3heatmap", "shiny", "rbokeh",
              "dynamicTreeCut","R2HTML","Rtsne", "zoo")
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    BiocManager::install(pkgs = pkg_name_tmp)
  }
}
## package ‘d3heatmap’ is not available (for R version 3.6.1) 

# 0. Load scRNA-seq dataset and gene sets
# Working directory ---------------------------------
dir.create("AUCell_tutorial")
setwd("AUCell_tutorial") # or in the first code chunk (kntr options), if running as Notebook


# Expression matrix -------------------------------------------------------
# (This may take a few minutes)
library(GEOquery)
attemptsLeft <- 20
while(attemptsLeft>0)
{
  geoFile <- tryCatch(getGEOSuppFiles("GSE60361", makeDirectory=FALSE), error=identity) 
  if(methods::is(geoFile,"error")) 
  {
    attemptsLeft <- attemptsLeft-1
    Sys.sleep(5)
  }
  else
    attemptsLeft <- 0
}

gzFile <- grep(".txt.gz", basename(rownames(geoFile)), fixed=TRUE, value=TRUE)
message(gzFile)
txtFile <- gsub(".gz", "", gzFile, fixed=TRUE)
message(txtFile)
gunzip(filename=gzFile, destname=txtFile, remove=TRUE)

library(data.table)
geoData <- fread(txtFile, sep="\t")
geneNames <- unname(unlist(geoData[,1, with=FALSE]))
exprMatrix <- as.matrix(geoData[,-1, with=FALSE])
rm(geoData)
dim(exprMatrix)
rownames(exprMatrix) <- geneNames
exprMatrix[1:5,1:4]

# Remove file
file.remove(txtFile)

# Save for future use
mouseBrainExprMatrix <- exprMatrix
save(mouseBrainExprMatrix, file="exprMatrix_AUCellVignette_MouseBrain.RData")

# load("exprMatrix_AUCellVignette_MouseBrain.RData")
set.seed(333)
exprMatrix <- mouseBrainExprMatrix[sample(rownames(mouseBrainExprMatrix), 5000),]
exprMatrix[1:5,1:4]


# Gene sets ---------------------------------------------------------------
library(AUCell)
library(GSEABase)
gmtFile <- paste(file.path(system.file('examples', package='AUCell')), "geneSignatures.gmt", sep="/")
geneSets <- getGmt(gmtFile)
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))
# Random
set.seed(321)
extraGeneSets <- c(
  GeneSet(sample(rownames(exprMatrix), 50), setName="Random (50g)"),
  GeneSet(sample(rownames(exprMatrix), 500), setName="Random (500g)"))

countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
# Housekeeping-like
extraGeneSets <- c(extraGeneSets,
                   GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK-like (100g)"))

geneSets <- GeneSetCollection(c(geneSets,extraGeneSets))
names(geneSets)

# 1. Build gene-expression rankings for each cell -------------------------
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)
save(cells_rankings, file="cells_rankings.RData")

# 2. Calculate enrichment for the gene signatures (AUC) -------------------
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file="cells_AUC.RData")

# 3. Determine the cells with the given gene signatures or active  --------
set.seed(123)
par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]
cells_assignment[[1]]$aucThr
cells_assignment[[1]]$assignment
oligodencrocytesAssigned <- cells_assignment$Oligodendrocyte_Cahoy$assignment

cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)

library(DT)
datatable(assignmentTable, options = list(pageLength = 10), filter="top")

# Explore cells/clusters based on the signature score ---------------------
load(paste(file.path(system.file('examples', package='AUCell')), "cellsTsne.RData", sep="/"))
cellsTsne <- cellsTsne$Y
plot(cellsTsne, pch=16, cex=.3)

selectedThresholds <- getThresholdSelected(cells_assignment)

test_obj <- getAUC(cells_AUC)
nBreaks <- 5 # Number of levels in the color palettes
# Color palette for the cells that do not pass the threshold
colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
# Color palette for the cells that pass the threshold
colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)

passThreshold <- getAUC(cells_AUC)[1,] >  selectedThresholds[1]
aucSplit <- split(getAUC(cells_AUC)[1,], passThreshold)
cut(aucSplit[[1]], breaks=nBreaks)
colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)]


