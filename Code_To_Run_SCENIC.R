ls()
rm(list=ls())
ls()
memory.size() ### Checking your memory size
memory.limit() ## Checking the set limit
memory.limit(size=56000) ### expanding your memory
memory.limit()
library(SeuratDisk)
library(SeuratData)
library(Seurat)

## Load in H5 Seurat File ##
Male_Seurat1 <- LoadH5Seurat("E:/MaleHigh_Social_January_19_2024/MaleHigh.h5Seurat")
Male_Seurat1

## Extract Log Normalized Counts From Seurat Object ##
counts <- Male_Seurat1@assays$RNA@data

## Look At First Few Rows Of Log Normalized Counts ##
head(counts)

## Get class of count object ##
class(counts)

## Convert counts object to a matrix called exprMatrix ##
exprMatrix <- as.matrix(counts)

## Examine exprMatrix object ##
head(exprMatrix) # first 6 rows
class(exprMatrix) # object class
is.atomic(exprMatrix) # True
dim(exprMatrix) # number of rows & cols
nrow(exprMatrix) # number of rows
ncol(exprMatrix) # number of cols
colnames(exprMatrix) # column names
exprMatrix <- exprMatrix[unique(rownames(exprMatrix)),]
exprMatrix[1:5,1:4]
dim(exprMatrix) # of rows & cols

## Cell Info/Phenodata, Should Be in Format Required by SCENIC (barcode, CellType, nGene, nUMI) ##
## Read In Phenodata Table as a csv file (extract needed info from Seurat metadata table) ##
cellInfo <- read.table(file.choose(), sep = ",", row.names = 1, header = TRUE)

## Cell info/phenodata ##
head(cellInfo)
cellInfo$nGene <- colSums(exprMatrix>0)
head(cellInfo)
dir.create("cellinfo")
saveRDS(cellInfo, file="E:/MaleHigh_Social_January_19_2024/cellinfo/cellInfo.Rds")

## Assign colors to variables ##
colVars <- list(CellType=c("Microglia"="forestgreen", 
                           "Mature Oligos"="darkorange", 
                           "L5/6 NP CTX"="magenta4", 
                           "Astrocytes"="hotpink", 
                           "L6b CT CTX"="red3", 
                           "VLMC"="skyblue", 
                           "Sncg GABA"="darkblue",
                           "L6 CT ENT"="grey50",
                           "Lamp5-Meis2-GABA"="coral1",
                           "Sst GABA"="gold2",
                           "OPCs"="olivedrab1",
                           "Pericytes"="moccasin",
                           "Meis2-Penk-GABA"="yellow"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
colVars
cellInfo
saveRDS(colVars, file="E:/MaleHigh_Social_January_19_2024/cellinfo/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

## Initialize SCENIC settings ##
library(SCENIC)
library(arrow)
library(feather)
org <- "mgi"
cisTarget_databases <- "/home/users/gd2/SCENIC3"
dbDir <- cisTarget_databases
myDatasetTitle <- "SCENIC Male High Social Mouse brain" # name for analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)

## Modify if needed ##
scenicOptions@inputDatasetInfo$cellInfo <- "E:/MaleHigh_Social_January_19_2024/cellinfo/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "E:/MaleHigh_Social_January_19_2024/cellinfo/colVars.Rds"

## Save For Use At A Later Time ##
saveRDS(scenicOptions, file="E:/MaleHigh_Social_January_19_2024/scenicOptions.Rds")

## Co-Expression Network ##
# (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMatrix, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMatrix),
                           minSamples=ncol(exprMatrix)*.01)
interestingGenes <- c("Plp1", "Sox10", "Mag")
interestingGenes[which(!interestingGenes %in% genesKept)]
exprMatrix_filtered <- exprMatrix[genesKept, ]
dim(exprMatrix_filtered)
write.table(exprMatrix_filtered, sep=",", file="exprMatrixFilt.csv", row.names=TRUE, col.names=NA, quote=FALSE)

## Correlation ##
runCorrelation(exprMatrix_filtered, scenicOptions)
exprMatrix_filtered[1:5,1:4]

# Run GENIE3 (10 steps, computationally intensive, 5,000 cells can take 4 to 5 days to run)
library(GENIE3)
runGenie3(exprMatrix_filtered, scenicOptions) ## won't run because no motifAnnotations_mgi so run next line then run again #
motifAnnotations_mgi <- motifAnnotations
runGenie3(exprMatrix_filtered, scenicOptions)

# Save GENIE3 Output #
saveRDS(scenicOptions, file="E:/MaleHigh_Social_January_19_2024/scenicOptions.Rds") 

# Run wrapper functions 1 thru 3 (runSCENIC_2 takes about 2-3 hours, runSCENIC_1 & runSCENIC_3 are quick) #
scenicOptions <- readRDS("E:/MaleHigh_Social_January_19_2024/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMatrix)

# Save current SCENIC status #
saveRDS(scenicOptions, file="E:/MaleHigh_Social_January_19_2024/scenicOptions.Rds")

## Prepare For Analysis Of Data ##
## Read in scenicOptions object ##
scenicOptions <- readRDS("scenicOptions.Rds")

## Enter colVars And cellInfo Data As A Variable ##
colVars <- readRDS("cellinfo/colVars.Rds")
head(colVars)
cellInfo <- readRDS("cellinfo/cellInfo.Rds")
head(cellInfo)

## Can Modify cellInfo And colVars Info In scenicOptions Object ##
scenicOptions@inputDatasetInfo$cellInfo <- "E:/MaleHigh_Social_January_19_2024/cellinfo/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "E:/MaleHigh_Social_January_19_2024/cellinfo/colVars.Rds"

# Check scenicOptions Stored Information #
getSettings(scenicOptions) # arguments for functions/steps
getDatasetInfo(scenicOptions) # info about dataset
getIntName(scenicOptions) # intermediate files
getDatabases(scenicOptions) # cis target databases used
getOutName(scenicOptions) # give output file names
tsneFileName(scenicOptions) # gives tsne default to use for plots and analysis
getDatasetInfo(scenicOptions, "datasetTitle") # give title of dataset SCENIC was run on
dbVersion(getDatabases(scenicOptions)) # gives version of cis target databases

## Prior to runSCENIC_4 I run the following steps. Then I run runSCENIC_4 ##
## This steps lets you select the tsne plot you think the cell clusters look the best in ##

## Here you pick the number of PCs (principal components) you want to examine in the tsne plots ##
nPcs <- c(5,15,50)

## set the seed for all of the tsne plots that will be examined ##
scenicOptions@settings$seed <- 123 

## create a variable called fileNames to store the file names of the different tsne plots that will be generated ##
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")

## Plot as pdf (individual files stored in int/): ##
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

## look at names of tsne files that were created ##
fileNames

## View Compare t-SNEs And Select Appropriate nPcs and Perpl For Use Later ##
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)
head(scenicOptions@inputDatasetInfo$cellInfo)

## Using only "high-confidence" regulons (normally similar to non-high-confidence tsnes generated above) ##
par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

## Pick nPcs And Perplexity To Use As Defaults For Future Plots. Enter Data Into scenicOptions Object And Save ##
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 15
scenicOptions@settings$defaultTsne$perpl <- 50

## To Run Binarize The Network Activity ##
## Read In exprMatrix And scenicOptions ##
exprMatrix <- read.table(file.choose(), sep = ",", row.names = 1, header = TRUE)

## Select Final AUC Thresholds ##
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMatrix)
savedSelections <- shiny::runApp(aucellApp)

## Save the modified thresholds ##
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
scenicOptions@settings$devType="png"

## Following threshold optimization, Run runSCENIC_4_aucell_binarize To Binarize Network ##
runSCENIC_4_aucell_binarize(scenicOptions, exprMat = exprMatrix)

## Save Output Generated To This Point ##
saveRDS(scenicOptions, file="scenicOptions.Rds") 

## Cell-type specific regulators ##
## (based on the Regulon Specificity Score (RSS) proposed by Suo et al. ##
## for the Mouse Cell Atlas in 2018). Useful for big analysis with many cell types, to identify the ##
## cell-type specific regulons. ##
library(AUCell)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"])
head(rss)
data <- as.data.frame(rss)
write.table(data, sep=",", file="Rss_Per_Cluster.csv", row.names=TRUE, col.names=NA, quote=FALSE)
rssPlot <- plotRSS(rss)
## Heatmap with regulons specific to cell types ##
plotly::ggplotly(rssPlot$plot)
## Plot showing top 5 regulons in each cell type ##
plotRSS_oneSet(rss, setName = "Mature Oligos")


## Now look at Average Regulon Activity Per Cluster ##
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled
data2 <- as.data.frame(regulonActivity_byCellType_Scaled)
write.table(data2, sep=",", file="Average Regulon Activity by Cell type.csv", row.names=TRUE, col.names=NA, quote=FALSE)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

## This looks at Top regulators per Cluster ##
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
meta2 <- as.data.frame(topRegulators)
write.table(meta2, sep=",", file="TopRegulators_Per_Cluster.csv", row.names=TRUE, col.names=NA, quote=FALSE)

## Binarized Version Or % Of Cell With Regulons Active ##
minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
binary <- as.data.frame(binaryActPerc_subset)
head(binary)
write.table(binary, sep=",", file="Heatmap_Cell_Binarized.csv", row.names=TRUE, col.names=NA, quote=FALSE)
binary2 <- read.table(file.choose(), sep = ",", row.names = 1, header = TRUE)
ComplexHeatmap::Heatmap(binary, name="Regulon activity (%)", col = c("white","red"))

## Top Regulators For Binary Examination ##
topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
viewTable(topRegulators)
meta3 <- as.data.frame(topRegulators)
write.table(meta3, sep=",", file="TopRegulators_BinaryActive_Per_Cluster.csv", row.names=TRUE, col.names=NA, quote=FALSE)

## Analysis Step 3 Graphing ##
## First Graph of AUC and Binary Activity In Each Cell Type ##
scenicOptions@settings$seed <- 123 # same seed for all of them
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=15, perpl=50)
fileNames <- tsneAUC(scenicOptions, aucType="Binary", nPcs=15, perpl=50)
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=15, perpl=50, onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
fileNames <- tsneAUC(scenicOptions, aucType="Binary", nPcs=15, perpl=50, onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=1)

## Density plot to detect most likely stable states (higher-density areas in the t-SNE) ##
## Get Density Plot For AUC Regulon ##
library(KernSmooth)
library(RColorBrewer)
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

## Get Density Plot For Binary Regulon ##
scenicOptions@settings$defaultTsne$aucType <- "binary"
scenicOptions@settings$defaultTsne$dims <- 15
scenicOptions@settings$defaultTsne$perpl <- 50
saveRDS(scenicOptions, file="scenicOptions.Rds")
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

## Project AUC and Transcription Factors on to tSNEs ##
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
tfs <- c("Olig1","Dbx2_extended", "Mafb", "Esrrg", "Sox10", "Fli1")
dev.off()
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "binaryAUC")

## Genes In The Regulons ##
regulons <- loadInt(scenicOptions, "regulons")
regulons
x <- regulons[c("Irf8")]
x
class(x)
regulon1 <- as.data.frame(x)
write.table(regulon1, sep=",", file="Irf8_Genes.csv", row.names=TRUE, col.names=NA, quote=FALSE)

## Details On TF targets ##
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Sox8" & highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=5))

## TF motifs ##
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Nkx6-2"]
viewMotifs(tableSubset)

## This Works #### View Regulon Activities In Specific Seurat Clusters ##
Male_Seurat1 <- LoadH5Seurat("E:/MaleHigh_Social_January_19_2024/MaleHigh.h5Seurat")
Male_Seurat1
dr_coords <- Embeddings(Male_Seurat1, reduction="umap")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
tfs <- c("Sox2","Elf2_extended","Olig2_extended", "Nkx6-2_extended", "Sox8", "Nfe2i3")
dev.off()
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")
AUCell::AUCell_plotTSNE(dr_coords, exprMat=NULL,
                        cellsAUC=selectRegulons(regulonAUC, tfs), 
                        plots = c("histogram", "binaryAUC", "AUC"))





