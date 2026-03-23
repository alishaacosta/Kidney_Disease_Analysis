##### WGCNA ####

#load libaries
library(WGCNA)
library(limma)
library(igraph)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(writexl)
library(enrichR)
library(ggplot2)
NDKD = readRDS(file = "NDKD_processed.RDS")
DKD = readRDS(file = "DKD_processed.RDS")

#### NDKD dataset ####
#retrieve normalized expression matrix, change to gene symbol
exp_NDKD <- NDKD$v$E
#change entrez id to gene name, drop the rows that are unmatched
genes_NDKD <- data.frame('ENTREZID' = rownames(exp_NDKD))
rownames(genes_NDKD) <- genes_NDKD$ENTREZID
genes_NDKD$SYMBOL <- mapIds(x = org.Hs.eg.db,
                           keys = rownames(exp_NDKD),
                           column = 'SYMBOL',
                           keytype = 'ENTREZID',
                           multiVals = 'first')
genes_NDKD <- genes_NDKD[!is.na(genes_NDKD$SYMBOL),]

#subset expression by gene and change row names
exp_NDKD <- exp_NDKD[rownames(genes_NDKD),]
rownames(exp_NDKD) <- genes_NDKD$SYMBOL

#keep top 10,000 genes with highest mean values
gtk_NDKD <- rank(-rowMeans(exp_NDKD)) <= 10000
exp_NDKD <- exp_NDKD[gtk_NDKD,]

#transpose the matrix
exp_NDKD <- t(exp_NDKD)

#make sure gene choices are good 
gsg_NDKD <- goodSamplesGenes(exp_NDKD)
gsg_NDKD$allOK #returns true so we can move forward

#cluster sample to detect possible out liars
pdf(file = 'NDKD_hclust.pdf', width = 8, height = 5)
sampleTree_NDKD = hclust(dist(exp_NDKD), method = "average")
plot(sampleTree_NDKD, main = "NDKD Sample clustering to detect outliers", sub="", xlab="")
dev.off()

#scale independence 
powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft_NDKD = pickSoftThreshold(exp_NDKD, powerVector = powers) 

pdf(file = 'NDKD_scale_indp.pdf', width = 8, height = 8)
plot(sft_NDKD$fitIndices[,1], -sign(sft_NDKD$fitIndices[,3])*sft_NDKD$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_NDKD$fitIndices[,1], -sign(sft_NDKD$fitIndices[,3])*sft_NDKD$fitIndices[,2],
     labels=powers,cex=0.9,col="red");
#this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
dev.off()

#mean connectivity as a function of the soft-thresholding power
pdf(file = 'NDKD_mean_conn.pdf', width = 8, height = 8)
plot(sft_NDKD$fitIndices[,1], sft_NDKD$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity NDKD")) 
text(sft_NDKD$fitIndices[,1], sft_NDKD$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()

#a power of 10 was chosen
final_power = 10
cor <- WGCNA::cor
net_NDKD = blockwiseModules(exp_NDKD, 
                           power = final_power,
                           TOMType = "signed",
                           maxBlockSize = 10005, #one block
                           minModuleSize = 100, 
                           reassignThreshold = 0, #no reassigning genes
                           mergeCutHeight = 0.15,
                           deepSplit = 2,
                           numericLabels = FALSE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = 'NDKD',
                           verbose = 3)
cor<- stats::cor
save(net_NDKD, file = "net_NDKD_singleblock.RData")

table(net_NDKD$colors)
pdf(file = 'NDKD_Module_Dendrogram.pdf', width = 8, height = 5)
par(mfrow=c(1,1)) 
par(omi = c(0.3,0.3,0.3,0.3))
par(mai=c(1, 2, 1, 1)) #bottom, left, top, right
plotDendroAndColors(net_NDKD$dendrograms[[1]], net_NDKD$colors,
                    "Module colors NDKD",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#save a file of genes and their assigned module
genes_NDKD <- colnames(exp_NDKD)
modColors_NDKD <- net_NDKD[["colors"]]
modules_NDKD <- cbind(genes_NDKD, modColors_NDKD)
modules_NDKD <- as.data.frame(modules_NDKD)
write_xlsx(modules_NDKD, 'NDKD_modules.xlsx')

#calculate the module eigengenes 
MEs0_NDKD <- moduleEigengenes(exp_NDKD, net_NDKD$colors)$eigengenes
MEs_NDKD <- orderMEs(MEs0_NDKD)

pdf(file = 'NDKD_ME.pdf', width = 8, height = 5)
plotEigengeneNetworks(MEs_NDKD, colorLabels = TRUE,
                      setLabels = 'MEs for NDKD',
                      plotHeatmaps = FALSE,
) + abline(h = 0.15, col = "red")
dev.off()

#### DKD data set ####
exp_DKD <- DKD$v$E

genes_DKD <- data.frame('ENTREZID' = rownames(exp_DKD))
rownames(genes_DKD) <- genes_DKD$ENTREZID
genes_DKD$SYMBOL <- mapIds(x = org.Hs.eg.db,
                           keys = rownames(exp_DKD),
                           column = 'SYMBOL',
                           keytype = 'ENTREZID',
                           multiVals = 'first')
genes_DKD <- genes_DKD[!is.na(genes_DKD$SYMBOL) ,]

exp_DKD <- exp_DKD[rownames(genes_DKD),]
rownames(exp_DKD) <- genes_DKD$SYMBOL

gtk_DKD <- rank(-rowMeans(exp_DKD)) <= 10000
exp_DKD <- exp_DKD[gtk_DKD,]

exp_DKD <- t(exp_DKD)

gsg_DKD <- goodSamplesGenes(exp_DKD)
gsg_DKD$allOK 

sampleTree_DKD = hclust(dist(exp_DKD), method = "average")
pdf(file = 'DKD_hclust.pdf', width = 8, height = 5)
plot(sampleTree_DKD, main = "DKD Sample clustering to detect outliers", sub="", xlab="")
dev.off()

sft_DKD = pickSoftThreshold(exp_DKD, powerVector = powers) 
pdf(file = 'DKD_scale_indp.pdf', width = 8, height = 8)
plot(sft_DKD$fitIndices[,1], -sign(sft_DKD$fitIndices[,3])*sft_DKD$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_DKD$fitIndices[,1], -sign(sft_DKD$fitIndices[,3])*sft_DKD$fitIndices[,2],
     labels=powers,cex=0.9,col="red");
#this line corresponds to using an R^2 cut-off of h
abline(h=0.7,col="red")
dev.off()

#mean connectivity as a function of the soft-thresholding power
pdf(file = 'DKD_mean_conn.pdf', width = 8, height = 8)
plot(sft_DKD$fitIndices[,1], sft_DKD$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft_DKD$fitIndices[,1], sft_DKD$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()

cor <- WGCNA::cor
net_DKD = blockwiseModules(exp_DKD, 
                           power = final_power, #same as NDKD
                           TOMType = "signed", 
                           maxBlockSize = 10005, 
                           minModuleSize = 100, 
                           reassignThreshold = 0, #no reassigning genes
                           deepSplit = 2,
                           mergeCutHeight = 0.15,
                           numericLabels = FALSE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = 'DKD',
                           verbose = 3)
cor<- stats::cor
save(net_DKD, file = "net_DKD_singleblock.RData")
table(net_DKD$colors)

pdf(file = 'DKD_Module_Dendrogram.pdf', width = 8, height = 5)
par(mfrow=c(1,1)) 
par(omi = c(0.3,0.3,0.3,0.3))
par(mai=c(1, 2, 1, 1)) #bottom, left, top, right
plotDendroAndColors(net_DKD$dendrograms[[1]], net_DKD$colors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#save a file of genes and their assigned module
genes_D <- colnames(exp_DKD)
modColors_D <- net_DKD[["colors"]]
modules_D <- cbind(genes_D, modColors_D)
modules_D <- as.data.frame(modules_D)
write_xlsx(modules_D, 'DKD_modules.xlsx')

#calculate the module eigengenes 
MEs0_D <- moduleEigengenes(exp_DKD, net_DKD$colors)$eigengenes
MEs_D <- orderMEs(MEs0_D)

pdf(file = 'DKD_ME.pdf', width = 8, height = 5)
plotEigengeneNetworks(MEs_D, colorLabels = TRUE,
                      setLabels = 'MEs for DKD',
                      plotHeatmaps = FALSE,
) + abline(h = 0.15, col = "red")
dev.off()


#### Module Correspondence ####

#transpose dendrogram of NDKD onto DKD to visualize the module preservation
pdf(file = 'NDKD_onto_DKD.pdf', width = 8, height = 5)
par(mfrow=c(1,1)) 
par(omi = c(0.3,0.3,0.3,0.3))
par(mai=c(1, 2, 1, 1)) #bottom, left, top, right
plotDendroAndColors(net_DKD$dendrograms[[1]], net_NDKD$colors,
                    "Module colors NDKD",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#transpose dendrogram of DKD onto NDKD to visualize the module preservation
pdf(file = 'DKD_onto_NDKD.pdf', width = 8, height = 5)
par(mfrow=c(1,1)) 
par(omi = c(0.3,0.3,0.3,0.3))
par(mai=c(1, 2, 1, 1)) #bottom, left, top, right
plotDendroAndColors(net_NDKD$dendrograms[[1]], net_DKD$colors,
                    "Module colors DKD",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#retrieve the module labels
NDKD_Mod_Labels <- substring(names(net_NDKD$MEs), 3) 
DKD_Mod_Labels <- substring(names(net_DKD$MEs), 3) 

#define the number of modules per data set
nNDKD <-length(NDKD_Mod_Labels) 
nDKD <- length(DKD_Mod_Labels)

#the variables modColors_NDKD and modColors_D contain each gene and its associated module
#get the common genes only + subset
common_genes <- intersect(names(modColors_D), names(modColors_NDKD)) #results in 8973 genes
#geneColor for subset
gene_Color_D<- modColors_D[common_genes]
gene_Color_NDKD<- modColors_NDKD[common_genes]

#sort by alphabetical order
gene_Color_D <- gene_Color_D[sort(names(gene_Color_D))]
gene_Color_NDKD <- gene_Color_NDKD[sort(names(gene_Color_NDKD))]

#ensure they are in the same order
table(names(gene_Color_D) == names(gene_Color_NDKD)) 

#initialize tables of p-values and of the corresponding counts, fill with 0 to begin
pTable = matrix(0, nrow = nNDKD, ncol = nDKD)
CountTbl = matrix(0, nrow = nNDKD, ncol = nDKD)

#pair wise comparisons
for (NDKDmod in 1:nNDKD)
  for (DKDmod in 1:nDKD)
  {
    NDKDMembers = (gene_Color_NDKD == NDKD_Mod_Labels[NDKDmod]);
    DKDMembers = (gene_Color_D == DKD_Mod_Labels[DKDmod]);
    pTable[NDKDmod, DKDmod] = -log10(fisher.test(NDKDMembers, DKDMembers, alternative = "greater")$p.value);
    CountTbl[NDKDmod, DKDmod] = sum(gene_Color_NDKD == NDKD_Mod_Labels[NDKDmod] & 
                                      gene_Color_D == DKD_Mod_Labels[DKDmod] 
                                      )
  }

#heatmap of pvalues and number of overlapping genes per module
NDKDModTotals = apply(CountTbl, 1, sum) 
DKDModTotals = apply(CountTbl, 2, sum) 

pdf(file = "NDKD vs DKD modules.pdf", width = 10, height = 7)
par(mfrow=c(1,1)) 
par(mai=c(1.75, 2.25, 1, 1)) #bottom, left, top, right
par(omi=c(0.25,0.25,0.25,0.25))
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", DKD_Mod_Labels),
               yLabels = paste(" ", NDKD_Mod_Labels),
               colorLabels = TRUE,
               xSymbols = paste("DKD ", DKD_Mod_Labels, ": ", DKDModTotals, sep=""),
               ySymbols = paste("NDKD ", NDKD_Mod_Labels, ": ", NDKDModTotals, sep=""),
               textMatrix = CountTbl,
               colors = blueWhiteRed(100)[50:100],
               main = "Correspondence of NDKD and DKD modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE,
               zlim = c(0,70))
dev.off()

#save pTable and CountTbl
pTable <- as.data.frame(pTable)
rownames(pTable) <- NDKD_Mod_Labels
colnames(pTable) <- DKD_Mod_Labels
write_xlsx(pTable, 'pTable.xlsx')

CountTbl <- as.data.frame(CountTbl)
rownames(CountTbl) <- NDKD_Mod_Labels
colnames(CountTbl) <- DKD_Mod_Labels
write_xlsx(pTable, 'CountTbl.xlsx')

#### Module preservation statistics between data sets ####
#Will use each dataset as the reference
#NDKD as reference 
multi_expr <- list(ndkd = list(data = exp_NDKD), dkd = list(data = exp_DKD))
multi_color <- list(ndkd = modColors_NDKD)

mod_preservation <- modulePreservation(multiData = multi_expr,
                                       multiColor = multi_color,
                                       dataIsExpr = TRUE,
                                       networkType = 'signed',
                                       referenceNetworks = 1, #this is the ndkd dataset
                                       maxGoldModuleSize = 1000,
                                       maxModuleSize = 3000,
                                       nPermutations = 100)
#took ~10mins
mod_stats_NDKDref <-mod_preservation$preservation$Z$ref.ndkd$inColumnsAlsoPresentIn.dkd
mod_stats_NDKDref <- mod_stats_NDKDref[order(-mod_stats_NDKDref$Zsummary.pres),]
write_xlsx(mod_stats_NDKDref, 'module_preservation_statistics_NDKDref.xlsx')

#DKD as reference
multi_color_D <- list(dkd = modColors_D)
mod_preservation_dkd <- modulePreservation(multiData = multi_expr,
                                           multiColor = multi_color_D,
                                           dataIsExpr = TRUE,
                                           networkType = 'signed',
                                           referenceNetworks = 2, #this is the dkd dataset
                                           maxGoldModuleSize = 1000,
                                           maxModuleSize = 3000,
                                           nPermutations = 100
)

mod_stats_DKDref <-mod_preservation_dkd$preservation$Z$ref.dkd$inColumnsAlsoPresentIn.ndkd
mod_stats_DKDref <- mod_stats_DKDref[order(-mod_stats_DKDref$Zsummary.pres),]
write_xlsx(mod_stats_DKDref, 'module_preservation_statistics_DKDref.xlsx')

#NDKD modules lightcyan and midnightblue do not have DKD counterparts
#potentially they are specific to NDKD?

#### Module Trait relationship ####

ndkd_meta <- NDKD$v$targets
ndkd_meta <- ndkd_meta[,c(1,6,8)]
colnames(ndkd_meta)[2] <- 'age'
colnames(ndkd_meta)[3] <- 'sex'

ndkd_meta$sex <- ifelse(ndkd_meta$sex == 'Male', 0, 1)
ndkd_meta$age <- as.numeric(ndkd_meta$age)
ndkd_meta$group <- as.numeric(ndkd_meta$group)

#making sure rownames match
table(rownames(MEs_NDKD)==rownames(ndkd_meta))

#calculate pearson coefficent between traits and MEs
nSamples_NDKD = nrow(exp_NDKD)
moduleTraitCorr_NDKD <- cor(MEs_NDKD, ndkd_meta, use = 'p')
moduleTraitP_NDKD <- corPvalueStudent(moduleTraitCorr_NDKD, nSamples_NDKD)

text_matrix = paste(signif(moduleTraitCorr_NDKD,2), "\n(",
                    signif(moduleTraitP_NDKD, 1), ")", sep = "")

pdf(file = "Module-trait relationships NDKD.pdf", width = 10, height = 8)
par(mfrow=c(1,1))
par(mai=c(1,2,1,0)+0.3) #bottom, left, top, right
par(omi=c(0.25,0.5,0.25,0))
labeledHeatmap(Matrix = moduleTraitCorr_NDKD,
               xLabels = colnames(ndkd_meta),
               yLabels = colnames(MEs_NDKD),
               ySymbols = colnames(MEs_NDKD),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_matrix,
               setStdMargins = TRUE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships NDKD"))
dev.off()
#Lightcyan has a significant negative correlation with disease state, 
#suggesting it could be protective of disease progression?

dkd_meta <- DKD$v$targets
dkd_meta <- dkd_meta[,c(1,6)]
dkd_meta$state <- ifelse(dkd_meta$state == 'Advanced', 1, 0)
dkd_meta$group <- as.numeric(dkd_meta$group)

#making sure rownames match
table(rownames(MEs_D)==rownames(dkd_meta))
nSamples_DKD = nrow(exp_DKD)
moduleTraitCorr_D <- cor(MEs_D, dkd_meta, use = 'p')
moduleTraitP_D <- corPvalueStudent(moduleTraitCorr_D, nSamples_DKD)

text_matrix2 = paste(signif(moduleTraitCorr_D,2), "\n(",
                     signif(moduleTraitP_D, 1), ")", sep = "")
pdf(file = "Module-trait relationships DKD.pdf", width = 8, height = 8)
par(mfrow=c(1,1))
par(mai=c(1,2,1,0)+0.3) #bottom, left, top, right
par(omi=c(0.25,0.5,0.25,0))
labeledHeatmap(Matrix = moduleTraitCorr_D,
               xLabels = colnames(dkd_meta),
               yLabels = colnames(MEs_D),
               ySymbols = colnames(MEs_D),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_matrix2,
               setStdMargins = TRUE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships DKD"))
dev.off()


#### Examining lightcyan NDKD module ####
LC_genes <- modules_NDKD[modules_NDKD$modColors_NDKD == 'lightcyan',]
LC_genes <- LC_genes$genes_NDKD

LC_enrich_GO <- enrichr(gene = LC_genes, 
                     databases = c('GO_Biological_Process_2025', 'GO_Cellular_Component_2025','GO_Molecular_Function_2025' ))

LC_enrich_KRW <- enrichr(gene = LC_genes, 
                   databases = c("KEGG_2021_Human", 
                                 "Reactome_2022",
                                 "WikiPathway_2023_Human"))

#no significant enrich results - module is likely artifact
