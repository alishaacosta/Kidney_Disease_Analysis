##### Preprocessing ####
##Downloading raw counts and metadata
##Limma-voom workflow


#load libraries 

library(GEOquery)
library(edgeR)
library(limma)
library(ggplot2)
library(ggfortify)

###########Retrieving raw count matrices and metadata###########################

#### Download non-diabetic data set, retrieve data from cohort 2 only
GSE137570 <- getGEO(GEO = 'GSE137570')
NDKD <- GSE137570$`GSE137570-GPL11154_series_matrix.txt.gz`
NDKD_meta <- pData(NDKD)
#retrieve necessary metadata only
NDKD_meta <- NDKD_meta[,-c(3:41)] 
NDKD_meta <- NDKD_meta[,-4]
#indicated progression level, 0 = not progressive, 1 = progressive
colnames(NDKD_meta)[4] <- 'group' 
NDKD_meta$state <- ifelse(NDKD_meta$group == 0, 'non-progressive', 'progressive')


raw_NDKD <- gzfile('GSE137570_raw_counts_GRCh38.p13_NCBI.tsv.gz')
raw_NDKD_counts <- read.table(raw_NDKD, header = TRUE, sep = '\t') 
rownames(raw_NDKD_counts) <- raw_NDKD_counts$GeneID
raw_NDKD_counts <- raw_NDKD_counts[,-1]
samples_NDKD <- rownames(NDKD_meta)
raw_NDKD_counts <- raw_NDKD_counts[, samples_NDKD]
#ensuring column names of the raw counts matrices are in the same order as the metadata
table(colnames(raw_NDKD_counts) == rownames(NDKD_meta))

#### Download diabetic data set
GSE142025 <- getGEO(GEO = 'GSE142025')
DKD <- GSE142025$GSE142025_series_matrix.txt.gz
DKD_meta <- pData(DKD)
DKD_meta <- DKD_meta[,-c(2:36)]
colnames(DKD_meta)[2] <- "group"
DKD_meta$group <- as.numeric(factor(DKD_meta$group,
                                         levels = c("Control",
                                                    "Early_DN",
                                                    "Advanced_DN"),
                                         ordered = TRUE))


raw_DKD <- gzfile('GSE142025_raw_counts_GRCh38.p13_NCBI.tsv.gz')
DKD_raw_counts <- read.table(raw_DKD, header = TRUE, sep = '\t')
rownames(DKD_raw_counts) <- DKD_raw_counts$GeneID
DKD_raw_counts <- DKD_raw_counts[,-1]

DKD_meta <- DKD_meta[DKD_meta$group != 1,] #remove non disease samples
DKD_meta$state <- ifelse(DKD_meta$group == 3, 'Advanced', 'Early')
DKD_raw_counts <- DKD_raw_counts[,rownames(DKD_meta)]

#ensuring column names of the raw counts matricies are in the same order as the metadata
table(colnames(DKD_raw_counts) == rownames(DKD_meta))

###########################Limma-voom workflow###########################

dge_DKD <- DGEList(counts = DKD_raw_counts,
                   samples = DKD_meta,
                   group = DKD_meta$group)
keep_DKD <- filterByExpr(dge_DKD)
dge_DKD <- dge_DKD[keep_DKD,,keep.lib.sizes = FALSE]
rm(keep_DKD)
dim(dge_DKD$counts) #check number of genes and samples
dge_DKD <- calcNormFactors(dge_DKD, method = "TMM") 

pdf(file = "DKD_voom.pdf", width = 8, height = 5)
v_DKD <- voom(dge_DKD, plot = TRUE)
dev.off()


fit_DKD <- lmFit(v_DKD) #calculate the log fold change and standard errors
fit_DKD <- eBayes(fit_DKD) #empirical Bayes smoothing to the standard errors

pca_DKD <- prcomp(t(v_DKD$E))
sum_DKD <- summary(pca_DKD)
variance_DKD <- sum_DKD$importance
variance_DKD[,1:2]

pdf(file = "PCA_DKD.pdf", width = 8, height = 5)
autoplot(pca_DKD, data = DKD_meta, color = 'state') +
  scale_color_manual(values = c('orange', 'purple'))+
  theme_minimal()+
  ggtitle("PCA DKD")
dev.off()

####Repeat for NDKD
dge_NDKD <- DGEList(counts = raw_NDKD_counts,
                   samples = NDKD_meta,
                   group = NDKD_meta$group)
keep_NDKD <- filterByExpr(dge_NDKD)
dge_NDKD <- dge_NDKD[keep_NDKD,,keep.lib.sizes = FALSE]
rm(keep_NDKD)
dim(dge_NDKD$counts)
dge_NDKD <- calcNormFactors(dge_NDKD, method = "TMM") 

pdf(file = 'NDKD_voom.pdf', width = 8, height = 5)
v_NDKD <- voom(dge_NDKD, plot = TRUE)
dev.off()

fit_NDKD <- lmFit(v_NDKD) 
fit_NDKD <- eBayes(fit_NDKD) 
pca_NDKD <- prcomp(t(v_NDKD$E))

sum_NDKD <- summary(pca_NDKD)
variance_NDKD <- sum_NDKD$importance
variance_NDKD[,1:2]

pdf(file = 'PCA_NDKD.pdf', width = 8, height = 5)
autoplot(pca_NDKD, data = NDKD_meta, color = 'state') +
  scale_color_manual(values = c('blue', 'green'))+
  theme_minimal()+
  ggtitle("PCA NDKD")
dev.off()

#save processed files
saveRDS(list(v = v_NDKD, fit = fit_NDKD), "NDKD_processed.RDS")
saveRDS(list(v = v_DKD, fit = fit_DKD), "DKD_processed.RDS")
