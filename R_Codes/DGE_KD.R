##### Differential Gene Expression ####

#load libraries 
library(AnnotationDbi) 
library(org.Hs.eg.db) 
library(ggplot2)
library(ggfortify)
library(ggVennDiagram)
library(writexl)
library(clusterProfiler)
library(limma)
library(dplyr)

NDKD = readRDS(file = "NDKD_processed.RDS")
DKD = readRDS(file = "DKD_processed.RDS")

#### DKD data set ####

deg_DKD_all <- topTable(DKD$fit, coef=2, sort.by="p", number = Inf)
deg_DKD_all$DF <- "NS" #make all not significant to start
deg_DKD_all$DF[deg_DKD_all$logFC > 1.5 & deg_DKD_all$adj.P.Val < 0.05] <- "UP"
deg_DKD_all$DF[deg_DKD_all$logFC < -1.5 & deg_DKD_all$adj.P.Val < 0.05] <- "DOWN"

#add in gene symbols
deg_DKD_all$symbol <- mapIds(x = org.Hs.eg.db,
                             keys = rownames(deg_DKD_all),
                             column = 'SYMBOL',
                             keytype = 'ENTREZID',
                             multiVals = 'first')
deg_DKD_all <- deg_DKD_all[!is.na(deg_DKD_all$symbol),] #remove genes unable to match to a symbol
table(deg_DKD_all$DF)

ggplot(data = deg_DKD_all, aes(x = logFC, y = -log10(P.Value), col = DF)) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  theme_minimal()+
  labs(x = 'Log2 Fold Change')+
  xlim(-5,5)+
  scale_color_manual(values = c("#00AFBB", "grey", "#FF4500"),
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  ggtitle("Progressive Diabetic Kidney Disease DEGs")+
  theme(plot.title = element_text(hjust = 0.5))+ #title centered
  theme(
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)
  )
ggsave('DKD_Volcano.pdf', width = 8, height = 5) 

#### NDKD data set ####

deg_NDKD_all <- topTable(NDKD$fit, coef=2, sort.by="p", number = Inf)
deg_NDKD_all$DF <- "NS" 
deg_NDKD_all$DF[deg_NDKD_all$logFC > 1.5 & deg_NDKD_all$adj.P.Val < 0.05] <- "UP"
deg_NDKD_all$DF[deg_NDKD_all$logFC < -1.5 & deg_NDKD_all$adj.P.Val < 0.05] <- "DOWN"


deg_NDKD_all$symbol <- mapIds(x = org.Hs.eg.db,
                              keys = rownames(deg_NDKD_all),
                              column = 'SYMBOL',
                              keytype = 'ENTREZID',
                              multiVals = 'first')
deg_NDKD_all <- deg_NDKD_all[!is.na(deg_NDKD_all$symbol),] 
table(deg_NDKD_all$DF)

ggplot(data = deg_NDKD_all, aes(x = logFC, y = -log10(P.Value), col = DF)) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  theme_minimal()+
  labs(x = 'Log2 Fold Change')+
  xlim(-5,5)+
  ylim(0,8)+
  scale_color_manual(values = c("#00AFBB", "grey", "#FF4500"),
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  ggtitle("Non-Diabetic Progressive Kidney Disease DEGs")+
  theme(plot.title = element_text(hjust = 0.5)) 
ggsave('NDKD_Volcano.pdf', width = 8, height = 5) 


####Comparing overlap of up-regulated and down-regulated genes
deg_dkd_sig <- deg_DKD_all[deg_DKD_all$DF != "NS",]
deg_ndkd_sig <- deg_NDKD_all[deg_NDKD_all$DF != "NS",]
DEG_comparision_UP <- list(DKD = deg_dkd_sig$symbol[deg_dkd_sig$DF == 'UP'],
                           NDKD = deg_ndkd_sig$symbol[deg_ndkd_sig$DF == 'UP'])

ggVennDiagram(DEG_comparision_UP) + coord_flip() +
  scale_fill_distiller(palette = "Reds", direction = 1,limits = c(200, 700))+
  ggtitle('Up-regulated Genes') +
  theme(plot.title = element_text(hjust = 0.5))
ggsave('Upreg_ven.pdf', width = 8, height = 6)

DEG_comparision_DOWN <- list(DKD = deg_dkd_sig$symbol[deg_dkd_sig$DF == 'DOWN'],
                             NDKD = deg_ndkd_sig$symbol[deg_ndkd_sig$DF == 'DOWN'])

ggVennDiagram(DEG_comparision_DOWN) + coord_flip() +
  scale_fill_distiller(palette = "Blues", direction = 1,limits = c(0, 400))+
  ggtitle('Down-regulated Genes') +
  theme(plot.title = element_text(hjust = 0.5))
ggsave('Downreg_ven.pdf', width = 8, height = 6)

write_xlsx(deg_NDKD_all, 'NDKD_DEGS.xlsx')
write_xlsx(deg_DKD_all, 'DKD_DEGS.xlsx')

#### Gene Ontology for common genes ####

#Up regulated genes
common_UP_DEGS <- intersect(DEG_comparision_UP$DKD, DEG_comparision_UP$NDKD)
common_UP_res <- enrichGO(gene =  common_UP_DEGS, OrgDb = "org.Hs.eg.db", 
                          keyType = "SYMBOL",
                          ont = "BP")
common_UP_res <- clusterProfiler::simplify(common_UP_res)
dotplot(common_UP_res,showCategory=20,font.size=10,label_format=70)+
  theme_minimal() +
  ggtitle("GO Enrichment of common up-regulated genes")
ggsave('GO_BP_common_up.pdf', width = 8, height = 8)
common_up_DF <- common_UP_res@result
write_xlsx(common_up_DF, 'GO_BP_common_up.xlsx')

#Down regulated genes
common_DOWN_DEGS <- intersect(DEG_comparision_DOWN$DKD, DEG_comparision_DOWN$NDKD)
common_DOWN_res <- enrichGO(gene =  common_DOWN_DEGS, OrgDb = "org.Hs.eg.db", 
                                keyType = "SYMBOL",
                                ont = "BP")
#Note this does not produce any results - not enough genes
#dotplot(common_DOWN_res,showCategory=20,font.size=10,label_format=70)+
#  theme_minimal() +
#  ggtitle("GO Enrichment of common down-regulated genes")
#ggsave('GO_BP_common_down.pdf', width = 8, height = 8)
#common_down_DF <- common_DOWN_res@result

#### Gene Ontology for DKD specific genes ####
DKD_UP_DEGS <- DEG_comparision_UP$DKD[!DEG_comparision_UP$DKD %in% common_UP_DEGS]
DKD_UP_res <- enrichGO(gene =  DKD_UP_DEGS, OrgDb = "org.Hs.eg.db", 
                       keyType = "SYMBOL",
                       ont = "BP")
DKD_UP_res <- clusterProfiler::simplify(DKD_UP_res)
dotplot(DKD_UP_res,showCategory=20,font.size=10,label_format=70)+
  theme_minimal() +
  ggtitle("GO Enrichment of DKD up-regulated genes")
ggsave('GO_BP_DKD_up.pdf', width = 8, height = 8)
DKD_up_df <- DKD_UP_res@result
write_xlsx(DKD_up_df, 'GO_BP_DKD_up.xlsx')

#create and save an additional plot that shows terms that are only present in DKD
#retrieve GO ID's
DKD_spec_GO <- DKD_up_df$ID
common_spec_GO <- common_up_DF$ID
#remove the common ID's to isolate DKD specific terms
DKD_only <- DKD_spec_GO[!DKD_spec_GO %in% common_spec_GO]
DKD_only_DF <- DKD_UP_res
DKD_only_DF@result <- DKD_UP_res@result %>%
  filter(ID %in% DKD_only)
dotplot(DKD_only_DF,showCategory=20,font.size=10,label_format=70)+
  theme_minimal() +
  ggtitle("GO BP Enrichment of DKD specific up-regulated genes")
ggsave('DKD_specific_GO_BP.pdf', width = 8, height = 8)

DKD_DOWN_DEGS <- DEG_comparision_DOWN$DKD[!DEG_comparision_DOWN$DKD %in% common_DOWN_DEGS]
DKD_DOWN_res <- enrichGO(gene =  DKD_DOWN_DEGS, OrgDb = "org.Hs.eg.db", 
                             keyType = "SYMBOL",
                             ont = "BP")
DKD_DOWN_res <- clusterProfiler::simplify(DKD_DOWN_res)
dotplot(DKD_DOWN_res,showCategory=20,font.size=10,label_format=70)+
  theme_minimal() +
  ggtitle("GO BP Enrichment of DKD down-regulated genes")
ggsave('GO_BP_DKD_DOWN.pdf', width = 8, height = 8)

DKD_down_df <- DKD_DOWN_res@result
write_xlsx(DKD_down_df, 'GO_BP_DKD_down.xlsx')


#### Gene Ontology for NDKD specific genes
NDKD_UP_DEGS <- DEG_comparision_UP$NDKD[!DEG_comparision_UP$NDKD %in% common_UP_DEGS]

NDKD_UP_res <- enrichGO(gene =  NDKD_UP_DEGS, OrgDb = "org.Hs.eg.db", 
                           keyType = "SYMBOL",
                           ont = "BP")

NDKD_UP_res <- clusterProfiler::simplify(NDKD_UP_res)
dotplot(NDKD_UP_res,showCategory=20,font.size=10,label_format=70)+
  theme_minimal() +
  ggtitle("GO BP Enrichment of NDKD up-regulated genes")
ggsave('GO_BP_NDKD_up.pdf', width = 8, height = 8)
NDKD_up_df <- NDKD_UP_res@result
write_xlsx(NDKD_up_df, 'GO_BP_NDKD_up.xlsx')

#create and save an additional plot that shows terms that are only present in NDKD
NDKD_spec_GO <- NDKD_up_df$ID
NDKD_only <- NDKD_spec_GO[!NDKD_spec_GO %in% common_spec_GO]
NDKD_only_DF <- NDKD_UP_res
NDKD_only_DF@result <- NDKD_UP_res@result %>%
  filter(ID %in% NDKD_only)
dotplot(NDKD_only_DF,showCategory=20,font.size=10,label_format=70)+
  theme_minimal() +
  ggtitle("GO BP Enrichment of NDKD specific up-regulated genes")
ggsave('NDKD_BP_specific_up.pdf', width = 8, height = 8)

NDKD_DOWN_DEGS <- DEG_comparision_DOWN$NDKD[!DEG_comparision_DOWN$NDKD %in% common_DOWN_DEGS]
NDKD_DOWN_res <- enrichGO(gene =  NDKD_DOWN_DEGS, OrgDb = "org.Hs.eg.db", 
                              keyType = "SYMBOL",
                              ont = "BP")
NDKD_DOWN_res <- clusterProfiler::simplify(NDKD_DOWN_res)
dotplot(NDKD_DOWN_res,showCategory=20,font.size=10,label_format=70)+
  theme_minimal() +
  ggtitle("GO BP Enrichment of NDKD down-regulated genes")
ggsave('GO_BP_NDKD_DOWN.pdf', width = 8, height = 8)

