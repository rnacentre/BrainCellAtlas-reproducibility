#edgeR
library(Libra)
library(MAST)
library(Seurat)
library(dplyr)
library(EnhancedVolcano)
library(stringr)
library(tidyverse)
library(edgeR)
library(readxl)
library(ggvenn)

edgeR_DE <- function(counts, group, design, cont, n_top_genes = Inf, method = 'TMM') {
  y <- DGEList(counts = counts, group = group) %>%
    calcNormFactors(method = method) %>%
    estimateDisp(design)
  
  fit <- glmFit(y, design = design)
  if (is.null(cont)) {
      test <- glmLRT(fit)
  } else {
        contrast <- makeContrasts(contrasts = cont, levels = design)
        test <- glmLRT(fit, contrast = contrast)
  }
  res <- topTags(test, n = n_top_genes) %>%
    as.data.frame() %>%
    rownames_to_column('gene')
  
  return(res)
}

npc_meta = readRDS('npc/data/annot.rds')
i = 'NPC'
npc$sample_group = ifelse(npc$state==i, i, "Others")
npc$sample_group = factor(npc$sample_group, levels = c('Others', i))
mat = to_pseudobulk(npc, cell_type_col = "group", label_col = "sample_group", replicate_col = 'sampleID')

group = npc@meta.data %>% distinct(sampleID, sample_group, .keep_all = T) %>% select(sampleID, sample_group) %>% remove_rownames()
rownames(group) = paste0(group$sampleID,':',group$sample_group)
group

counts = mat[[1]][rownames(group)]
counts                 
                 
group$sampleID = as.character(group$sampleID)
group$sampleID %>% unique
                 
# for TMM
design <- model.matrix(~
                        sampleID+
                       sample_group, data = group)

resTmm = edgeR_DE(counts, group = group$group, design, method = 'TMM', cont = 'sample_groupNPC')
write.csv(resTmm, 'npc/results/npc_rmdonor_DEG.csv')

                 





                 
                 


library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
library(purrr)
library(Seurat)

library(sceasy)
ad = convertFormat("~/projA/data/hipngs/hm3f4a_noBr_hip_integration_rmPyr_Micro_qc1_reharmony_anno_a2n_sp100.h5ad", from="anndata", to="seurat", outFile=NULL)
obs<-ad@meta.data



#all NPC
deg_all_npc<-FindMarkers(ad,ident.1 = "NPC", group.by = "anno",test.use = "wilcox",logfc.threshold = 0.1)
write.csv(deg_all_npc,file="~/projA/data/hipngs/hm3f4a_noBr_hip_integration_rmPyr_Micro_qc1_reharmony_anno_a2n_sp100_DEG.csv")
sum((deg_all_npc$p_val_adj<0.05) & (abs(deg_all_npc$avg_log2FC)>0.5))
write.csv(deg_all_npc[((deg_all_npc$p_val_adj<0.05) & (abs(deg_all_npc$avg_log2FC)>0.5)),],file="~/projA/data/hipngs/hm3f4a_noBr_hip_integration_rmPyr_Micro_qc1_reharmony_anno_a2n_sp100_DEG1.csv")


library(EnhancedVolcano)

EnhancedVolcano(deg_all_npc,
    lab = rownames(deg_all_npc),
    selectLab = c('UBE2C', 'HMGB2', 'TOP2A'),
    x = 'avg_log2FC',
    y = 'p_val',
    #drawConnectors = TRUE,
    #pointSize = 1,
    #labSize = 2.5,
    xlim = c(-2, 2),
    gridlines.major = FALSE,
    gridlines.minor = FALSE)



options(repr.plot.width = 4, repr.plot.height = 4)
EnhancedVolcano(deg_all_npc,
    lab = rownames(deg_all_npc),
    selectLab = c('UBE2C', 'HMGB2', 'TOP2A', 'EOMES', 'CKS2', 'SMC4', 'HES6', 'TUBAIB'),
    x = 'avg_log2FC',
    y = 'p_val',
    #drawConnectors = TRUE,
    #pointSize = 1,
    #labSize = 2.5,
    xlim = c(-2, 2),
    gridlines.major = FALSE,
    gridlines.minor = FALSE)


deg_stage_npc<-lapply(names(table(obs$stage)),function(g){
  ad1<-subset(ad,subset=stage==g)
  degdf <- FindMarkers(ad1,ident.1 = "NPC", group.by = "anno",test.use = "wilcox")
  return(degdf)
})
names(deg_stage_npc)<- names(table(obs$stage))



tmp<-deg_homo_neu[-log10(deg_homo_neu$p_val_adj)<300,]
EnhancedVolcano(tmp,
    lab = rownames(tmp),
    #selectLab = c('UBE2C', 'HMGB2', 'TOP2A'),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    labhjust = 1,
    #drawConnectors = TRUE,
    #pointSize = 1,
    #labSize = 2.5,
    FCcutoff = 0.5,
    xlim = c(-1, 1),
    gridlines.major = FALSE,
    gridlines.minor = FALSE)



display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}



display_venn(
  lapply(names(table(obs$stage)),function(i)rownames(deg_stage_npc[[i]])),
  category.names = names(table(obs$stage)),
  fill = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
  )


degwcx <- lapply(names(table(obs$select)),function(s)FindMarkers(ad,ident.1 = s, group.by = "select",test.use = "wilcox",logfc.threshold = 0.1))
names(degwcx)<-names(table(obs$select))


