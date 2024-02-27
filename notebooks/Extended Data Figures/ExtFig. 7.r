rm(list=ls())
gc()
library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)
library(Seurat)
#Read the object in Figure 6
PFC<-readRDS("PFC.cellchat.rds")
PFC <- updateCellChat(PFC)
Hippo<-readRDS("Hippo.cellchat.rds")
Hippo <- updateCellChat(Hippo)
levels(Hippo@idents)
My_levels <- c("Oligodendrocyte precursor", "Oligodendrocyte", "Newly formed oligodendrocytes", "Astrocyte", "Microglia", "Microglia(PCDH9high)", "Vascular leptomeningeal cell", "Fibroblast", "Endothelial", "GABAergic", "Glutamatergic", "Glutamatergic(Pyr)")
Hippo@idents <- factor(Hippo@idents, levels= My_levels)
levels(Hippo@idents)
levels(PFC@idents)
My_levels <- c("Oligodendrocyte precursor", "Oligodendrocyte", "Newly formed oligodendrocytes", "Astrocyte", "Microglia", "Microglia(PCDH9high)", "Vascular leptomeningeal cell", "Fibroblast", "Endothelial", "GABAergic", "Glutamatergic", "Glutamatergic(Pyr)")
PFC@idents <- factor(PFC@idents, levels= My_levels)
levels(PFC@idents)
PFC <- updateCellChat(PFC)
Hippo <- updateCellChat(Hippo)
PFC <- updateCellChat(PFC)
Hippo <- updateCellChat(Hippo)
group.new = levels(Hippo@idents)
PFC <- liftCellChat(PFC, group.new)
object.list <- list(PFC = PFC, Hippo = Hippo)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

##extended_data_fig7a
#Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 16,color.use=c("#7DF9FF", "#006091", "#5C89CC", "#74caff", "#AFE1AF", "#F0E68C", "#8B8000", "#D2B48C", "#50C878", "#E30B5C", "#CC5500", "#FA5F55"))
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 16,color.use=c("#7DF9FF", "#006091", "#5C89CC", "#74caff", "#AFE1AF", "#F0E68C", "#8B8000", "#D2B48C", "#50C878", "#E30B5C", "#CC5500", "#FA5F55"))
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

##extended_data_fig7b
# incoming signaling
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 16, color.heatmap = "GnBu",color.use=c("#7DF9FF", "#006091", "#5C89CC", "#74caff", "#AFE1AF", "#F0E68C", "#8B8000", "#D2B48C", "#50C878", "#E30B5C", "#CC5500", "#FA5F55"))
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 16, color.heatmap = "GnBu",color.use=c("#7DF9FF", "#006091", "#5C89CC", "#74caff", "#AFE1AF", "#F0E68C", "#8B8000", "#D2B48C", "#50C878", "#E30B5C", "#CC5500", "#FA5F55"))
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

##extended_data_fig7c
pathways.show<-c("NRG","CADM", "NEGR", "LAMININ")
par(mfrow=c(1,1))
netVisual_aggregate(Hippo, signaling = pathways.show, layout = "circle")

##extended_data_fig7d
library(Seurat)
library(ggplot2)
fig4a<-readRDS("fig4a.rds")
sce.all.list <- SplitObject(fig4a, split.by = "region")
seu<-sce.all.list[["Hippocampus"]]
##Left
L_genes<-c("CADM1", "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMB1", "LAMC1", "LAMC3", "NEGR1", "NRG1", "NRG2", "NRG3")
p1 <- DotPlot(object = seu, features = L_genes, assay = "RNA",cols = c("white", "red")) +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
    coord_flip() +
    theme(legend.position = "left",
        plot.margin = margin(r = 0),
        legend.title  = element_text (size = 10),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid.major = element_line(size = 0.2),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.1)) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right")
p1
ggsave(filename="dot-left.pdf",plot=p1, width = 7, height = 9)

#Right
R_genes<-c("CADM1", "ITGAV", "ITGB8", "ITGA1", "ITGA6", "ITGA7", "ITGA9", "ITGB1", "ITGB4", "CD44", "DAG1", "SV2A", "SV2B", "NEGR1", "ERBB3", "ERBB4")
p3 <- DotPlot(object = seu, features = R_genes, assay = "RNA",cols = c("white", "red")) +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
    coord_flip() +
    theme(axis.title = element_blank(),
        plot.margin = margin(l = 0),
        panel.border = element_rect(color = "black", fill = NA,size = 1.5),
        panel.grid.major = element_line(size = 0.2),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.1)) +
    scale_y_discrete(position = "right")  + NoLegend()
p3
ggsave(filename="dot-right.pdf",plot=p3, width = 5, height = 9)

library(patchwork)
p1  + p3
p2<-p1  + p3
ggsave(filename="extended_data_fig7d.pdf",plot=p2, width = 8, height = 6)