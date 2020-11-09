####need to do this so that negative logFC corresponds to downregulation in mutants. These are the differential result tables containing logFC, pval etc
res_B_KS1$logFC <- -res_B_KS1$log2FoldChange
res_B_KS2$logFC <- -res_B_KS2$log2FoldChange
res_B_RT$logFC <- -res_B_RT$log2FoldChange


qobj_KS2_KS1_B <- qvalue(p = res_B_KS2$P.Value[which(rownames(res_B_KS2) %in% rownames(res_B_KS1)[which(res_B_KS1$padj < 0.1)])], 
                         fdr.level = 0.1, pi0.method = "bootstrap")

KS2_KS1_B <- res_B_KS2[which(rownames(res_B_KS2) %in% rownames(res_B_KS1)[which(res_B_KS1$padj < 0.1)]), ][
  which(qobj_KS2_KS1_B$significant == TRUE), ]
KS2_KS1_B$logFC <- sapply(rownames(KS2_KS1_B), function(xx) 
  median(res_B_KS1$logFC[which(rownames(res_B_KS1) == xx)], 
         res_B_KS2$logFC[which(rownames(res_B_KS2) == xx)]))

#directionality
common_genes_KS1_KS2_only <- data.frame(gene_id = rownames(KS2_KS1_B), 
                                        gene_name = as.character(KS2_KS1_B$gene_symbol), 
                                        logFC_KS1 = sapply(rownames(KS2_KS1_B), 
                                                           function(xx) -res_B_KS1$log2FoldChange[which(rownames(res_B_KS1) == xx)]), 
                                        logFC_KS2 = sapply(rownames(KS2_KS1_B), 
                                                           function(xx) -res_B_KS2$log2FoldChange[which(rownames(res_B_KS2) == xx)]), 
                                        stringsAsFactors = FALSE)

common_genes_KS1_KS2_only <- common_genes_KS1_KS2_only[
  which(common_genes_KS1_KS2_only$logFC_KS1*common_genes_KS1_KS2_only$logFC_KS2 > 0), ]
KS2_KS1_B <- KS2_KS1_B[which(rownames(KS2_KS1_B) %in% common_genes_KS1_KS2_only$gene_id), ]

qobj_KS2_KS1_RT <- qvalue(p = res_B_RT$P.Value[which(rownames(res_B_RT) %in% rownames(KS2_KS1_B))], fdr.level = 0.1, pi0.method = "bootstrap")
KS2_KS1_RT_B <- res_B_RT[which(rownames(res_B_RT) %in% rownames(KS2_KS1_B)), ][
  which(qobj_KS2_KS1_RT$significant == TRUE), ]

common_genes_B <- data.frame(gene_id = rownames(KS2_KS1_RT_B), 
                             gene_name = as.character(KS2_KS1_RT_B$gene_symbol), 
                             logFC_KS1 = sapply(rownames(KS2_KS1_RT_B), 
                                                function(xx) res_B_KS1$logFC[which(rownames(res_B_KS1) == xx)]), 
                             logFC_KS2 = sapply(rownames(KS2_KS1_RT_B), 
                                                function(xx) res_B_KS2$logFC[which(rownames(res_B_KS2) == xx)]), 
                             logFC_RT = sapply(rownames(KS2_KS1_RT_B), 
                                               function(xx) res_B_RT$logFC[which(rownames(res_B_RT) == xx)]), 
                             stringsAsFactors = FALSE)

gene_ids_to_use1 <- common_genes_B$gene_id[which(common_genes_B$logFC_KS1 < 0 & 
                                                   common_genes_B$logFC_KS2 < 0 & 
                                                   common_genes_B$logFC_RT < 0)]
gene_ids_to_use2 <- common_genes_B$gene_id[which(common_genes_B$logFC_KS1 > 0 & 
                                                   common_genes_B$logFC_KS2 > 0 & 
                                                   common_genes_B$logFC_RT > 0)]

common_genes_B <- common_genes_B[which(common_genes_B$gene_id %in% c(gene_ids_to_use1, gene_ids_to_use2)), ]
write_csv(common_genes_KS1_KS2_only, "common_genes_B_KS1_KS2_only.csv")
write_csv(common_genes_B, "common_genes_B_KS1_KS2_and_RT.csv")

####Figures
quartz(file = "KS1_vs_KS2_rnaseq_B_cells.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2), mar = c(5,4,1,1) + 0.1)
compareDEAnalyses(res_B_KS2, res_B_KS1, max(res_B_KS1$P.Value[which(res_B_KS1$padj < 0.1)]), 
                  "p-value_density", "RNA-seq", "KS2 p-values", "Density", 
                  "KS1 FDR < 0.1", "KS1 FDR >= 0.1", y_lim = c(0, 3.5))

compareDEAnalyses(res_B_KS2, res_B_KS1, max(res_B_KS1$P.Value[which(res_B_KS1$padj < 0.1)]), "log2FC", 
                  "RNA-seq", "KS1 logFC (KS1 FDR < 0.1)", "KS2 logFC (KS1 FDR < 0.1)", 
                  "", "", x_lim = c(-3, 3), y_lim = c(-2,2))

dev.off()

quartz(file = "RT_vs_KS1_and_2_common_hits_rnaseq_B_cells.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
compareDEAnalyses(res_B_RT, KS2_KS1_B, max(KS2_KS1_B$P.Value), 
                  "p-value_density", "RNA-seq", "RT1 p-values", "Density", 
                  "KS1/2 common hits", "other", y_lim = c(0, 7), bwidth = 0.035)

compareDEAnalyses(res_B_RT, KS2_KS1_B, max(KS2_KS1_B$P.Value), "log2FC", 
                  "RNA-seq", "common KS1/2 hits\nmedian(logFC)", "RT1 logFC", 
                  "", "", x_lim = c(-3, 3), y_lim = c(-2,2))

dev.off()



##########PCA with common genes

quartz(file = "pca_from_common_genes_RT.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
vstab <- vst(dds_B_RT)

pcaData <- plotPCA(vstab[which(rownames(vstab) %in% c(gene_ids_to_use1, gene_ids_to_use2)), ], 
                   returnData=TRUE, ntop = length(c(gene_ids_to_use1, gene_ids_to_use2)))
colnames(pcaData)[4] <- "genotype"
percentVar <- round(100 * attr(pcaData, "percentVar"))

#ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
#  geom_point(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#  coord_fixed() + scale_fill_manual(values=c("WT"="blue", "R172K"="red")) + labs(fill="Type", title = "RNA-seq") + 
#  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme_classic()
plot(pcaData$PC1[which(pcaData$genotype == "RT1")], pcaData$PC2[which(pcaData$genotype == "RT1")], 
     cex = 1.15, col = "deep pink", xlab = paste0("PC1 (", percentVar[1], "%)"), ylab = paste0("PC2 (", percentVar[2], "%)"), 
     main = "", pch = 19, xlim = c(-3.3, 2), ylim = c(-1, 1.5), bty = 'l', xaxt = 'n', yaxt = 'n', cex.lab = 1.25)
points(pcaData$PC1[which(pcaData$genotype == "WT")], pcaData$PC2[which(pcaData$genotype == "WT")], 
     cex = 1.15, col = rgb(0, 0, 0, 0.8), pch = 19)
legend <- legend("top", legend = c("RT1", "WT"), bty = 'n', cex = 0.82, col = c("deep pink", rgb(0, 0, 0, 0.8)), pch = 19)
axis(1, at = c(-3,2), cex.axis = 1.2)
axis(2, at = c(-1, 1), cex.axis = 1.2)
dev.off()

quartz(file = "pca_from_common_genes_KS1_KS2.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
vstab <- vst(dds_B_KS1_KS2)

pcaData <- plotPCA(vstab[which(rownames(vstab) %in% c(gene_ids_to_use1, gene_ids_to_use2)), ], 
                   returnData=TRUE, ntop = length(c(gene_ids_to_use1, gene_ids_to_use2)))
colnames(pcaData)[4] <- "genotype"
percentVar <- round(100 * attr(pcaData, "percentVar"))

#ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
#  geom_point(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#  coord_fixed() + scale_fill_manual(values=c("WT"="blue", "R172K"="red")) + labs(fill="Type", title = "RNA-seq") + 
#  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme_classic()

plot(pcaData$PC1[which(pcaData$genotype == "KS1")], pcaData$PC2[which(pcaData$genotype == "KS1")], 
     cex = 1.15, col = "orange", xlab = paste0("PC1 (", percentVar[1], "%)"), ylab = paste0("PC2 (", percentVar[2], "%)"), 
     main = "", pch = 19, xlim = c(-3.3, 3.3), ylim = c(-2, 3), bty = 'l', xaxt = 'n', yaxt = 'n', cex.lab = 1.25)
points(pcaData$PC1[which(pcaData$genotype == "WT")], pcaData$PC2[which(pcaData$genotype == "WT")], 
       cex = 1.15, col = rgb(0, 0, 0, 0.8), pch = 19)
points(pcaData$PC1[which(pcaData$genotype == "KS2")], pcaData$PC2[which(pcaData$genotype == "KS2")], 
       cex = 1.15, col = "forest green", pch = 19)
legend <- legend("top", legend = c("KS1", "KS2", "WT"), bty = 'n', cex = 0.82, col = c("orange", "forest green", rgb(0, 0, 0, 0.8)), pch = 19)
axis(1, at = c(-3,3), cex.axis = 1.2)
axis(2, at = c(-2, 2), cex.axis = 1.2)
dev.off()

quartz(file = "pca_from_common_genes_KS2.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
vstab <- vst(dds_B_KS2)

pcaData <- plotPCA(vstab[which(rownames(vstab) %in% c(gene_ids_to_use1, gene_ids_to_use2)), ], 
                   returnData=TRUE, ntop = length(c(gene_ids_to_use1, gene_ids_to_use2)))
colnames(pcaData)[4] <- "genotype"
percentVar <- round(100 * attr(pcaData, "percentVar"))

#ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
#  geom_point(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#  coord_fixed() + scale_fill_manual(values=c("WT"="blue", "R172K"="red")) + labs(fill="Type", title = "RNA-seq") + 
#  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme_classic()

plot(pcaData$PC1[which(pcaData$genotype == "KS2")], pcaData$PC2[which(pcaData$genotype == "KS2")], 
     cex = 1.15, col = "forest green", xlab = paste0("PC1 (", percentVar[1], "%)"), ylab = paste0("PC2 (", percentVar[2], "%)"), 
     main = "", pch = 19, xlim = c(-2.2, 2.2), ylim = c(-1.8, 1.8), bty = 'l', xaxt = 'n', yaxt = 'n', cex.lab = 1.25)
points(pcaData$PC1[which(pcaData$genotype == "WT")], pcaData$PC2[which(pcaData$genotype == "WT")], 
       cex = 1.15, col = rgb(0, 0, 0, 0.8), pch = 19)
legend <- legend("top", legend = c("KS2", "WT"), bty = 'n', cex = 0.82, col = c("forest green", rgb(0, 0, 0, 0.8)), pch = 19)
axis(1, at = c(-2,2), cex.axis = 1.2)
axis(2, at = c(-1.5, 1.5), cex.axis = 1.2)
dev.off()


##############
quartz(file = "RNA_effect_direction.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(5.25, 5.25, 1.25, 1.25) + 0.1)
plot(1, prop.table(table(as.factor(res_B_KS1$logFC[which(res_B_KS1$padj < 0.1)] > 0)))[2], pch = 19, 
     col = alpha("red", 0.75), xlab = "", ylab = "% changes towards\ngreater expression", main = "", font.main = 1, 
     xlim = c(0.8, 5.2), ylim = c(0, 1), xaxt = 'n', bty = 'l', yaxt = 'n', bty = "l", cex = 1.25, cex.lab = 1.2)
points(2, prop.table(table(as.factor(res_B_KS2$logFC[which(res_B_KS2$padj < 0.1)] > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(3, prop.table(table(as.factor(res_B_RT$logFC[which(res_B_RT$padj < 0.1)] > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(4, prop.table(table(as.factor(KS2_KS1_B$logFC > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(5, prop.table(table(as.factor(KS2_KS1_RT_B$logFC > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
#points(7, 1-qobj_RT_KS2_KS1_B_atac_promoters$pi0, pch = 19, 
#       col = alpha("brown", 0.75), xlab = "", ylab = "% shared differentially accessible\nregulatory elements", main = "")

#points(8, 1-qobj_RT_KS2_KS1_B_atac_non_promoters$pi0, pch = 19, 
#       col = alpha("forest green", 1), xlab = "", ylab = "% shared differentially accessible\nregulatory elements", main = "")

axis(1, at = c(1, 2, 3, 4, 5), labels = c("KS1", "KS2", "RT", "KS1/KS2\ncommon", "KS1/KS2/RT\ncommon"), cex.axis = 0.84, las = 2)
axis(2, at = c(0, 0.5, 1), labels = c(0, 50, 100))
abline(v = c(1.5, 2.5, 3.5, 4.5, 5.5), lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()



###pathway analysis
expression_pathways <- getReactomeEnrichedPathways(unique(common_genes_B$gene_id), 
                                                   unique(unlist(Reduce(intersect, list(rownames(res_B_KS1), 
                                                                                        rownames(res_B_KS2), rownames(res_B_RT))))))







