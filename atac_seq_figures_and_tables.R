res_B_KS1_atac$logFC <- -res_B_KS1_atac$log2FoldChange
res_B_KS2_atac$logFC <- -res_B_KS2_atac$log2FoldChange
res_B_RT_atac$logFC <- -res_B_RT_atac$log2FoldChange

atac_B_KS1_granges$logFC <- -atac_B_KS1_granges$log2FoldChange
atac_B_KS2_granges$logFC <- -atac_B_KS2_granges$log2FoldChange
atac_B_RT_granges$logFC <- -atac_B_RT_granges$log2FoldChange

###promoters
atac_B_KS1_granges2 <- atac_B_KS1_granges[-which(is.na(atac_B_KS1_granges$prom_overlapping_id))]
atac_B_RT_granges2 <- atac_B_RT_granges[-which(is.na(atac_B_RT_granges$prom_overlapping_id))]
atac_B_KS2_granges2 <- atac_B_KS2_granges[-which(is.na(atac_B_KS2_granges$prom_overlapping_id))]


###KS1 vs KS2
qobj_KS2_KS1_B_atac_promoters <- qvalue(p = atac_B_KS2_granges2$P.Value[which(atac_B_KS2_granges2$GeneID %in% 
                                                                      atac_B_KS1_granges2$GeneID[which(atac_B_KS1_granges2$padj < 0.1)])], 
                              fdr.level = 0.1, pi0.method = "bootstrap")
KS2_KS1_B_atac_promoters <- atac_B_KS2_granges2[which(atac_B_KS2_granges2$GeneID %in% 
                                                        atac_B_KS1_granges2$GeneID[which(atac_B_KS1_granges2$padj < 0.1)])][
                                                          which(qobj_KS2_KS1_B_atac_promoters$significant == TRUE)]

common_promoters_KS1_KS2_only <- annoGR2DF(KS2_KS1_B_atac_promoters)
common_promoters_KS1_KS2_only$logFC_KS1 <- sapply(common_promoters_KS1_KS2_only$GeneID, function(xx) 
  res_B_KS1_atac$logFC[which(rownames(res_B_KS1_atac) == xx)])
common_promoters_KS1_KS2_only$logFC_KS2 <- sapply(common_promoters_KS1_KS2_only$GeneID, function(xx) 
  res_B_KS2_atac$logFC[which(rownames(res_B_KS2_atac) == xx)])
common_promoters_KS1_KS2_only <- common_promoters_KS1_KS2_only[
  which(common_promoters_KS1_KS2_only$logFC_KS1*common_promoters_KS1_KS2_only$logFC_KS2 > 0), ]


###RT vs KS1/2
KS2_KS1_B_atac_promoters <- KS2_KS1_B_atac_promoters[which(KS2_KS1_B_atac_promoters$GeneID %in% common_promoters_KS1_KS2_only$GeneID)]
tab1_2_B_common <- res_B_KS1_atac[which(rownames(res_B_KS1_atac) %in% KS2_KS1_B_atac_promoters$GeneID), ]
tab1_2_B_common$logFC <- sapply(rownames(tab1_2_B_common), 
                                function(xx) mean(c(res_B_KS1_atac$logFC[which(rownames(res_B_KS1_atac) == xx)], 
                                                    res_B_KS2_atac$logFC[which(rownames(res_B_KS2_atac) == xx)])))


qobj_RT_KS2_KS1_B_atac_promoters <- qvalue(p = atac_B_RT_granges2$P.Value[which(atac_B_RT_granges2$GeneID %in% 
                                                                                  KS2_KS1_B_atac_promoters$GeneID)], 
                                           fdr.level = 0.1, pi0.method = "bootstrap")
RT_KS2_KS1_B_atac_promoters <- atac_B_RT_granges2[which(atac_B_RT_granges2$GeneID %in% 
                                                          KS2_KS1_B_atac_promoters$GeneID)][
                                                            which(qobj_RT_KS2_KS1_B_atac_promoters$significant == TRUE)]



#####
common_promoters_KS1_KS2_only_for_csv <- common_promoters_KS1_KS2_only[, c("chr", "start", "end", 
                                                        "prom_overlapping_id", "prom_overlapping_name", "logFC_KS1", "logFC_KS2")]

#####
common_promoters <- annoGR2DF(RT_KS2_KS1_B_atac_promoters)
common_promoters$logFC_KS1 <- sapply(common_promoters$GeneID, function(xx) 
  res_B_KS1_atac$logFC[which(rownames(res_B_KS1_atac) == xx)])
common_promoters$logFC_KS2 <- sapply(common_promoters$GeneID, function(xx) 
  res_B_KS2_atac$logFC[which(rownames(res_B_KS2_atac) == xx)])
common_promoters$logFC_RT <- sapply(common_promoters$GeneID, function(xx) 
  res_B_RT_atac$logFC[which(rownames(res_B_RT_atac) == xx)])

gene_ids_to_use1 <- common_promoters$GeneID[which(common_promoters$logFC_KS1 < 0 & 
                                                    common_promoters$logFC_KS2 < 0 & 
                                                    common_promoters$logFC_RT < 0)]
gene_ids_to_use2 <- common_promoters$GeneID[which(common_promoters$logFC_KS1 > 0 & 
                                                    common_promoters$logFC_KS2 > 0 & 
                                                    common_promoters$logFC_RT > 0)]
common_promoters <- common_promoters[which(common_promoters$GeneID %in% c(gene_ids_to_use1, gene_ids_to_use2)), ]


common_promoters_for_csv <- common_promoters[, c("chr", "start", "end", 
                                        "prom_overlapping_id", "prom_overlapping_name", "logFC_KS1", "logFC_KS2", "logFC_RT")]


write_csv(common_promoters_KS1_KS2_only_for_csv, "common_promoters_B_KS1_KS2_only.csv")
write_csv(common_promoters_for_csv, "common_promoters_B_KS1_KS2_and_RT.csv")


####Figures
tab1_B_proms <- res_B_KS2_atac[which(rownames(res_B_KS2_atac) %in% atac_B_KS2_granges2$GeneID), ]
tab2_B_proms <- res_B_KS1_atac[which(rownames(res_B_KS1_atac) %in% atac_B_KS1_granges2$GeneID), ]
tab3_B_proms <- res_B_RT_atac[which(rownames(res_B_RT_atac) %in% atac_B_RT_granges2$GeneID), ]


quartz(file = "KS1_vs_KS2_atacseq_B_cells_promoters.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2), mar = c(5,4,1,1) + 0.1)
compareDEAnalyses(tab1_B_proms, tab2_B_proms, max(tab2_B_proms$P.Value[which(tab2_B_proms$padj < 0.1)]), 
                  "p-value_density", "ATAC-seq (promoters)", "KS2 p-values", "Density", 
                  "KS1 FDR < 0.1", "KS1 FDR >= 0.1", y_lim = c(0, 6), bwidth = 0.035)

compareDEAnalyses(tab1_B_proms, tab2_B_proms, max(tab2_B_proms$P.Value[which(tab2_B_proms$padj < 0.1)]), "log2FC", 
                  "ATAC-seq (promoters)", "KS1 logFC (KS1 FDR < 0.1)", "KS2 logFC (KS1 FDR < 0.1)", 
                  "", "", x_lim = c(-1.5, 2), y_lim = c(-1.5, 1.5))

dev.off()

quartz(file = "RT_vs_KS1_and_2_common_hits_atacseq_B_cells_promoters.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2), mar = c(5,4,1,1) + 0.1)
compareDEAnalyses(tab3_B_proms, tab1_2_B_common, max(tab1_2_B_common$P.Value), 
                  "p-value_density", "ATAC-seq (promoters)", "RT1 p-values", "Density", 
                  "KS1/2 common hits", "other", y_lim = c(0, 7), bwidth = 0.025)

compareDEAnalyses(tab3_B_proms, tab1_2_B_common, max(tab1_2_B_common$P.Value), "log2FC", 
                  "ATAC-seq (promoters)", "common KS1/2 hits\nmean(logFC)", "RT1 logFC", 
                  "", "", x_lim = c(-1, 1.5), y_lim = c(-1, 2.25))

dev.off()


###PCA from common promoters



quartz(file = "pca_from_common_genes_RT_atac_promoters.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
load(file = "rt_atacseq/count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_rt_cohort.rda")
peaks_by_samples <- getPeaksBySamplesMatrix("rt_atacseq/count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_rt_cohort.rda", 
                                            "rt_blood")

sample_info <- data.frame(genotype = getGenotypeVector(colnames(peaks_by_samples)), 
                          cell_type = getCellTypeVector(colnames(peaks_by_samples)))
sampleTable <- data.frame(condition = factor(sample_info$genotype[which(sample_info$cell_type == "B")]))
dds <- DESeqDataSetFromMatrix(round(peaks_by_samples[, grep("_B", colnames(peaks_by_samples))]), sampleTable, design = ~condition)
vstab <- vst(dds)

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
     main = "", pch = 19, xlim = c(-4.5, 3.5), ylim = c(-1.5, 2.5), bty = 'l', xaxt = 'n', yaxt = 'n', cex.lab = 1.25)
points(pcaData$PC1[which(pcaData$genotype == "WT")], pcaData$PC2[which(pcaData$genotype == "WT")], 
       cex = 1.15, col = rgb(0, 0, 0, 0.8), pch = 19)
legend <- legend("top", legend = c("RT1", "WT"), bty = 'n', cex = 0.82, col = c("deep pink", rgb(0, 0, 0, 0.8)), pch = 19)
axis(1, at = c(-4,3), cex.axis = 1.2)
axis(2, at = c(-1.5, 2.5), cex.axis = 1.2)
dev.off()

quartz(file = "pca_from_common_genes_KS1_KS2_atac_promoters.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
load(file = "rt_atacseq/count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_kabuki_cohort.rda")
peaks_by_samples <- getPeaksBySamplesMatrix("rt_atacseq/count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_kabuki_cohort.rda", 
                                            "blood_all_samples")

sample_info <- data.frame(genotype = getGenotypeVector(colnames(peaks_by_samples)), 
                          cell_type = getCellTypeVector(colnames(peaks_by_samples)))
sampleTable <- data.frame(condition = factor(sample_info$genotype[which(sample_info$cell_type == "B")]))
dds <- DESeqDataSetFromMatrix(round(peaks_by_samples[, grep("_B_", colnames(peaks_by_samples))]), sampleTable, design = ~condition)
vstab <- vst(dds)


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
     main = "", pch = 19, xlim = c(min(pcaData$PC1)+0.1, max(pcaData$PC1)+0.1), 
     ylim = c(-2.25, max(pcaData$PC2)+0.1), bty = 'l', xaxt = 'n', yaxt = 'n', cex.lab = 1.25)
points(pcaData$PC1[which(pcaData$genotype == "WT")], pcaData$PC2[which(pcaData$genotype == "WT")], 
       cex = 1.15, col = rgb(0, 0, 0, 0.8), pch = 19)
points(pcaData$PC1[which(pcaData$genotype == "KS2")], pcaData$PC2[which(pcaData$genotype == "KS2")], 
       cex = 1.15, col = "forest green", pch = 19)
legend <- legend("topleft", legend = c("KS1", "KS2", "WT"), bty = 'n', cex = 0.82, col = c("orange", "forest green", 
                                                                                           rgb(0, 0, 0, 0.8)), pch = 19)
axis(1, at = c(-6, 5), cex.axis = 1.2)
axis(2, at = c(-2, 3), cex.axis = 1.2)
dev.off()

quartz(file = "pca_from_common_genes_KS2_atac_promoters.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")

#ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
#  geom_point(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#  coord_fixed() + scale_fill_manual(values=c("WT"="blue", "R172K"="red")) + labs(fill="Type", title = "RNA-seq") + 
#  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme_classic()

plot(pcaData$PC1[which(pcaData$genotype == "KS2")], pcaData$PC2[which(pcaData$genotype == "KS2")], 
     cex = 1.15, col = "forest green", xlab = paste0("PC1 (", percentVar[1], "%)"), ylab = paste0("PC2 (", percentVar[2], "%)"), 
     main = "", pch = 19, xlim = c(-2, max(pcaData$PC1)+0.1), 
     ylim = c(min(pcaData$PC2)+0.1, max(pcaData$PC2)+0.1), bty = 'l', xaxt = 'n', yaxt = 'n', cex.lab = 1.25)
points(pcaData$PC1[which(pcaData$genotype == "WT")], pcaData$PC2[which(pcaData$genotype == "WT")], 
       cex = 1.15, col = rgb(0, 0, 0, 0.8), pch = 19)
legend <- legend("top", legend = c("KS2", "WT"), bty = 'n', cex = 0.82, col = c("forest green", rgb(0, 0, 0, 0.8)), pch = 19)
axis(1, at = c(-2,5), cex.axis = 1.2)
axis(2, at = c(-2, 3), cex.axis = 1.2)
dev.off()

###non promoters
atac_B_KS1_granges3 <- atac_B_KS1_granges[which(is.na(atac_B_KS1_granges$prom_overlapping_id))]
atac_B_RT_granges3 <- atac_B_RT_granges[which(is.na(atac_B_RT_granges$prom_overlapping_id))]
atac_B_KS2_granges3 <- atac_B_KS2_granges[which(is.na(atac_B_KS2_granges$prom_overlapping_id))]



qobj_KS2_KS1_B_atac_non_promoters <- qvalue(p = atac_B_KS2_granges3$P.Value[which(atac_B_KS2_granges3$GeneID %in% 
                                                                      atac_B_KS1_granges3$GeneID[which(atac_B_KS1_granges3$padj < 0.1)])], 
                              fdr.level = 0.1, pi0.method = "bootstrap")
KS2_KS1_B_atac_non_promoters <- atac_B_KS2_granges3[which(atac_B_KS2_granges3$GeneID %in% 
                                                        atac_B_KS1_granges3$GeneID[which(atac_B_KS1_granges3$padj < 0.1)])][
                                                          which(qobj_KS2_KS1_B_atac_non_promoters$significant == TRUE)]


common_non_promoters_KS1_KS2_only <- annoGR2DF(KS2_KS1_B_atac_non_promoters)
common_non_promoters_KS1_KS2_only$logFC_KS1 <- sapply(common_non_promoters_KS1_KS2_only$GeneID, function(xx) 
  res_B_KS1_atac$logFC[which(rownames(res_B_KS1_atac) == xx)])
common_non_promoters_KS1_KS2_only$logFC_KS2 <- sapply(common_non_promoters_KS1_KS2_only$GeneID, function(xx) 
  res_B_KS2_atac$logFC[which(rownames(res_B_KS2_atac) == xx)])
common_non_promoters_KS1_KS2_only <- common_non_promoters_KS1_KS2_only[
  which(common_non_promoters_KS1_KS2_only$logFC_KS1*common_non_promoters_KS1_KS2_only$logFC_KS2 > 0), ]

KS2_KS1_B_atac_non_promoters <- KS2_KS1_B_atac_non_promoters[which(KS2_KS1_B_atac_non_promoters$GeneID %in% common_non_promoters_KS1_KS2_only$GeneID)]


tab1_2_B_common_2 <- res_B_KS1_atac[which(rownames(res_B_KS1_atac) %in% KS2_KS1_B_atac_non_promoters$GeneID), ]
tab1_2_B_common_2$logFC <- sapply(rownames(tab1_2_B_common_2), 
                                function(xx) mean(c(res_B_KS1_atac$logFC[which(rownames(res_B_KS1_atac) == xx)], 
                                                    res_B_KS2_atac$logFC[which(rownames(res_B_KS2_atac) == xx)])))


qobj_RT_KS2_KS1_B_atac_non_promoters <- qvalue(p = atac_B_RT_granges3$P.Value[which(atac_B_RT_granges3$GeneID %in% 
                                                                                  KS2_KS1_B_atac_non_promoters$GeneID)], 
                                           fdr.level = 0.1, pi0.method = "bootstrap")

RT_KS2_KS1_B_atac_non_promoters <- atac_B_RT_granges3[which(atac_B_RT_granges3$GeneID %in% 
                                                          KS2_KS1_B_atac_non_promoters$GeneID)][
                                                            which(qobj_RT_KS2_KS1_B_atac_non_promoters$significant == TRUE)]

#####
common_non_promoters <- annoGR2DF(RT_KS2_KS1_B_atac_non_promoters)
common_non_promoters$logFC_KS1 <- sapply(common_non_promoters$GeneID, function(xx) 
  res_B_KS1_atac$logFC[which(rownames(res_B_KS1_atac) == xx)])
common_non_promoters$logFC_KS2 <- sapply(common_non_promoters$GeneID, function(xx) 
  res_B_KS2_atac$logFC[which(rownames(res_B_KS2_atac) == xx)])
common_non_promoters$logFC_RT <- sapply(common_non_promoters$GeneID, function(xx) 
  res_B_RT_atac$logFC[which(rownames(res_B_RT_atac) == xx)])

gene_ids_to_use1_non_promoters <- common_non_promoters$GeneID[which(common_non_promoters$logFC_KS1 < 0 & 
                                                    common_non_promoters$logFC_KS2 < 0 & 
                                                    common_non_promoters$logFC_RT < 0)]
gene_ids_to_use2_non_promoters <- common_non_promoters$GeneID[which(common_non_promoters$logFC_KS1 > 0 & 
                                                    common_non_promoters$logFC_KS2 > 0 & 
                                                    common_non_promoters$logFC_RT > 0)]
common_non_promoters <- common_non_promoters[which(common_non_promoters$GeneID %in% 
                                                     c(gene_ids_to_use1_non_promoters, gene_ids_to_use2_non_promoters)), ]


common_non_promoters_for_csv <- common_non_promoters[, c("chr", "start", "end", 
                                                 "prom_overlapping_id", "prom_overlapping_name", "logFC_KS1", "logFC_KS2", "logFC_RT")]


common_non_promoters_KS1_KS2_only_for_csv <- 

write_csv(common_non_promoters_KS1_KS2_only_for_csv, "common_non_promoters_B_KS1_KS2_only.csv")
write_csv(common_non_promoters_for_csv, "common_non_promoters_B_KS1_KS2_and_RT.csv")


####Figures
tab1_B_no_proms <- res_B_KS2_atac[-which(rownames(res_B_KS2_atac) %in% atac_B_KS2_granges2$GeneID), ]
tab2_B_no_proms <- res_B_KS1_atac[-which(rownames(res_B_KS1_atac) %in% atac_B_KS1_granges2$GeneID), ]
tab3_B_no_proms <- res_B_RT_atac[-which(rownames(res_B_RT_atac) %in% atac_B_RT_granges2$GeneID), ]




###Figures
quartz(file = "KS1_vs_KS2_atacseq_B_cells_non_promoters.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2), mar = c(5,4,1,1) + 0.1)
compareDEAnalyses(tab1_B_no_proms, tab2_B_no_proms, max(tab2_B_no_proms$P.Value[which(tab2_B_no_proms$padj < 0.1)]), 
                  "p-value_density", "ATAC-seq (distal reg. elem.)", "KS2 p-values", "Density", 
                  "KS1 FDR < 0.1", "KS1 FDR >= 0.1", y_lim = c(0, 4), bwidth = 0.035)

compareDEAnalyses(tab1_B_no_proms, tab2_B_no_proms, max(tab2_B_no_proms$P.Value[which(tab2_B_no_proms$padj < 0.1)]), "log2FC", 
                  "ATAC-seq (distal reg. elem.)", "KS1 logFC (KS1 FDR < 0.1)", "KS2 logFC (KS1 FDR < 0.1)", 
                  "", "", x_lim = c(-2, 2), y_lim = c(-1.5, 1.5))

dev.off()

quartz(file = "RT_vs_KS1_and_2_common_hits_atacseq_B_cells_non_promoters.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2), mar = c(5,4,1,1) + 0.1)
compareDEAnalyses(tab3_B_no_proms, tab1_2_B_common_2, max(tab1_2_B_common_2$P.Value), 
                  "p-value_density", "ATAC-seq (distal reg. elem.)", "RT1 p-values", "Density", 
                  "KS1/2 common hits", "other", y_lim = c(0, 3), bwidth = 0.035)

compareDEAnalyses(tab3_B_no_proms, tab1_2_B_common_2, max(tab1_2_B_common_2$P.Value), "log2FC", 
                  "ATAC-seq (distal reg. elem.)", "common KS1/2 hits\nmean(logFC)", "RT1 logFC", 
                  "", "", x_lim = c(-1.2, 1.2), y_lim = c(-2,2))

dev.off()

####pairwise comparisons






###intersection between common DA promoters and common DE genes
n11 <- length(intersect(unlist(common_promoters$prom_overlapping_id), common_genes_B$gene_id))
n12 <- length(common_genes_B$gene_id) - n11
n21 <- length(unlist(common_promoters$prom_overlapping_id)) - n11
n22 <- length(intersect(Reduce(intersect, list(unique(atac_B_KS1_granges2$prom_overlapping_id), 
                                               unique(atac_B_KS1_granges2$prom_overlapping_id), 
                                               unique(atac_B_KS1_granges2$prom_overlapping_id))), 
                        Reduce(intersect, list(rownames(res_B_KS1), rownames(res_B_KS2), rownames(res_B_RT))))) - n12 - n21 - n11

fisher.test(matrix(c(n11,n12,n21,n22), nrow = 2))




###enrichment of common hits in promoters
all_proms_tested <- length(Reduce(intersect, list(atac_B_KS1_granges2$GeneID, atac_B_KS2_granges2$GeneID, atac_B_RT_granges2$GeneID)))
all_dres_tested <- length(Reduce(intersect, list(atac_B_KS1_granges3$GeneID, atac_B_KS2_granges3$GeneID, atac_B_RT_granges3$GeneID)))

n11 <- dim(common_promoters)[1]
n12 <- all_proms_tested - n11
n21 <- dim(common_non_promoters)[1]
n22 <- all_dres_tested - n21
fisher.test(matrix(c(n11,n12,n21,n22), nrow = 2))


####permutation distributions for pi0 p values
permutation_dist_RT_KS1_KS2_proms <- replicate(10000, {
  pi0_random <- pi0est(p = sample(atac_B_RT_granges2$P.Value, length(which(atac_B_RT_granges2$GeneID %in% 
                                                                                    KS2_KS1_B_atac_promoters$GeneID))), 
                        pi0.method = "bootstrap")$pi0
  1- pi0_random
  #length(which(qobj_random$significant == TRUE))
})
length(which(permutation_dist_RT_KS1_KS2_proms >= 1-qobj_RT_KS2_KS1_B_atac_promoters$pi0))/length(permutation_dist_RT_KS1_KS2_proms)
#
permutation_dist_KS1_KS2_proms <- replicate(10000, {
  pi0_random <- pi0est(p = sample(atac_B_KS2_granges2$P.Value, 
                                  length(which(atac_B_KS2_granges2$GeneID 
                                               %in% atac_B_KS1_granges2$GeneID[which(atac_B_KS1_granges2$padj < 0.1)]))), 
                       pi0.method = "bootstrap")$pi0
  1- pi0_random
  #length(which(qobj_random$significant == TRUE))
})


length(which(permutation_dist_KS1_KS2_proms >= 1-qobj_KS2_KS1_B_atac_promoters$pi0))/length(permutation_dist_KS1_KS2_proms)
#
permutation_dist_RT_KS1_KS2_non_proms <- replicate(10000, {
  pi0_random <- pi0est(p = sample(atac_B_RT_granges3$P.Value, length(which(atac_B_RT_granges3$GeneID %in% 
                                                                             KS2_KS1_B_atac_non_promoters$GeneID))), 
                       pi0.method = "bootstrap")$pi0
  1- pi0_random
  #length(which(qobj_random$significant == TRUE))
})
length(which(permutation_dist_RT_KS1_KS2_non_proms >= 1-qobj_RT_KS2_KS1_B_atac_non_promoters$pi0))/length(permutation_dist_RT_KS1_KS2_non_proms)
#
permutation_dist_KS1_KS2_non_proms <- replicate(10000, {
  pi0_random <- pi0est(p = sample(atac_B_KS2_granges3$P.Value, 
                                  length(which(atac_B_KS2_granges3$GeneID 
                                               %in% atac_B_KS1_granges3$GeneID[which(atac_B_KS1_granges3$padj < 0.1)]))), 
                       pi0.method = "bootstrap")$pi0
  1- pi0_random
  #length(which(qobj_random$significant == TRUE))
})


length(which(permutation_dist_KS1_KS2_non_proms >= 1-qobj_KS2_KS1_B_atac_non_promoters$pi0))/length(permutation_dist_KS1_KS2_non_proms)



###compare empirical power with naive approach
length(Reduce(intersect, list(atac_B_KS1_granges2$GeneID[which(atac_B_KS1_granges2$padj < 0.1)], 
                       atac_B_KS2_granges2$GeneID[which(atac_B_KS2_granges2$padj < 0.1)], 
                       atac_B_RT_granges2$GeneID[which(atac_B_RT_granges2$padj < 0.1)])))

length(Reduce(intersect, list(atac_B_KS1_granges2$GeneID[which(p.adjust(atac_B_KS1_granges2$pvalue, method = "fdr") < 0.1)], 
                       atac_B_KS2_granges2$GeneID[which(p.adjust(atac_B_KS2_granges2$pvalue, method = "fdr") < 0.1)], 
                       atac_B_RT_granges2$GeneID[which(p.adjust(atac_B_RT_granges2$pvalue, method = "fdr") < 0.1)])))


###pathway enrichment analysis
promoter_pathways <- getReactomeEnrichedPathways(unique(common_promoters$prom_overlapping_id), 
          unique(unlist(Reduce(intersect, list(atac_B_KS1_granges2$prom_overlapping_id, 
                                               atac_B_KS2_granges2$prom_overlapping_id, atac_B_RT_granges2$prom_overlapping_id)))))

expression_pathways <- getReactomeEnrichedPathways(unique(common_genes_B$gene_id), 
          unique(unlist(Reduce(intersect, list(rownames(res_B_KS1), 
                                               rownames(res_B_KS2), rownames(res_B_RT))))))


####
quartz(file = "ATAC_promoters_vs_non_promoters_pairwise_comparisons.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4, 5.7, 2, 1) + 0.1)
plot(1, 1 - pi0est(p = atac_B_KS2_granges2$P.Value[which(atac_B_KS2_granges2$GeneID %in% 
                                                       atac_B_KS1_granges2$GeneID[which(atac_B_KS1_granges2$padj < 0.1)])], 
               pi0.method = "bootstrap")$pi0, pch = 19, 
     col = alpha("orange", 1), xlab = "", ylab = "% shared diff access\nregulatory elements", main = "", 
     xlim = c(0.8, 6.2), ylim = c(0, 0.8), xaxt = 'n', bty = 'l', yaxt = 'n', bty = "l", cex = 1.25, cex.lab = 1.2)
points(2, 1 - pi0est(p = atac_B_KS2_granges3$P.Value[which(atac_B_KS2_granges3$GeneID %in% 
                                                         atac_B_KS1_granges3$GeneID[which(atac_B_KS1_granges3$padj < 0.1)])], 
                 pi0.method = "bootstrap")$pi0, pch = 19, 
       col = alpha("forest green", 1), cex = 1.25)
points(3, 1 - pi0est(p = atac_B_KS1_granges2$P.Value[which(atac_B_KS1_granges2$GeneID %in% 
                                                         atac_B_RT_granges2$GeneID[which(atac_B_RT_granges2$padj < 0.1)])], 
                 pi0.method = "bootstrap")$pi0, pch = 19, 
       col = alpha("orange", 1), cex = 1.25)
points(4, 1 - pi0est(p = atac_B_KS1_granges3$P.Value[which(atac_B_KS1_granges3$GeneID %in% 
                                                         atac_B_RT_granges3$GeneID[which(atac_B_RT_granges3$padj < 0.1)])], 
                 pi0.method = "bootstrap")$pi0, pch = 19, 
       col = alpha("forest green", 1), cex = 1.25)
points(5, 1 - pi0est(p = atac_B_KS2_granges2$P.Value[which(atac_B_KS2_granges2$GeneID %in% 
                                                         atac_B_RT_granges2$GeneID[which(atac_B_RT_granges2$padj < 0.1)])], 
                 pi0.method = "bootstrap")$pi0, pch = 19, 
       col = alpha("orange", 1), cex = 1.25)
points(6, 1 - pi0est(p = atac_B_KS2_granges3$P.Value[which(atac_B_KS2_granges3$GeneID %in% 
                                                         atac_B_RT_granges3$GeneID[which(atac_B_RT_granges3$padj < 0.1)])], 
                 pi0.method = "bootstrap")$pi0, pch = 19, 
       col = alpha("forest green", 1), cex = 1.25)
#points(7, 1-qobj_RT_KS2_KS1_B_atac_promoters$pi0, pch = 19, 
#       col = alpha("brown", 0.75), xlab = "", ylab = "% shared differentially accessible\nregulatory elements", main = "")

#points(8, 1-qobj_RT_KS2_KS1_B_atac_non_promoters$pi0, pch = 19, 
#       col = alpha("forest green", 1), xlab = "", ylab = "% shared differentially accessible\nregulatory elements", main = "")

axis(1, at = c(1.5, 3.5, 5.5), labels = c("KS1/KS2", "KS1/RT", "KS2/RT"), cex.axis = 0.9, las = 2)
axis(2, at = c(0, 0.4, 0.8))
legend <- legend("topright", legend = c("promoters", "distal reg elements"), col = c(alpha("orange", 1), alpha("forest green", 1)), 
                 bty = 'n', pch = 19, cex = 0.57)
abline(v = c(2.5, 4.5), lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()



######
quartz(file = "ATAC_effect_direction_promoters.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(5.25, 5.25, 1.25, 1.25) + 0.1)
plot(1, prop.table(table(as.factor(atac_B_KS1_granges2$logFC[which(atac_B_KS1_granges2$padj < 0.1)] > 0)))[2], pch = 19, 
     col = alpha("red", 0.75), xlab = "", ylab = "% changes towards\ngreater accessibility", main = "promoters", font.main = 1,  
     xlim = c(0.8, 5.2), ylim = c(0, 1), xaxt = 'n', bty = 'l', yaxt = 'n', bty = "l", cex = 1.25, cex.lab = 1.2)
points(2, prop.table(table(as.factor(atac_B_KS2_granges2$logFC[which(atac_B_KS2_granges2$padj < 0.1)] > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(3, prop.table(table(as.factor(atac_B_RT_granges2$logFC[which(atac_B_RT_granges2$padj < 0.1)] > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(4, prop.table(table(as.factor(KS2_KS1_B_atac_promoters$logFC > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(5, prop.table(table(as.factor(RT_KS2_KS1_B_atac_promoters$logFC > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
#points(7, 1-qobj_RT_KS2_KS1_B_atac_promoters$pi0, pch = 19, 
#       col = alpha("brown", 0.75), xlab = "", ylab = "% shared differentially accessible\nregulatory elements", main = "")

#points(8, 1-qobj_RT_KS2_KS1_B_atac_non_promoters$pi0, pch = 19, 
#       col = alpha("forest green", 1), xlab = "", ylab = "% shared differentially accessible\nregulatory elements", main = "")

axis(1, at = c(1, 2, 3, 4, 5), labels = c("KS1", "KS2", "RT", "KS1/KS2\ncommon", "KS1/KS2/RT\ncommon"), cex.axis = 0.84, las = 2)
axis(2, at = c(0, 0.5, 1), labels = c(0, 50, 100))
abline(v = c(1.5, 2.5, 3.5, 4.5, 5.5), lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()






###
quartz(file = "ATAC_effect_direction_non_promoters.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(5.25, 5.25, 1.25, 1.25) + 0.1)
plot(1, prop.table(table(as.factor(atac_B_KS1_granges3$logFC[which(atac_B_KS1_granges3$padj < 0.1)] > 0)))[2], pch = 19, 
     col = alpha("red", 0.75), xlab = "", ylab = "% changes towards\ngreater accessibility", main = "distal reg elements", font.main = 1, 
     xlim = c(0.8, 5.2), ylim = c(0, 1), xaxt = 'n', bty = 'l', yaxt = 'n', bty = "l", cex = 1.25, cex.lab = 1.2)
points(2, prop.table(table(as.factor(atac_B_KS2_granges3$logFC[which(atac_B_KS2_granges3$padj < 0.1)] > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(3, prop.table(table(as.factor(atac_B_RT_granges3$logFC[which(atac_B_RT_granges3$padj < 0.1)] > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(4, prop.table(table(as.factor(KS2_KS1_B_atac_non_promoters$logFC > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(5, prop.table(table(as.factor(RT_KS2_KS1_B_atac_non_promoters$logFC > 0)))[2], pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
#points(7, 1-qobj_RT_KS2_KS1_B_atac_promoters$pi0, pch = 19, 
#       col = alpha("brown", 0.75), xlab = "", ylab = "% shared differentially accessible\nregulatory elements", main = "")

#points(8, 1-qobj_RT_KS2_KS1_B_atac_non_promoters$pi0, pch = 19, 
#       col = alpha("forest green", 1), xlab = "", ylab = "% shared differentially accessible\nregulatory elements", main = "")

axis(1, at = c(1, 2, 3, 4, 5), labels = c("KS1", "KS2", "RT", "KS1/KS2\ncommon", "KS1/KS2/RT\ncommon"), cex.axis = 0.84, las = 2)
axis(2, at = c(0, 0.5, 1), labels = c(0, 50, 100))
abline(v = c(1.5, 2.5, 3.5, 4.5, 5.5), lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()

###
###
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


