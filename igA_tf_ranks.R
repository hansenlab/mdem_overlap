igA_deficiency <- c("Adgrg3", "Aicda", "Atp11c", "Bach2", "Batf", "Bst1", "Btk", "Cd40", "Card11", 
                    "Ccl28", "Ccr6", "Ccr10", "Cd19", "Ell2", "Tgfbr2", "Cd40", "Cd40lg", "Cd79a", "Huwe1", 
                    "Slc39a10", "Cd79b", "Cd83", "Cd207", "Cyp7b1", "Dll1", "Ednrb", "H2az2", "Fas", "Ifng", 
                    "Ifngr1", "Il4", "Gcnt3", "Ighg1", "Map3k14", "Hoxc4", "Icos", "Igh", "Igha", "Ighm", 
                    "Il5ra", "Il6", "Il6st", "Il21r", "Il23a", "Lig4", "Lta", "Tnf", "Ltb", "Madcam1", 
                    "Nfkb1", "Nfkbia", "Nfkb2", "Niban3", "Parp14", "Pax5", "Pigr", "Pik3cd", "Pou2af1", 
                    "Ikzf3", "Pvr", "Rag2", "Rc3h1", "Commd10", "Rela", "Rif1", "Rnf31", "Sh3bp2", "Sharpin", 
                    "Sppl2a", "St6galnac2", "Stat6", "Tcrb", "Tnfrsf11a", "Tnfrsf13b", "Tnfsf13", "Usp12")
igA_peyers_combined <- igA_deficiency #the name reflects that initially we also included peyer's patch genes. We ultimately ended up not including them for the final analysis

rank_KS1 <- wilcox.test(res_B_KS1$P.Value[which(res_B_KS1$gene_symbol %in% igA_peyers_combined)], 
                        res_B_KS1$P.Value[-which(res_B_KS1$gene_symbol %in% igA_peyers_combined)])$statistic
rank_KS2 <- wilcox.test(res_B_KS2$P.Value[which(res_B_KS2$gene_symbol %in% igA_peyers_combined)], 
                        res_B_KS2$P.Value[-which(res_B_KS2$gene_symbol %in% igA_peyers_combined)])$statistic
rank_RT <- wilcox.test(res_B_RT$P.Value[which(res_B_RT$gene_symbol %in% igA_peyers_combined)], 
                       res_B_RT$P.Value[-which(res_B_RT$gene_symbol %in% igA_peyers_combined)])$statistic

length1 <- length(which(res_B_KS1$gene_symbol %in% igA_peyers_combined))
length2 <- length(which(res_B_KS2$gene_symbol %in% igA_peyers_combined))
length3 <- length(which(res_B_RT$gene_symbol %in% igA_peyers_combined))

permutation_rank_KS1 <- replicate(10000, {
  indices <- sample(1:length(rownames(res_B_KS1)), length1)
  wilcox.test(res_B_KS1$P.Value[indices], res_B_KS1$P.Value[-indices])$statistic
})

permutation_rank_KS2 <- replicate(10000, {
  indices <- sample(1:length(rownames(res_B_KS2)), length2)
  wilcox.test(res_B_KS2$P.Value[indices], res_B_KS2$P.Value[-indices])$statistic
})

permutation_rank_RT <- replicate(10000, {
  indices <- sample(1:length(rownames(res_B_RT)), length3)
  wilcox.test(res_B_RT$P.Value[indices], res_B_RT$P.Value[-indices])$statistic
})


quartz(file = "igA_ranks.pdf", width = 6.5, height = 2.2, pointsize = 8, type = "pdf")
par(mfrow = c(1,3))
hist(permutation_rank_KS1, col = "cornflowerblue", lty = 0, 
     breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat", cex.lab = 1.2, yaxt = 'n',
     main = "KS1", cex.main = 1.2, font.main = 1, xlim = c(min(permutation_rank_KS1)-0.1, max(permutation_rank_KS1)+0.1), xaxt = 'n')
axis(1, at = c(350000, 500000), cex.axis = 1.35)
axis(2, at = c(0, 0.000012), cex.axis = 1.35)
abline(v = rank_KS1, col = alpha("red", 0.6), lwd = 2.5)
legend <- legend("topright", legend = c("random", "observed"), col = c("cornflowerblue", alpha("red", 0.6)), bty = 'n', 
                 cex = 1.25, lty = "solid", lwd = 2.5)

hist(permutation_rank_KS2, col = "cornflowerblue", lty = 0, 
     breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat", cex.lab = 1.2, yaxt = 'n',
     main = "KS2", cex.main = 1.2, font.main = 1, xlim = c(min(permutation_rank_KS2)-0.1, max(permutation_rank_KS2)+0.1), xaxt = 'n')
axis(1, at = c(250000, 350000), cex.axis = 1.35)
axis(2, at = c(0, 0.000015), cex.axis = 1.35)
abline(v = rank_KS2, col = alpha("red", 0.6), lwd = 2.5)


hist(permutation_rank_RT, col = "cornflowerblue", lty = 0,
     breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat", cex.lab = 1.2, yaxt = 'n',
     main = "RT", cex.main = 1.2, font.main = 1, xlim = c(min(permutation_rank_RT)-0.1, max(permutation_rank_RT)+0.1), xaxt = 'n')
axis(1, at = c(300000, 450000), cex.axis = 1.35)
axis(2, at = c(0, 0.000012), cex.axis = 1.35)
abline(v = rank_RT, col = alpha("red", 0.6), lwd = 2.5)

dev.off()





rank_KS1_2 <- as.numeric(sapply(igA_peyers_combined, function(xx) 
  rank(res_B_KS1$P.Value)[which(res_B_KS1$gene_symbol == xx)]))
rank_RT_2 <- as.numeric(sapply(igA_peyers_combined, function(xx) 
  rank(res_B_RT$P.Value)[which(res_B_RT$gene_symbol == xx)]))
rank_KS1_2 <- rank_KS1[-which(is.na(rank_KS1_2))]
rank_RT_2 <- rank_RT[-which(is.na(rank_RT_2))]


quartz(file = "igA_ranks_comparisons.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4,5,1,1)+0.1)
plot(1, length(which(rank_KS1_2[which(rank_RT_2 <= quantile(rank_RT_2, 0.25))] <= quantile(rank_KS1_2, 0.25))) / 
       length(rank_KS1_2[which(rank_RT_2 <= quantile(rank_RT_2, 0.25))]), pch = 19, cex = 1.25, bty = 'l', 
     xlim = c(0.8, 4.2), ylim = c(0, 0.5), ylab = "% in top 25% diff expr\nIgA genes (KS1)", col = "orange", 
     xlab = "RT diff expr p-val quartile(IgA only)", xaxt = 'n', yaxt = 'n')


points(2, length(which(rank_KS1_2[which(rank_RT_2 > quantile(rank_RT_2, 0.25) & 
                                              rank_RT_2 <= quantile(rank_RT_2, 0.5))] <= quantile(rank_KS1_2, 0.25))) / 
         length(rank_KS1_2[which(rank_RT_2 > quantile(rank_RT_2, 0.25) & rank_RT_2 <= quantile(rank_RT_2, 0.5))]), pch = 19, cex = 1.25, 
       col = "orange")


points(3, length(which(rank_KS1_2[which(rank_RT_2 > quantile(rank_RT_2, 0.5) & 
                                              rank_RT_2 <= quantile(rank_RT_2, 0.75))] <= quantile(rank_KS1_2, 0.25))) / 
         length(rank_KS1_2[which(rank_RT_2 > quantile(rank_RT_2, 0.5) & rank_RT_2 <= quantile(rank_RT_2, 0.75))]), pch = 19, cex = 1.25, 
       col = "orange")


points(4, length(which(rank_KS1_2[which(rank_RT_2 > quantile(rank_RT_2, 0.75))] <= quantile(rank_KS1_2, 0.25))) / 
         length(rank_KS1_2[which(rank_RT_2 > quantile(rank_RT_2, 0.75))]), pch = 19, cex = 1.25, col = "orange")

axis(1, at = c(1, 2, 3, 4))
axis(2, at = c(0, 0.25, 0.5))
abline(v = c(1.5,2.5,3.5), lty = "longdash", col = rgb(0,0,0,0.7))

dev.off()




####TFs
tf_genes <- read_csv('tf_gene_list.csv')
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human_mouse_homologs <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                             "mmusculus_homolog_ensembl_gene", 
                                             "mmusculus_homolog_associated_gene_name", 
                                             "mmusculus_homolog_orthology_confidence", 
                                             "mmusculus_homolog_perc_id_r1"),
                              filters = "ensembl_gene_id",
                              values = unique(tf_genes$`Enseml Gene ID`),
                              mart = human)

human_mouse_homologs <- human_mouse_homologs[which(human_mouse_homologs$mmusculus_homolog_orthology_confidence == 1), ]

mouse_tfs <- unique(human_mouse_homologs$mmusculus_homolog_ensembl_gene)

mouse_tfs_combined <- mouse_tfs[which(mouse_tfs %in% rownames(res_B_KS1) & 
                                        mouse_tfs %in% rownames(res_B_KS2) & 
                                        mouse_tfs %in% rownames(res_B_RT))]

rank_KS1_tfs_2 <- sapply(mouse_tfs_combined, function(xx) 
  rank(res_B_KS1$P.Value)[which(rownames(res_B_KS1) == xx)])
rank_KS2_tfs_2 <- sapply(mouse_tfs_combined, function(xx) 
  rank(res_B_KS2$P.Value)[which(rownames(res_B_KS2) == xx)])
rank_RT_tfs_2 <- sapply(mouse_tfs_combined, function(xx) 
  rank(res_B_RT$P.Value)[which(rownames(res_B_RT) == xx)])


rank_KS1_tfs <- wilcox.test(res_B_KS1$P.Value[which(rownames(res_B_KS1) %in% mouse_tfs_combined)], 
                            res_B_KS1$P.Value[-which(rownames(res_B_KS1) %in% mouse_tfs_combined)])$statistic
rank_KS2_tfs <- wilcox.test(res_B_KS2$P.Value[which(rownames(res_B_KS2) %in% mouse_tfs_combined)], 
                            res_B_KS2$P.Value[-which(rownames(res_B_KS2) %in% mouse_tfs_combined)])$statistic
rank_RT_tfs <- wilcox.test(res_B_RT$P.Value[which(rownames(res_B_RT) %in% mouse_tfs_combined)], 
                           res_B_RT$P.Value[-which(rownames(res_B_RT) %in% mouse_tfs_combined)])$statistic


length1_tfs <- length(which(rownames(res_B_KS1) %in% mouse_tfs_combined))
length2_tfs <- length(which(rownames(res_B_KS2) %in% mouse_tfs_combined))
length3_tfs <- length(which(rownames(res_B_RT) %in% mouse_tfs_combined))

permutation_rank_KS1_tfs <- replicate(10000, {
  indices <- sample(1:length(rownames(res_B_KS1)), length1_tfs)
  wilcox.test(res_B_KS1$P.Value[indices], res_B_KS1$P.Value[-indices])$statistic
})

permutation_rank_KS2_tfs <- replicate(10000, {
  indices <- sample(1:length(rownames(res_B_KS2)), length2_tfs)
  wilcox.test(res_B_KS2$P.Value[indices], res_B_KS2$P.Value[-indices])$statistic
})

permutation_rank_RT_tfs <- replicate(10000, {
  indices <- sample(1:length(rownames(res_B_RT)), length3_tfs)
  wilcox.test(res_B_RT$P.Value[indices], res_B_RT$P.Value[-indices])$statistic
})

###
quartz(file = "tf_ranks.pdf", width = 6.5, height = 2.2, pointsize = 8, type = "pdf")
par(mfrow = c(1,3))
hist(permutation_rank_KS1_tfs, col = "cornflowerblue", lty = 0, xaxt = 'n', yaxt = 'n', font.main = 1, cex.lab = 1.2,
     breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat", 
     main = "KS1", cex.main = 1.2, xlim = c(rank_KS1_tfs, max(permutation_rank_KS1_tfs)+0.2))
axis(1, at = c(3100000, 3650000), cex.axis = 1.35)
axis(2, at = c(0, 0.000004), cex.axis = 1.35)
abline(v = rank_KS1_tfs, col = alpha("red", 0.6), lwd = 2.5)
legend <- legend("topright", legend = c("random", "observed"), col = c("cornflowerblue", alpha("red", 0.6)), bty = 'n', 
                 cex = 1.25, lty = "solid", lwd = 2.5)

hist(permutation_rank_KS2_tfs, col = "cornflowerblue", lty = 0, xaxt = 'n', yaxt = 'n', font.main = 1, cex.lab = 1.2,
     breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat", 
     main = "KS2", cex.main = 1.2, xlim = c(rank_KS2_tfs, max(permutation_rank_KS2_tfs)+0.2))
axis(1, at = c(2500000, 2900000), cex.axis = 1.35)
axis(2, at = c(0, 0.000004), cex.axis = 1.35)
abline(v = rank_KS2_tfs, col = alpha("red", 0.6), lwd = 2.5)

hist(permutation_rank_RT_tfs, col = "cornflowerblue", lty = 0,xaxt = 'n', yaxt = 'n', font.main = 1, cex.lab = 1.2,
     breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat",
     main = "RT", cex.main = 1.2, xlim = c(rank_RT_tfs, max(permutation_rank_RT_tfs)+0.2))
axis(1, at = c(3000000, 3500000), cex.axis = 1.35)
axis(2, at = c(0, 0.000004), cex.axis = 1.35)
abline(v = rank_RT_tfs, col = alpha("red", 0.6), lwd = 2.5)

dev.off()

quartz(file = "tf_ranks_comparisons.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4,5,1,1)+0.1)
plot(1, length(which(rank_KS1_tfs_2[which(rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.25))] <= quantile(rank_KS1_tfs_2, 0.25))) / 
       length(rank_KS1_tfs_2[which(rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.25))]), pch = 19, cex = 1.25, bty = 'l', 
     xlim = c(0.8, 8.2), ylim = c(0, 0.4), ylab = "% in top 25% diff expr\nTF genes", col = "orange", 
     xlab = "RT diff expr p-val quartile(TFs only)", xaxt = 'n', yaxt = 'n')

points(2, length(which(rank_KS2_tfs_2[which(rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.25))] <= quantile(rank_KS2_tfs_2, 0.25))) / 
         length(rank_KS2_tfs_2[which(rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.25))]), pch = 19, cex = 1.25, 
       col = "forest green")

points(3, length(which(rank_KS1_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.25) & 
                                              rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.5))] <= quantile(rank_KS1_tfs_2, 0.25))) / 
         length(rank_KS1_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.25) & rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.5))]), pch = 19, cex = 1.25, 
       col = "orange")

points(4, length(which(rank_KS2_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.25) & 
                                              rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.5))] <= quantile(rank_KS2_tfs_2, 0.25))) / 
         length(rank_KS2_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.25) & rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.5))]), pch = 19, cex = 1.25, 
       col = "forest green")

points(5, length(which(rank_KS1_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.5) & 
                                              rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.75))] <= quantile(rank_KS1_tfs_2, 0.25))) / 
         length(rank_KS1_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.5) & rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.75))]), pch = 19, cex = 1.25, 
       col = "orange")

points(6, length(which(rank_KS2_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.25) & 
                                              rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.5))] <= quantile(rank_KS2_tfs_2, 0.25))) / 
         length(rank_KS2_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.25) & rank_RT_tfs_2 <= quantile(rank_RT_tfs_2, 0.5))]), pch = 19, cex = 1.25, 
       col = "forest green")

points(7, length(which(rank_KS1_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.75))] <= quantile(rank_KS1_tfs_2, 0.25))) / 
         length(rank_KS1_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.75))]), pch = 19, cex = 1.25, col = "orange")

points(8, length(which(rank_KS2_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.75))] <= quantile(rank_KS2_tfs_2, 0.25))) / 
         length(rank_KS2_tfs_2[which(rank_RT_tfs_2 > quantile(rank_RT_tfs_2, 0.75))]), pch = 19, cex = 1.25, col = "forest green")
axis(1, at = c(1.5,3.5,5.5,7.5), labels = c(1, 2, 3, 4))
axis(2, at = c(0, 0.2, 0.4))
abline(v = c(2.5,4.5,6.5), lty = "longdash", col = rgb(0,0,0,0.7))
legend("bottomleft", legend = c("KS1", "KS2"), col = c("orange", "forest green"), pch = 19, bty = 'n', cex = 0.75)

dev.off()



#####ATAC promoter ranks
rank_KS1 <- wilcox.test(atac_B_KS1_granges2$P.Value[which(atac_B_KS1_granges2$prom_overlapping_name %in% igA_peyers_combined)], 
                        atac_B_KS1_granges2$P.Value[-which(atac_B_KS1_granges2$prom_overlapping_name %in% igA_peyers_combined)])$statistic
rank_KS2 <- wilcox.test(atac_B_KS2_granges2$P.Value[which(atac_B_KS2_granges2$prom_overlapping_name %in% igA_peyers_combined)], 
                        atac_B_KS2_granges2$P.Value[-which(atac_B_KS2_granges2$prom_overlapping_name %in% igA_peyers_combined)])$statistic
rank_RT <- wilcox.test(atac_B_RT_granges2$P.Value[which(atac_B_RT_granges2$prom_overlapping_name %in% igA_peyers_combined)], 
                       atac_B_RT_granges2$P.Value[-which(atac_B_RT_granges2$prom_overlapping_name %in% igA_peyers_combined)])$statistic

length1 <- length(which(atac_B_KS1_granges2$prom_overlapping_name %in% igA_peyers_combined))
length2 <- length(which(atac_B_KS2_granges2$prom_overlapping_name %in% igA_peyers_combined))
length3 <- length(which(atac_B_RT_granges2$prom_overlapping_name %in% igA_peyers_combined))

permutation_rank_KS1 <- replicate(10000, {
  indices <- sample(1:length(atac_B_KS1_granges2), length1)
  wilcox.test(atac_B_KS1_granges2$P.Value[indices], atac_B_KS1_granges2$P.Value[-indices])$statistic
})

permutation_rank_KS2 <- replicate(10000, {
  indices <- sample(1:length(atac_B_KS2_granges2), length2)
  wilcox.test(atac_B_KS2_granges2$P.Value[indices], atac_B_KS2_granges2$P.Value[-indices])$statistic
})

permutation_rank_RT <- replicate(10000, {
  indices <- sample(1:length(atac_B_RT_granges2), length3)
  wilcox.test(atac_B_RT_granges2$P.Value[indices], atac_B_RT_granges2$P.Value[-indices])$statistic
})




####
rank_KS1_tfs <- wilcox.test(atac_B_KS1_granges2$P.Value[which(atac_B_KS1_granges2$prom_overlapping_id %in% igA_peyers_combined)], 
                        atac_B_KS1_granges2$P.Value[-which(atac_B_KS1_granges2$prom_overlapping_name %in% igA_peyers_combined)])$statistic
rank_KS2_tfs <- wilcox.test(atac_B_KS2_granges2$P.Value[which(atac_B_KS2_granges2$prom_overlapping_name %in% igA_peyers_combined)], 
                        atac_B_KS2_granges2$P.Value[-which(atac_B_KS2_granges2$prom_overlapping_name %in% igA_peyers_combined)])$statistic
rank_RT_tfs <- wilcox.test(atac_B_RT_granges2$P.Value[which(atac_B_RT_granges2$prom_overlapping_name %in% igA_peyers_combined)], 
                       atac_B_RT_granges2$P.Value[-which(atac_B_RT_granges2$prom_overlapping_name %in% igA_peyers_combined)])$statistic

length1_tfs <- length(which(atac_B_KS1_granges2$prom_overlapping_name %in% igA_peyers_combined))
length2_tfs <- length(which(atac_B_KS2_granges2$prom_overlapping_name %in% igA_peyers_combined))
length3_tfs <- length(which(atac_B_RT_granges2$prom_overlapping_name %in% igA_peyers_combined))

permutation_rank_KS1_tfs <- replicate(10000, {
  indices <- sample(1:length(atac_B_KS1_granges2), length1_tfs)
  wilcox.test(atac_B_KS1_granges2$P.Value[indices], atac_B_KS1_granges2$P.Value[-indices])$statistic
})

permutation_rank_KS2_tfs <- replicate(10000, {
  indices <- sample(1:length(atac_B_KS2_granges2), length2_tfs)
  wilcox.test(atac_B_KS2_granges2$P.Value[indices], atac_B_KS2_granges2$P.Value[-indices])$statistic
})

permutation_rank_RT_tfs <- replicate(10000, {
  indices <- sample(1:length(atac_B_RT_granges2), length3_tfs)
  wilcox.test(atac_B_RT_granges2$P.Value[indices], atac_B_RT_granges2$P.Value[-indices])$statistic
})










