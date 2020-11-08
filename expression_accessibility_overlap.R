##########################

pi0_vector_KS1 <- sapply(seq(1000, 5000, by = 250), function(xx) 
  pi0est(res_B_KS1$P.Value[which(rownames(res_B_KS1) %in% atac_B_KS1_granges2$prom_overlapping_id[
    which(rank(atac_B_KS1_granges2$P.Value) < xx)])], pi0.method = "bootstrap")$pi0)

pi0_vector_KS2 <- sapply(seq(1000, 5000, by = 250), function(xx) 
  pi0est(res_B_KS2$P.Value[which(rownames(res_B_KS2) %in% atac_B_KS2_granges2$prom_overlapping_id[
    which(rank(atac_B_KS2_granges2$P.Value) < xx)])], pi0.method = "bootstrap")$pi0)

pi0_vector_RT <- sapply(seq(1000, 5000, by = 250), function(xx) 
  pi0est(res_B_RT$P.Value[which(rownames(res_B_RT) %in% atac_B_RT_granges2$prom_overlapping_id[
    which(rank(atac_B_RT_granges2$P.Value) < xx)])], pi0.method = "bootstrap")$pi0)

quartz(file = "accessibility_expression_comparison.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4, 4, 1, 1)+0.1)
plot(1:17, 1-pi0_vector_KS1, pch = 19, col = "orange", bty = 'l', ylab = "% diff expr genes", main = "", xaxt = 'n', yaxt = 'n', 
     xlim = c(0.8, 54), ylim = c(0.15, 0.52), xlab = "promoter diff accessibility ranking", cex = 0.84)
points(19:35, 1-pi0_vector_KS2, pch = 19, col = "forest green", cex = 0.84)
points(37:53, 1-pi0_vector_RT, pch = 19, col = "deep pink", cex = 0.84)
#abline(v = c(17.5, 34.5), lty = "longdash", col = rgb(0,0,0,0.7))
axis(2, at = c(0.20, 0.35, 0.50), labels = c(20, 35, 50))
text(1, 1-pi0_vector_KS1[1], "top\n1000", cex = 0.84, pos = 4)
text(17, 1-pi0_vector_KS1[17], "top\n5000", cex = 0.84, pos = 4)
text(18, 1-pi0_vector_KS2[1], "top\n1000", cex = 0.84, pos = 4)
text(34, 1-pi0_vector_KS2[17], "top\n5000", cex = 0.84, pos = 4)
text(35, 1-pi0_vector_RT[1], "top\n1000", cex = 0.84, pos = 4)
text(51, 1-pi0_vector_RT[17], "top\n5000", cex = 0.84, pos = 1)
legend <- legend("topleft", legend = c("KS1", "KS2", "RT"), pch = 19, bty = 'n', cex = 0.57, 
                 col = c("orange", "forest green", "deep pink"))
dev.off()


####
top_1000_proms_KS1_de_genes <- qvalue(res_B_KS1$P.Value[which(rownames(res_B_KS1) %in% atac_B_KS1_granges2$prom_overlapping_id[
  which(rank(atac_B_KS1_granges2$P.Value) < 1000)])], pi0.method = "bootstrap", fdr.level = 0.1)
top_1000_proms_KS1_de_genes <- res_B_KS1[which(rownames(res_B_KS1) %in% atac_B_KS1_granges2$prom_overlapping_id[
  which(rank(atac_B_KS1_granges2$P.Value) < 1000)])[which(top_1000_proms_KS1_de_genes$significant == TRUE)], ]

top_1000_proms_KS2_de_genes <- qvalue(res_B_KS2$P.Value[which(rownames(res_B_KS2) %in% atac_B_KS2_granges2$prom_overlapping_id[
  which(rank(atac_B_KS2_granges2$P.Value) < 1000)])], pi0.method = "bootstrap", fdr.level = 0.1)
top_1000_proms_KS2_de_genes <- res_B_KS2[which(rownames(res_B_KS2) %in% atac_B_KS2_granges2$prom_overlapping_id[
  which(rank(atac_B_KS2_granges2$P.Value) < 1000)])[which(top_1000_proms_KS2_de_genes$significant == TRUE)], ]

top_1000_proms_RT_de_genes <- qvalue(res_B_RT$P.Value[which(rownames(res_B_RT) %in% atac_B_RT_granges2$prom_overlapping_id[
  which(rank(atac_B_RT_granges2$P.Value) < 1000)])], pi0.method = "bootstrap", fdr.level = 0.1)
top_1000_proms_RT_de_genes <- res_B_RT[which(rownames(res_B_RT) %in% atac_B_RT_granges2$prom_overlapping_id[
  which(rank(atac_B_RT_granges2$P.Value) < 1000)])[which(top_1000_proms_RT_de_genes$significant == TRUE)], ]


proms_logFC_KS1 <- sapply(rownames(top_1000_proms_KS1_de_genes), function(xx) 
  median(atac_B_KS1_granges2$logFC[grep(xx, atac_B_KS1_granges2$prom_overlapping_id)]))

proms_logFC_KS2 <- sapply(rownames(top_1000_proms_KS2_de_genes), function(xx) 
  median(atac_B_KS2_granges2$logFC[grep(xx, atac_B_KS2_granges2$prom_overlapping_id)]))

proms_logFC_RT <- sapply(rownames(top_1000_proms_RT_de_genes), function(xx) 
  median(atac_B_RT_granges2$logFC[grep(xx, atac_B_RT_granges2$prom_overlapping_id)]))

quartz(file = "expression_vs_accessibility_logFC.pdf", height = 2.2, width = 6.59, type = "pdf")
par(mfrow = c(1, 3), mar = c(5,4,1,1) + 0.1)
plot(proms_logFC_KS1, 
     top_1000_proms_KS1_de_genes$logFC, 
     pch = 19, col = alpha("red", 0.57), 
     xlab = "promoter accessibility logFC", ylab = "gene expression logFC", main = "KS1", bty = 'l', yaxt = 'n', xaxt = 'n', 
     font.main = 1, ylim = c(-3, 3))
abline(h = 0, lty = "longdash", col = rgb(0,0,0,0.7))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
axis(2, at = c(-3, 3))
axis(1, at = c(-0.5, 0.5))

plot(proms_logFC_KS2, 
     top_1000_proms_KS2_de_genes$logFC, 
     pch = 19, col = alpha("red", 0.57), 
     xlab = "promoter accessibility logFC", ylab = "gene expression logFC", main = "KS2", bty = 'l', yaxt = 'n', xaxt = 'n', 
     font.main = 1, ylim = c(-1.5, 2))
abline(h = 0, lty = "longdash", col = rgb(0,0,0,0.7))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
axis(2, at = c(-1.5, 2))
axis(1, at = c(-0.5, 0.75))

plot(proms_logFC_RT, 
     top_1000_proms_RT_de_genes$logFC, 
     pch = 19, col = alpha("red", 0.57), 
     ylab = "promoter accessibility logFC", xlab = "gene expression logFC", main = "RT", bty = 'l', yaxt = 'n', xaxt = 'n', 
     font.main = 1, ylim = c(-2, 3))
abline(h = 0, lty = "longdash", col = rgb(0,0,0,0.7))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
axis(2, at = c(-2, 3))
axis(1, at = c(-0.75, 1.75))
dev.off()


###unique vs common hits
pi0_vector_KS1 <- sapply(seq(1000, 5000, by = 250), function(xx) 
  pi0est(res_B_KS1$P.Value[which(rownames(res_B_KS1) %in% atac_B_KS1_granges2$prom_overlapping_id[
    which(rank(atac_B_KS1_granges2$P.Value) < xx)] & !(rownames(res_B_KS1) %in% common_ids))], pi0.method = "bootstrap")$pi0)

pi0_vector_KS2 <- sapply(seq(1000, 5000, by = 250), function(xx) 
  pi0est(res_B_KS2$P.Value[which(rownames(res_B_KS2) %in% atac_B_KS2_granges2$prom_overlapping_id[
    which(rank(atac_B_KS2_granges2$P.Value) < xx)] & !(rownames(res_B_KS2) %in% common_ids))], pi0.method = "bootstrap")$pi0)

pi0_vector_RT <- sapply(seq(1000, 5000, by = 250), function(xx) 
  pi0est(res_B_RT$P.Value[which(rownames(res_B_RT) %in% atac_B_RT_granges2$prom_overlapping_id[
    which(rank(atac_B_RT_granges2$P.Value) < xx)] & !(rownames(res_B_RT) %in% common_ids))], pi0.method = "bootstrap")$pi0)

quartz(file = "accessibility_expression_comparison_2.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(1, 1-pi0_vector_KS1[1], pch = 19, col = rgb(0,0,0,0.7), bty = 'l', ylab = "% diff expr genes", main = "", xaxt = 'n', yaxt = 'n', 
     xlim = c(0.8, 6.2), ylim = c(0.30, 0.7), xlab = "", cex = 1.25)
points(2, 1-pi0est(res_B_KS1$P.Value[which(rownames(res_B_KS1) %in% common_ids)], pi0.method = "bootstrap")$pi0, pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(3, 1-pi0_vector_KS2[1], pch = 19, col = rgb(0,0,0,0.7), cex = 1.25)
points(4, 1-pi0est(res_B_KS2$P.Value[which(rownames(res_B_KS2) %in% common_ids)], pi0.method = "bootstrap")$pi0, pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
points(5, 1-pi0_vector_RT[1], pch = 19, col = rgb(0,0,0,0.7), cex = 1.25)
points(6, 1-pi0est(res_B_RT$P.Value[which(rownames(res_B_RT) %in% common_ids)], pi0.method = "bootstrap")$pi0, pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
abline(v = c(2.5, 4.5), lty = "longdash", col = rgb(0,0,0,0.7))
axis(2, at = c(0.30, 0.50, 0.70), labels = c(30, 50, 70))
axis(1, at = c(1.5, 3.5, 5.5), labels = c("KS1", "KS2", "RT"), cex.axis = 0.8, las = 1)
legend <- legend("topleft", legend = c("genes w/ commonly diff \naccessible promoters", 
                                       "top 1000 genes w/ uniquely diff \naccessible promoters"), 
                 col = c(alpha("red", 0.75), rgb(0,0,0,0.7)), 
                 bty = 'n', pch = 19, cex = 0.57)
dev.off()






###########
pi0_vector_KS1 <- sapply(seq(1000, 5000, by = 250), function(xx) 
  pi0est(atac_B_KS1_granges2$pvalue[which(atac_B_KS1_granges2$prom_overlapping_id %in% rownames(res_B_KS1)[
    which(rank(res_B_KS1$P.Value) < xx)])], pi0.method = "bootstrap")$pi0)

pi0_vector_KS2 <- sapply(seq(1000, 5000, by = 250), function(xx) 
  pi0est(atac_B_KS2_granges2$pvalue[which(atac_B_KS2_granges2$prom_overlapping_id %in% rownames(res_B_KS2)[
    which(rank(res_B_KS2$P.Value) < xx)])], pi0.method = "bootstrap")$pi0)

pi0_vector_RT <- sapply(seq(1000, 5000, by = 250), function(xx) 
  pi0est(atac_B_RT_granges2$pvalue[which(atac_B_RT_granges2$prom_overlapping_id %in% rownames(res_B_RT)[
    which(rank(res_B_RT$P.Value) < xx)])], pi0.method = "bootstrap")$pi0)

quartz(file = "accessibility_expression_comparison_3.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4, 4, 1, 1)+0.1)
plot(1:17, 1-pi0_vector_KS1, pch = 19, col = "orange", bty = 'l', ylab = "% diff accessible prom. peaks", main = "", xaxt = 'n', yaxt = 'n', 
     xlim = c(0.8, 54), ylim = c(0.15, 0.52), xlab = "diff expression ranking", cex = 0.84)
points(20:36, 1-pi0_vector_KS2, pch = 19, col = "forest green", cex = 0.84)
points(39:55, 1-pi0_vector_RT, pch = 19, col = "deep pink", cex = 0.84)
#abline(v = c(17.5, 34.5), lty = "longdash", col = rgb(0,0,0,0.7))
axis(2, at = c(0.20, 0.35, 0.50), labels = c(20, 35, 50))
text(1, 1-pi0_vector_KS1[1], "top\n1000", cex = 0.84, pos = 4)
text(17, 1-pi0_vector_KS1[17], "top\n5000", cex = 0.84, pos = 1)
text(20, 1-pi0_vector_KS2[1], "top\n1000", cex = 0.84, pos = 4)
text(36, 1-pi0_vector_KS2[17], "top\n5000", cex = 0.84, pos = 4)
text(39, 1-pi0_vector_RT[1], "top\n1000", cex = 0.84, pos = 2)
text(54, 1-pi0_vector_RT[17], "top\n5000", cex = 0.84, pos = 1)
legend <- legend("topleft", legend = c("KS1", "KS2", "RT"), pch = 19, bty = 'n', cex = 0.57, 
                 col = c("orange", "forest green", "deep pink"))
dev.off()

