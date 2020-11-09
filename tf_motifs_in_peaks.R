library(BSgenome.Mmusculus.UCSC.mm10)
load(file = "/dcl01/hansen/data/kabuki_mice/HiSeq325_ATAC_new/ATAC/dedup_files/merged_files/combined_granges_from_peaks_called_from_condition_specific_bams_with_all_fragments_B.rda")

setwd('/dcl01/hansen/data/archetypal_tf_motifs_mm10')
library(readr)
df <- read_tsv('mm10.archetype_motifs.v1.0.bed', col_names = FALSE)


colnames(df)[c(1,2,3)] <- c("chr", "start", "end")
df <- df[which(df$chr %in% seqnames(Mmusculus)[1:21]), ]

df_granges <- makeGRangesFromDataFrame(df, keep.extra.columns = FALSE, starts.in.df.are.0based = TRUE)
overlap_indices <- unique(queryHits(findOverlaps(df_granges, combined_granges_from_peaks_called_from_all_fragments_B)))
motifs_df <- df[overlap_indices, ]
colnames(motifs_df)[c(4:8)] <- c("motif_cluster_name", "match_score", "strand", "best_motif_match", "number_of_models")

save(motifs_df, file = "motifs_overlapping_B_cell_peaks_df.rda")





###tf motif analysis. Looks for TF motif
load(file = "motifs_overlapping_B_cell_peaks_df.rda")

getTFMotifEnrichmentDF <- function(tf_motif_df, atac_grange){
  output_df <- data.frame(motif_name = unique(tf_motif_df$motif_cluster_name))
  output_df$odds_ratio <- NA
  output_df$p_value <- NA

  for (i in 1:length(output_df$motif_name)){
    tf_df <- tf_motif_df[which(tf_motif_df$motif_cluster_name == output_df$motif_name[i]), ]
    g_range_tf_df <- makeGRangesFromDataFrame(tf_df, starts.in.df.are.0based = TRUE)
    atac_grange$contains_motif <- "no"
    overlaps <- findOverlaps(atac_grange, g_range_tf_df)
    indices <- unique(queryHits(overlaps))
    atac_grange$contains_motif[indices] <- "yes"
    n11 <- length(which(atac_grange$differential == "yes" & atac_grange$contains_motif == "yes"))
    n12 <- length(which(atac_grange$differential == "yes" & atac_grange$contains_motif == "no"))
    n21 <- length(which(atac_grange$differential == "no" & atac_grange$contains_motif == "yes"))
    n22 <- length(which(atac_grange$differential == "no" & atac_grange$contains_motif == "no"))
    fisher_output <- fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))
    output_df$odds_ratio[i] <- fisher_output$estimate
    output_df$p_value[i] <- fisher_output$p.value
    output_df
  }
  output_df$fdr <- p.adjust(output_df$p_value, method = "fdr")
  output_df[order(output_df$odds_ratio, decreasing = TRUE), ]
}

####enrichment in diff acc peaks at promoters of diff expr genes
atac_B_KS1_granges2$differential <- "no"
atac_B_KS1_granges2$differential[which(atac_B_KS1_granges2$prom_overlapping_id %in% 
                                         rownames(top_1000_proms_KS1_de_genes))] <- "yes"

atac_B_KS2_granges2$differential <- "no"
atac_B_KS2_granges2$differential[which(atac_B_KS2_granges2$prom_overlapping_id %in% 
                                         rownames(top_1000_proms_KS2_de_genes))] <- "yes"

atac_B_RT_granges2$differential <- "no"
atac_B_RT_granges2$differential[which(atac_B_RT_granges2$prom_overlapping_id %in% 
                                        rownames(top_1000_proms_RT_de_genes))] <- "yes"



motif_enrichment_df_KS1_DE_proms <- getTFMotifEnrichmentDF(motifs_df, atac_B_KS1_granges2)
motif_enrichment_df_KS2_DE_proms <- getTFMotifEnrichmentDF(motifs_df, atac_B_KS2_granges2)
motif_enrichment_df_RT_DE_proms <- getTFMotifEnrichmentDF(motifs_df, atac_B_RT_granges2)


####enrichment in peaks of diff expr genes regardless of whether the peaks are diff accessible
atac_B_KS1_granges2$differential <- "no"
atac_B_KS1_granges2$differential[which(atac_B_KS1_granges2$prom_overlapping_id %in% 
                                         rownames(res_B_KS1)[which(res_B_KS1$padj < 0.1)])] <- "yes"

atac_B_KS2_granges2$differential <- "no"
atac_B_KS2_granges2$differential[which(atac_B_KS2_granges2$prom_overlapping_id %in% 
                                         rownames(res_B_KS2)[which(res_B_KS2$padj < 0.1)])] <- "yes"

atac_B_RT_granges2$differential <- "no"
atac_B_RT_granges2$differential[which(atac_B_RT_granges2$prom_overlapping_id %in% 
                                        rownames(res_B_RT)[which(res_B_RT$padj < 0.1)])] <- "yes"


motif_enrichment_df_KS1_de_proms_2 <- getTFMotifEnrichmentDF(motifs_df, atac_B_KS1_granges2)
motif_enrichment_df_KS2_de_proms_2 <- getTFMotifEnrichmentDF(motifs_df, atac_B_KS2_granges2)
motif_enrichment_df_RT_de_proms_2 <- getTFMotifEnrichmentDF(motifs_df, atac_B_RT_granges2)















