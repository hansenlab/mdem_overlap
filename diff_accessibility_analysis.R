#####ATAC analysis 
peaks_by_samples <- getPeaksBySamplesMatrix("rt_atacseq/count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_kabuki_cohort.rda", 
                                                   "blood_all_samples")
#or for RT
peaks_by_samples <- getPeaksBySamplesMatrix("rt_atacseq/count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_rt_cohort.rda", 
                                                   "rt_blood")
sample_info <- data.frame(genotype = getGenotypeVector(colnames(peaks_by_samples)), 
                          cell_type = getCellTypeVector(colnames(peaks_by_samples)))


genotype <- "RT1" #or KS2, or RT
celltype <- "B"  #or B cells

features_by_samples_mat <- peaks_by_samples[, which(sample_info$genotype %in% c("WT", genotype) & 
                                                             sample_info$cell_type == celltype)]

sample_info2 <- sample_info[which(sample_info$genotype %in% c("WT", genotype) & 
                                    sample_info$cell_type == celltype), , drop = FALSE]
sample_info2$genotype <- droplevels(sample_info2$genotype) #this is not needed for RT

sampleTable <- data.frame(condition = factor(sample_info2$genotype))
dds <- DESeqDataSetFromMatrix(round(features_by_samples_mat), sampleTable, design = ~condition)

idx <- rowMedians(counts(dds)) > 10
dat <- counts(dds)[idx,]
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svseq <- svaseq(dat, mod, mod0)

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
design(ddssva) <- formula(~ SV1 + SV2 + SV3 + condition)
ddssva <- DESeq(ddssva)

dds <- ddssva
res <- results(dds)
if (length(which(is.na(res$padj))) > 0){res <- res[-which(is.na(res$padj)), ]}
res$P.Value <- res$pvalue
res_B_KS1_atac <- res



