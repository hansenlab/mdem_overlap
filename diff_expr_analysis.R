##KS1 and KS2

####read quantification files
files <- paste0("KS1_KS2_analysis/analysis_objects/quants/", 
                list.files("KS1_KS2_analysis/analysis_objects/quants"), "/quant.sf")

names(files) <- c("100.B.KDM6A_1", "100.B.KDM6A_2", "100.B.KDM6A_3", 
                  "100.T.KDM6A_1", "100.T.KDM6A_2", "100.T.KDM6A_3", 
                  "101.B.WT_1", "101.B.WT_2", "101.B.WT_3", 
                  "101.T.WT_1", "101.T.WT_2", "101.T.WT_3",
                  "102.B.KDM6A_1", "102.B.KDM6A_2", "102.B.KDM6A_3", 
                  "102.T.KDM6A_1", "102.T.KDM6A_2", "102.T.KDM6A_3",
                  "103.B.KDM6A_1", "103.B.KDM6A_2", "103.B.KDM6A_3", 
                  "103.T.KDM6A_1", "103.T.KDM6A_2", "103.T.KDM6A_3",
                  "104.B.WT_1", "104.B.WT_2", "104.B.WT_3", 
                  "104.T.WT_1", "104.T.WT_2", "104.T.WT_3", 
                  "106.B.KDM6A_1", "106.B.KDM6A_2", "106.B.KDM6A_3", 
                  "106.T.KDM6A_1", "106.T.KDM6A_2", "106.T.KDM6A_3",
                  "23.B.KMT2D_1", "23.B.KMT2D_2", "23.B.KMT2D_3", 
                  "23.T.KMT2D_1", "23.T.KMT2D_2", "23.T.KMT2D_3",
                  "27.B.KMT2D_1", "27.B.KMT2D_2", "27.B.KMT2D_3", 
                  "27.T.KMT2D_1", "27.T.KMT2D_2", "27.T.KMT2D_3",
                  "32.B.WT_1", "32.B.WT.2", "32.B.WT_3", 
                  "32.T.WT_1", "32.T.WT_2", "32.T.WT_3",
                  "40.B.WT_1", "40.B.WT_2", "40.B.WT_3", 
                  "40.T.WT_1", "40.T.WT_2", "40.T.WT_3",
                  "42.B.KDM6A_1", "42.B.KDM6A_2", "42.B.KDM6A_3", 
                  "42.T.KDM6A_1", "42.T.KDM6A_2", "42.T.KDM6A_3", 
                  "63.B.WT_1", "63.B.WT_2", "63.B.WT_3", 
                  "63.T.WT_1", "63.T.WT_2", "63.T.WT_3",
                  "6.B.KMT2D_1", "6.B.KMT2D_2", "6.B.KMT2D_3", 
                  "6.T.KMT2D_1", "6.T.KMT2D_2", "6.T.KMT2D_3",
                  "79.B.KMT2D_1", "79.B.KMT2D_2", "79.B.KMT2D_3", 
                  "79.T.KMT2D_1", "79.T.KMT2D_2", "79.T.KMT2D_3",
                  "84.B.KMT2D_1", "84.B.KMT2D_2", "84.B.KMT2D_3", 
                  "84.T.KMT2D_1", "84.T.KMT2D_2", "84.T.KMT2D_3")


######
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1] 

#get expression matrix
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

counts_mat <- txi$counts
combined_mat_cols <- lapply(seq(1, 90, by = 3), function(xx) 
  rowSums(counts_mat[, c(xx, xx+1, xx+2), drop = FALSE]))

combined_mat <- matrix(unlist(combined_mat_cols), ncol = 30)
rownames(combined_mat) <- names(combined_mat_cols[[1]])
colnames(combined_mat) <- gsub("_.*", "", colnames(counts_mat)[seq(1, 90, by = 3)])

colnames(combined_mat) <- sub("[.]", "_", colnames(combined_mat))
colnames(combined_mat) <- sub("[.]", "_", colnames(combined_mat)) #for some reason I need to run this twice. there must be something wrong



####run differential analysis
genotype <- "KS1" #analogously for KS2
celltype <- "B"
features_by_samples_mat <- combined_mat

sample_info <- data.frame(genotype = getGenotypeVector(colnames(features_by_samples_mat)), 
                          cell_type = getCellTypeVector(colnames(features_by_samples_mat)))

features_by_samples_mat <- features_by_samples_mat[, which(sample_info$genotype %in% c("WT", genotype) & 
                                                             sample_info$cell_type == celltype)]

sample_info2 <- sample_info[which(sample_info$genotype %in% c("WT", genotype) & 
                                    sample_info$cell_type == celltype), , drop = FALSE]

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
res_B_KS1 <- res



###RT
files <- paste0("rt_rnaseq/quants/", 
                list.files("rt_rnaseq/quants"), "/quant.sf")

names(files) <- c("8.B.CBP_1", "8.B.CBP_2", "8.B.CBP_3",
                  "8.T.CBP_1", "8.T.CBP_2", "8.T.CBP_3",
                  "5.B.WT_1", "5.B.WT_2", "5.B.WT_3",
                  "5.T.WT_1","5.T.WT_2", "5.T.WT_3",
                  "16.B.WT_1", "16.B.WT_2", "16.B.WT_3",
                  "16.T.WT_1", "16.T.WT_2", "16.T.WT_3",
                  "15.B.WT_1", "15.B.WT_2", "15.B.WT_3",
                  "15.T.WT_1", "15.T.WT_2", "15.T.WT_3",
                  "21.B.CBP_1", "21.B.CBP_2", "21.B.CBP_3", 
                  "21.T.CBP_1", "21.T.CBP_2", "21.T.CBP_3", 
                  "6.B.WT_1", "6.B.WT_2", "6.B.WT_3",
                  "6.T.WT_1", "6.T.WT_2", "6.T.WT_3",
                  "10.B.CBP_1", "10.B.CBP_2", "10.B.CBP_3", 
                  "10.T.CBP_1", "10.T.CBP_2", "10.T.CBP_3", 
                  "66.B.CBP_1", "66.B.CBP_2", "66.B.CBP_3",
                  "66.T.CBP_1", "66.T.CBP_2", "66.T.CBP_3",
                  "44.B.WT_1", "44.B.WT_2", "44.B.WT_3",
                  "44.T.WT_1", "44.T.WT_2", "44.T.WT_3",
                  "2.B.WT_1", "2.B.WT_2", "2.B.WT_3",
                  "2.T.WT_1", "2.T.WT_2", "2.T.WT_3",
                  ".B.WT_1", ".B.WT_2", ".B.WT_3",
                  ".T.WT_1",  ".T.WT_2", ".T.WT_3",
                  "42.B.CBP_1", "42.B.CBP_2", "42.B.CBP_3", 
                  "42.T.CBP_1", "42.T.CBP_2", "42.T.CBP_3")


#and the rest of the steps as for KS1/2



