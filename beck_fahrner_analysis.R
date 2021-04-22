meth_signature <- read_csv('tet3_signature.csv')
load(file = "expression_properties_df_gene_level.rda") #contains gtex data

library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)
edb <- EnsDb.Hsapiens.v75
proms_all <- promoters(edb, upstream = 2000, downstream = 2000, 
                       columns = c("gene_name", "tx_id", "tx_cds_seq_start", "tx_cds_seq_end",
                                   "tx_biotype", "gene_id"))
names(proms_all) <- proms_all$tx_id
genome(seqinfo(proms_all)) <- "hg19"
seqlevelsStyle(proms_all) <- "ucsc"

chrs <- names(Hsapiens)[1:22]
proms_all <- proms_all[which(seqnames(proms_all) %in% chrs[1:22])] #exclude promoters on the sex chromosomes because it is hard to talk about LoF-intolerance for those
seqlevels(proms_all) <- seqlevels(proms_all)[1:22]

meth_signature$gene_ids <- sapply(meth_signature$`Overlapping Genes`, function(xx) 
  unique(proms_all$gene_id[which(proms_all$gene_name %in% unlist(strsplit(xx, ",")))]))

meth_signature$brain_expr <- sapply(meth_signature$gene_ids, function(xx) 
  expr_properties_df$median_expr_in_brain[which(rownames(expr_properties_df) %in% xx)])
  
meth_signature$brain_expr_max <- sapply(meth_signature$brain_expr, function(xx) xx[which.max(xx)])
meth_signature$brain_expr_max <- as.numeric(meth_signature$brain_expr_max) 

#all_canonical_tx_ids <- all_canonical_tx_ids[which(all_canonical_tx_ids %in% names(proms_all))]
all_canonical_gene_ids <- unique(proms_all$gene_id[which(proms_all$tx_biotype == "protein_coding")])
all_canonical_gene_ids <- all_canonical_gene_ids[-which(all_canonical_gene_ids %in% unlist(meth_signature$gene_ids))]
all_canonical_brain_expr <- expr_properties_df[all_canonical_gene_ids, "median_expr_in_brain"]



####
fetal_expr <- read_csv('fetal_single_cell_expression/gene_expression_celltype.txt') 
fetal_expr <- read_csv('fetal_single_cell_expression/gene_fraction_celltype.txt') #if I want proportion of cells with >0 UMI

fetal_expr$RowID <- gsub("[.].*", "", fetal_expr$RowID)
fetal_expr_gene_ids <- fetal_expr$RowID
fetal_expr <- as.matrix(fetal_expr[, -1])
rownames(fetal_expr) <- fetal_expr_gene_ids

#
meth_signature$cerebrum_excitatory_neurons_expr <- sapply(meth_signature$gene_ids, function(xx) 
  fetal_expr[which(rownames(fetal_expr) %in% xx), "Cerebrum-Excitatory neurons"])
meth_signature$cerebrum_excitatory_neurons_expr_max <- sapply(meth_signature$cerebrum_excitatory_neurons_expr, 
                                                              function(xx) xx[which.max(xx)])
meth_signature$cerebrum_excitatory_neurons_expr_max <- as.numeric(meth_signature$cerebrum_excitatory_neurons_expr_max)

meth_signature$cerebrum_inhibitory_neurons_expr <- sapply(meth_signature$gene_ids, function(xx) 
  fetal_expr[which(rownames(fetal_expr) %in% xx), "Cerebrum-Inhibitory neurons"])
meth_signature$cerebrum_inhibitory_neurons_expr_max <- sapply(meth_signature$cerebrum_inhibitory_neurons_expr, 
                                                              function(xx) xx[which.max(xx)])
meth_signature$cerebrum_inhibitory_neurons_expr_max <- as.numeric(meth_signature$cerebrum_inhibitory_neurons_expr_max)
  
all_canonical_excitatory <- fetal_expr[all_canonical_gene_ids, "Cerebrum-Excitatory neurons"]
all_canonical_inhibitory <- fetal_expr[all_canonical_gene_ids, "Cerebrum-Inhibitory neurons"]


###
plot_df <- data.frame(cell_type = c(rep("fetal cerebral\nexcitatory neurons", 19286), 
                                    rep("fetal cerebral\ninhibitory neurons", 19286)), 
                     expr = c(meth_signature$cerebrum_excitatory_neurons_expr_max, 
                              all_canonical_excitatory, 
                              meth_signature$cerebrum_inhibitory_neurons_expr_max, 
                              all_canonical_inhibitory), 
                     group = c(rep("DMR\ngenes", 50), rep("other", 19236), 
                               rep("DMR\ngenes", 50), rep("other", 19236)))

quartz(file = "signature_genes_expression_comparison.pdf", height = 3, width = 3, pointsize = 8, type = "pdf")
ggplot(plot_df, aes(x = group, y = expr)) + 
  geom_boxplot(fill = c(alpha("red", 0.59), "cornflowerblue", alpha("red", 0.59), "cornflowerblue"), 
               color = rep(rgb(0,0,0,0.7), 4), outlier.shape=NA) + 
  coord_cartesian(ylim = c(0, 200)) + labs(y = "expression (TPM)", x = "") + 
  theme(axis.text=element_text(size=1)) + theme_classic() +
  facet_wrap(~cell_type) 
dev.off()

wilcox.test(meth_signature$cerebrum_excitatory_neurons_expr_max, all_canonical_excitatory, alt = "g")
wilcox.test(meth_signature$cerebrum_inhibitory_neurons_expr_max, all_canonical_inhibitory, alt = "g")
  
