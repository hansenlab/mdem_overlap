###
###
###compare our peaks with publicly available encode peaks

encode_dnase_B_cells_rep1 <- getGrangesFromNarrowPeaks('B_cells_DNase_seq_mm9_Rep1.narrowPeak')
write_tsv(data.frame(chr = seqnames(encode_dnase_B_cells_rep1), 
                     start = start(encode_dnase_B_cells_rep1), 
                     end = end(encode_dnase_B_cells_rep1)), "B_cells_DNase_seq_mm9_Rep1.txt")


write_csv(data.frame(chrom = paste0(seqnames(encode_dnase_B_cells_rep1), ":",
                                    start = start(encode_dnase_B_cells_rep1), "-",
                                    end = end(encode_dnase_B_cells_rep1))), "B_cells_DNase_seq_mm9_Rep1.bed")

encode_dnase_B_cells_rep2 <- getGrangesFromNarrowPeaks('B_cells_DNase_seq_mm9_Rep2.narrowPeak')
write_tsv(data.frame(chr = seqnames(encode_dnase_B_cells_rep2), 
                     start = start(encode_dnase_B_cells_rep2), 
                     end = end(encode_dnase_B_cells_rep2)), "B_cells_DNase_seq_mm9_Rep2.txt")


write_csv(data.frame(coord =  paste0(seqnames(encode_dnase_B_cells_rep2), ":",
                                     start = start(encode_dnase_B_cells_rep2), "-", 
                                     end = end(encode_dnase_B_cells_rep2))), "B_cells_DNase_seq_mm9_Rep2.bed")


##covert to mm10 using liftover and then import again
encode_dnase_B_cells_rep1 <- read_tsv('B_cells_DNase_seq_mm10_Rep1.txt', col_names = FALSE)
colnames(encode_dnase_B_cells_rep1) <- "coord"
encode_dnase_B_cells_rep1$coord <- gsub(":", "_", encode_dnase_B_cells_rep1$coord)
encode_dnase_B_cells_rep1$coord <- gsub("-", "_", encode_dnase_B_cells_rep1$coord)
column_split <- strsplit(encode_dnase_B_cells_rep1$coord, "_")
encode_dnase_B_cells_rep1$chr <- unlist(lapply(column_split, function(xx) xx[1]))
encode_dnase_B_cells_rep1$start <- unlist(lapply(column_split, function(xx) xx[2]))
encode_dnase_B_cells_rep1$end <- unlist(lapply(column_split, function(xx) xx[3]))
encode_dnase_B_cells_rep1 <- encode_dnase_B_cells_rep1[, -1]
encode_dnase_B_cells_rep1 <- makeGRangesFromDataFrame(encode_dnase_B_cells_rep1)


encode_dnase_B_cells_rep2 <- read_tsv('B_cells_DNase_seq_mm10_Rep2.txt', col_names = FALSE)
colnames(encode_dnase_B_cells_rep2) <- "coord"
encode_dnase_B_cells_rep2$coord <- gsub(":", "_", encode_dnase_B_cells_rep2$coord)
encode_dnase_B_cells_rep2$coord <- gsub("-", "_", encode_dnase_B_cells_rep2$coord)
column_split <- strsplit(encode_dnase_B_cells_rep2$coord, "_")
encode_dnase_B_cells_rep2$chr <- unlist(lapply(column_split, function(xx) xx[1]))
encode_dnase_B_cells_rep2$start <- unlist(lapply(column_split, function(xx) xx[2]))
encode_dnase_B_cells_rep2$end <- unlist(lapply(column_split, function(xx) xx[3]))
encode_dnase_B_cells_rep2 <- encode_dnase_B_cells_rep2[, -1]
encode_dnase_B_cells_rep2 <- makeGRangesFromDataFrame(encode_dnase_B_cells_rep2)


combined_B_cell_encode_peaks <- reduce(unlist(GRangesList(encode_dnase_B_cells_rep1, encode_dnase_B_cells_rep2)))



