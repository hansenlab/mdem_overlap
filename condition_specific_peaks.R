
###
library(Rsubread)
getGrangesFromNarrowPeaks <- function(narrowPeak_filepath){
  peaks <- read.delim(paste0(narrowPeak_filepath), header = FALSE)
  colnames(peaks) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "fold.enrichment",
                       "log10.pval", "log10.qval", "peak")
  peaks_granges <- GRanges(seqnames = peaks$chrom, IRanges(peaks$chromStart+1, peaks$chromEnd), 
                           fold.enrichment = peaks$fold.enrichment, log10.pval = peaks$log10.pval, 
                           log10.qval = peaks$log10.qval)
  peaks_granges
}

#
B_wt_granges <- getGrangesFromNarrowPeaks('all_wt_samples_B_peaks.narrowPeak')
B_KS1_granges <- getGrangesFromNarrowPeaks('all_KS1_samples_B_peaks.narrowPeak')
B_KS2_granges <- getGrangesFromNarrowPeaks('all_KS2_samples_B_peaks.narrowPeak')

B_wt_granges_cbp_cohort <- getGrangesFromNarrowPeaks('/dcl01/hansen/data/kabuki_mice/HiSeq360/ATAC/all_wt_samples_B_rt_cohort_peaks.narrowPeak')
B_RT_granges_cbp_cohort <- getGrangesFromNarrowPeaks('/dcl01/hansen/data/kabuki_mice/HiSeq360/ATAC/all_RT_samples_B_rt_cohort_peaks.narrowPeak')

combined_granges_from_peaks_called_from_all_fragments_B <- reduce(unlist(GRangesList(B_wt_granges, B_KS1_granges, B_KS2_granges, 
                                                                                     B_wt_granges_cbp_cohort, 
                                                                                     B_RT_granges_cbp_cohort)))


combined_granges_from_peaks_called_from_all_fragments_B <- combined_granges_from_peaks_called_from_all_fragments_B[-unique(
  queryHits(findOverlaps(combined_granges_from_peaks_called_from_all_fragments_B, mm10_regions_to_exclude_granges)))]

save(combined_granges_from_peaks_called_from_all_fragments_B, 
     file = "combined_granges_from_peaks_called_from_condition_specific_bams_with_all_fragments_B.rda")

#######counts reads overlapping peaks in each of the bam files
combined_granges <- combined_granges_from_peaks_called_from_all_fragments_B

combined_granges$id <- paste0(seqnames(combined_granges), "_", 
                              as.character(start(combined_granges)), "_", as.character(end(combined_granges)))
#combined_peaks_ann <- createAnnotationFile(combined_granges) ###this function is longer supported so we run the following

combined_peaks_ann <- data.frame(GeneID = combined_granges$id, 
                                 Chr = seqnames(combined_granges), 
                                 Start = start(combined_granges), 
                                 End = end(combined_granges), 
                                 strand = strand(combined_granges), 
                                 stringsAsFactors = FALSE)


#Kabuki cohort
bamFiles_filepaths <- list.files(pattern = "noMito_merged.bam")
bamFiles_filepaths <- bamFiles_filepaths[-grep("bam.bai", bamFiles_filepaths)]

count_reads_in_peaks <- featureCounts(bamFiles_filepaths, annot.ext = combined_peaks_ann, 
                                      isPairedEnd = TRUE, 
                                      requireBothEndsMapped = TRUE, countChimericFragments = FALSE, 
                                      nthreads = 30, countMultiMappingReads = FALSE)

save(count_reads_in_peaks, file = "count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_kabuki_cohort.rda")


#RT cohort
setwd('/dcl01/hansen/data/kabuki_mice/HiSeq360/ATAC')
bamFiles_filepaths <- list.files(pattern = "noMito_merged.bam")

count_reads_in_peaks <- featureCounts(bamFiles_filepaths, annot.ext = combined_peaks_ann, 
                                      isPairedEnd = TRUE, 
                                      requireBothEndsMapped = TRUE, countChimericFragments = FALSE, 
                                      nthreads = 30, countMultiMappingReads = FALSE)

save(count_reads_in_peaks, file = "count_reads_in_combined_peaks_from_condition_specific_bam_files_B_cells_only_rt_cohort.rda")




