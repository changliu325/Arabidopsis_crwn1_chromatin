
#Combine ATAC-seq peaks identified from different genotypes and
#generate a master list of ATAC-seq peaks

options(stringsAsFactors = F, scipen=999)
setwd("~/Documents/ZMBP/DNA_methylation_and_NP/ATAC_seq/ATAC_peaks_2018_Aug/")
WT_2C <- read.table("ATAC_seq_SE_WT_2C_peaks.xls",header=T)
crwn1_2C <- read.table("ATAC_seq_SE_crwn1_2C_peaks.xls",header=T)
kaku4_2C <- read.table("ATAC_seq_SE_kaku4_2_2C_peaks.xls",header=T)

#merge peaks in different datasets
suppressPackageStartupMessages(library(IRanges))
peaks_vec <- vector(mode="list",length=5)
for(i in 1:5){
  peaks_set1 <- IRanges(start=WT_2C[WT_2C[,1]==i, "start"], end=WT_2C[WT_2C[,1]==i, "end"])
  peaks_set2 <- IRanges(start=crwn1_2C[crwn1_2C[,1]==i, "start"], end=crwn1_2C[crwn1_2C[,1]==i, "end"])
  peaks_set3 <- IRanges(start=kaku4_2C[kaku4_2C[,1]==i, "start"], end=kaku4_2C[kaku4_2C[,1]==i, "end"])
  
  peaks_vec[[i]] <-reduce(c(peaks_set1, peaks_set2, peaks_set3))
}
#creat a file describing merged peaks:
merged_peaks <- data.frame()
for(i in 1:5){
  merged_peaks_per_chr <- data.frame(chr=i, start=peaks_vec[[i]]@start, end=peaks_vec[[i]]@start+peaks_vec[[i]]@width-1,
                                     peak_name=paste("chr",i, 1:length(peaks_vec[[i]]), sep="_"))
  merged_peaks <- rbind(merged_peaks, merged_peaks_per_chr)
}

write.table(merged_peaks, file="merged_ATAC_peaks_2C_wt_crwn1_kaku4.txt", quote=F, row.names=F, col.names=F, sep="\t")
