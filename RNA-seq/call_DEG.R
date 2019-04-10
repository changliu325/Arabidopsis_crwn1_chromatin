options(stringsAsFactors=F)
options(width=200)
library(Rsamtools)
library(GenomicAlignments)

##############################################################
gtf <- read.delim("/ebio/abt6_projects7/small_projects/wzhu/bin/TAIR10_index/Arabidopsis_thaliana.TAIR10.24_tophat.gtf",
                  header=FALSE,stringsAsFactors=FALSE)
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame","attributes")
chronly <- 1:5    #chronly: which chromosome to extract?
gtf <- gtf[as.character(gtf$seqname) %in% chronly, ] # Cleanup to remove non-chromosome rows
gtf <- gtf[gtf$feature =="exon",] #only exon 
gene.ids <-  gsub(".*gene_id (.*?);.*", "\\1", gtf$attributes) #get gene_id from attributes column
#gene.name <- gsub(".*gene_name (.*?);.*", "\\1", gtf$attributes) #get gene_name from attributes column
#transcript.ids<- gsub(".*transcript_id (.*?);.*", "\\1", gtf$attributes) #get transcript_id from attributes column
gene.index<-gene.ids!=""             #skip those have no value
#transcript.index<-transcript.ids!=""   #skip those have no value
gene.gr<-GRanges(seqnames=gtf$seqname[gene.index],
                   ranges=IRanges(gtf$start[gene.index],gtf$end[gene.index]),
                   #strand="*",
                   strand=gtf$strand[gene.index],
                   #tx_id=transcript.ids[gene.index],
                   gene_id=gene.ids[gene.index]
                   #gene_name=gene.name[gene.index]
                   )
gene.gr.list <- split(gene.gr,gene.ids[gene.index]) # feature list

bam_files <- c("/ebio/abt6_projects7/small_projects/wzhu/data/chang/col_rep1_thout/accepted_hits.bam",
               "/ebio/abt6_projects7/small_projects/wzhu/data/chang/col_rep2_thout/accepted_hits.bam",
               "/ebio/abt6_projects7/small_projects/wzhu/data/chang/col_rep3_thout/accepted_hits.bam",
               "/ebio/abt6_projects7/small_projects/wzhu/data/chang/crwn_rep1_thout/accepted_hits.bam",
               "/ebio/abt6_projects7/small_projects/wzhu/data/chang/crwn_rep2_thout/accepted_hits.bam",
               "/ebio/abt6_projects7/small_projects/wzhu/data/chang/crwn_rep3_thout/accepted_hits.bam")

names(bam_files) <- c("Col_rep1","Col_rep2", "Col_rep3", "crwn_rep1","crwn_rep2", "crwn_rep3")

reads_per_gene <- list()
for(i in names(bam_files)){
   reads_per_gene[[i]] <- summarizeOverlaps(gene.gr.list, readGAlignments(bam_files[i], use.names=TRUE))
   print(paste(i,"--done"))
}

reads_count_table <- data.frame(Col_rep1=assay(reads_per_gene[["Col_rep1"]]),
                                Col_rep2=assay(reads_per_gene[["Col_rep2"]]),
                                Col_rep3=assay(reads_per_gene[["Col_rep3"]]),
                                crwn_rep1=assay(reads_per_gene[["crwn_rep1"]]),
                                crwn_rep2=assay(reads_per_gene[["crwn_rep2"]]),
                                crwn_rep3=assay(reads_per_gene[["crwn_rep3"]])
                                )

colnames(reads_count_table) <- names(bam_files)


x <- suppressWarnings(unlist(seqnames(gene.gr.list)))
chrs <- as.numeric(as.character(x[gene.gr.list@partitioning@end]))
reads_count_table <- data.frame(chr=chrs,reads_count_table)
### this table would be used for DESeq2
############################################################

library(DESeq2)

#crwn vs Col

countData <- data.matrix(reads_count_table[reads_count_table[,"chr"]%in%c(1:5),c(2,3,4,5,6,7)])
colData <- data.frame(genotype=factor(c("Col","Col", "Col", "crwn", "crwn", "crwn")))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ genotype)
dds$genotype <- relevel(dds$genotype, "Col")
dds_crwn_vs_Col <- DESeq(dds)
res_crwn_vs_Col <- results(dds_crwn_vs_Col,independentFiltering=F)
