##make sure that the table about enzyme cutting site has been prepared
args <- commandArgs(TRUE)

mapped_reads <- args[1] #input file, with mapped HiC reads, 8 columns, no header
#format: chr_a; pos_a; direction_a; chr_b; pos_b; direction_b,
#  NOTE!!!!  where directions are "0" or "16", meaning + or - strand
cutting_table <- args[2] #input file, with 4 columns, with header
#format: chr; start; end; enzyme_id

####################
#
#  First filter, this filter checks both intra- and inter-chromosomal reads
#
####################
cutoff <- 700 # the sum of distances shall be smaller than it, this is subject to change depending on 
# how library was constructed. e.g. if select sonicated DNA as 300-500bp, the cutoff shall be 500


print(paste("Loading----",mapped_reads))
options(stringsAsFactors = F)
options(scipen = 999)
reads <- read.table(file=mapped_reads)
cutting_sites <- read.table(file=cutting_table, header=T)
print(paste("Loading----",mapped_reads, "---completed"))

if(nrow(reads)==0){ # this could happen if above filter is turned on
  filename <- paste(mapped_reads, ".result", sep="")
  write.table(reads, file=filename, row.names=F, col.names=F, quote=F, sep="\t")
  print(paste("filtering----",mapped_reads, "----completed"))
  q("no")
}

#######################################


reads <- cbind(reads, distance_a=cutoff, distance_b=cutoff, flag=0) 
# distance_a/b, distance to next cutting site; flag: 0 becomes 1 if a read pair passes filter
#seperate the job into n parts (n is the chromosome number)
n <- table(cutting_sites[,1])

suppressPackageStartupMessages(library(IRanges))
position_vec_plus <- list() 
position_vec_minus <- list()
for(i in names(n)){
              temp <- cutting_sites[cutting_sites[,1]==i,]
              vector_len <- nrow(temp)
              position_vec_plus[[i]] <- IRanges(start=temp[,2][1:(vector_len-1)], 
                                           width=temp[,2][-1]-temp[,2][1:(vector_len-1)] )
              #(for reads mapped to + strand) 
              
              
            position_vec_minus[[i]] <- IRanges(start=temp[,3][1:(vector_len-1)], 
                                           width=temp[,3][-1]-temp[,3][1:(vector_len-1)])
              #(for reads mapped to - strand)
              
              
}

reads_split <- split(reads[,1:4], reads[,1]) # analyse read_a


suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

registerDoParallel(cores=5)
distance_a <- foreach(i=names(reads_split)) %dopar% {
  distance_chr <- rep(cutoff, nrow(reads_split[[i]])) #distance to be calculated
  
  fragments_plus <- findOverlaps(IRanges(reads_split[[i]][reads_split[[i]][,3]==0,2]+reads_split[[i]][reads_split[[i]][,3]==0,4]-20, width=1), 
                                         position_vec_plus[[i]], select="first")
  fragments_minus <- findOverlaps(IRanges(reads_split[[i]][reads_split[[i]][,3]==16,2]-reads_split[[i]][reads_split[[i]][,3]==16,4]+20, width=1), 
                                          position_vec_minus[[i]], select="first")
  #these two lines give you the fragment id where checking position belongs to
  #it searches the 15th last nt of each read
  
  site_plus <- position_vec_plus[[i]][fragments_plus]@"start"+position_vec_plus[[i]][fragments_plus]@"width"
  site_minus <- position_vec_minus[[i]][fragments_minus]@"start"
  #these two lines give you the location of the first cutting site downstream of mapping direction
  
  distance_plus <- site_plus-reads_split[[i]][reads_split[[i]][,3]==0,2]  #The distance between 5' of reads and the next possible cutting site
  distance_minus <- reads_split[[i]][reads_split[[i]][,3]==16,2]-site_minus
  #so here are the distances
  
  distance_chr[reads_split[[i]][,3]==0] <- distance_plus
  distance_chr[reads_split[[i]][,3]==16] <- distance_minus
  
  return(distance_chr)
}
names(distance_a) <- names(reads_split) #names of each member in the list corresponds to chromosome number
#############################################
reads_split <- split(reads[,5:8], reads[,5]) # analyse read_b
distance_b <- foreach(i=names(reads_split)) %dopar% {
  distance_chr <- rep(cutoff, nrow(reads_split[[i]])) #distance to be calculated
  
  fragments_plus <- findOverlaps(IRanges(reads_split[[i]][reads_split[[i]][,3]==0,2]+reads_split[[i]][reads_split[[i]][,3]==0,4]-20, width=1), 
                                         position_vec_plus[[i]], select="first")
  fragments_minus <- findOverlaps(IRanges(reads_split[[i]][reads_split[[i]][,3]==16,2]-reads_split[[i]][reads_split[[i]][,3]==16,4]+20, width=1), 
                                          position_vec_minus[[i]], select="first")

  site_plus <- position_vec_plus[[i]][fragments_plus]@"start"+position_vec_plus[[i]][fragments_plus]@"width"
  site_minus <- position_vec_minus[[i]][fragments_minus]@"start"
  
  distance_plus <- site_plus-reads_split[[i]][reads_split[[i]][,3]==0,2]  #The distance between 5' of reads and the next possible cutting site
  distance_minus <- reads_split[[i]][reads_split[[i]][,3]==16,2]-site_minus
 
  distance_chr[reads_split[[i]][,3]==0] <- distance_plus
  distance_chr[reads_split[[i]][,3]==16] <- distance_minus
  
  return(distance_chr)
}
names(distance_b) <- names(reads_split)
registerDoParallel(cores=1)
### Now filling distance_a and distance_b:
for(i in names(distance_a)){
  reads[reads[,1]==i,"distance_a"] <- distance_a[[i]]
}
for(i in names(distance_b)){
  reads[reads[,5]==i,"distance_b"] <- distance_b[[i]]
}


rm(distance_a); rm(distance_b)
reads[abs(reads[,"distance_a"]+reads[,"distance_b"])<=cutoff,"flag"] <- 1

reads_filtered_inter <- reads[reads[,"flag"]==1&(reads[,1]!=reads[,5]),c(1:8)]
####################
#
#  Second and Third filter, these two filters further check intra-chromosomal reads
#
####################

reads_intra <- reads[reads[,"flag"]==1&(reads[,1]==reads[,5]),c(1:8,11)]
if(nrow(reads_intra)>0){
  
#Second filter:
###############
#remove dangling reads (as an intact genomic DNA fragment with size less than cutoff)
reads_intra[((reads_intra[,2]<reads_intra[,6]&reads_intra[,3]<reads_intra[,7])|
              (reads_intra[,2]>reads_intra[,6]&reads_intra[,3]>reads_intra[,7]))&
              (abs(reads_intra[,2]-reads_intra[,6])<cutoff),"flag"] <- 0 # they are potential self ligation products
reads_intra <- reads_intra[reads_intra[,"flag"]==1,] 

#Third filter:
##############
# for MP mode (back-to-back mapping direction) intra-chromosomal interaction, if the mapping distance is less than 5000, make sure 
# that there is at least one cutting site
reads_intra[((reads_intra[,2]<reads_intra[,6]&reads_intra[,3]>reads_intra[,7])|
              (reads_intra[,2]>reads_intra[,6]&reads_intra[,3]<reads_intra[,7]))&
              (abs(reads_intra[,2]-reads_intra[,6])<=5000),"flag"] <- 0 # they are potential self ligation products

chromosome_names <- names(table(reads_intra[reads_intra[,"flag"]==0,1])) #chromosomes to be checked
for(i in chromosome_names){
mat_temp <- cbind(a=reads_intra[reads_intra[,1]==i&reads_intra[,"flag"]==0,2],
                  b=reads_intra[reads_intra[,1]==i&reads_intra[,"flag"]==0,6])
mat_temp <- t(apply(mat_temp,1,FUN=function(x){return(x[order(x)])}))
  
enzyme_site <- countOverlaps(IRanges(start=mat_temp[,1], end=mat_temp[,2]), position_vec_plus[[i]], type="any")-1
reads_intra[reads_intra[,1]==i&reads_intra[,"flag"]==0,"flag"][enzyme_site>0] <- 1
}

####
#Fourth filter:
##
## remove a read pair in both MP mode and having distance < 1500bp, this is to reduce the chance of picking up incompletely digested self ligation products.

reads_intra[((reads_intra[,2]<reads_intra[,6]&reads_intra[,3]>reads_intra[,7])|
              (reads_intra[,2]>reads_intra[,6]&reads_intra[,3]<reads_intra[,7]))&
              (abs(reads_intra[,2]-reads_intra[,6])<=1500),"flag"] <- 0 


### combined filtered reads:

reads_filtered_intra <- reads_intra[reads_intra[,"flag"]==1,c(1:8)]
}else{reads_filtered_intra <- reads_intra[,1:8]}


reads_filtered <- rbind(reads_filtered_inter, reads_filtered_intra)
filename <- paste(mapped_reads, ".result", sep="")
write.table(reads_filtered[,c(1:3,5:7)], file=filename, row.names=F, col.names=F, quote=F, sep="\t") # read length is discarded.

print(paste("filtering----",mapped_reads, "----completed"))


q("no")
