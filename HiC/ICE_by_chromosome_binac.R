
options(stringsAsFactors=F, scipen=999)
args <- commandArgs(TRUE)
mapped_reads <- args[1] #filtered reads, each row has 6 columns
sample_id <- args[2]
#small windows, run ICE seperately on each chromosome:
bin <- 10000
bin_no_chr <- ceiling(c(30427671,19698289,23459830,18585056,26975502)/bin)

library(data.table)
HiC_reads <- fread(input=mapped_reads, header=F, sep='auto', data.table=F)

HiC_reads[,2] <- ceiling(HiC_reads[,2]/bin); HiC_reads[,5] <- ceiling(HiC_reads[,5]/bin)

suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
HiC_reads_split <- split(HiC_reads, HiC_reads[,1])
rm(HiC_reads); garbage <- gc(verbose = F)


IterativeCorNormalization <- function(x, max_iter=100, eps=1e-4){
  m <- dim(x)[1]
  
  ## Initialization    
  sum_ss <- matrix(rep(0, m), ncol=1)
  bias <- matrix(rep(1, m), ncol=1)
  old_dbias <- NULL
  ## Remove Diagonal ?
  
  for (it in 1:max_iter){
    message("it=",it," ", Sys.time())
    
    ## 1- calculate sum of W over all rows ++
    sum_ds <- rowSums(x, na.rm=TRUE)
    ##sum_ds <- sqrt(rowSums(x^2))
    
    ## 2- Calculate a vector of corrected ss reads
    ## NOT DONE
    
    ## 3- Calculate vector of bias
    dbias <- as.matrix(sum_ds, ncol=1) + sum_ss
    
    ## 4 - Renormalize bias by its mean valude over non-zero bins to avoid numerical instabilities
    dbias <- dbias/mean(dbias[dbias!=0])
    
    ## 5- Set zero values of bias to 1 to avoid 0/0 error
    dbias[dbias==0] <- 1
    
    ## 6- Divide W by bias BiBj for all (i,j) ++++
    x <- x/(dbias %*% t(dbias))
    
    ## 7- Multiple total vector of bias by additional biases
    ##bias <- bias * dbias
    
    if(it>1){message("current cycle: ", it, "; eps= ",sum(abs(old_dbias - dbias)))}
    
    if (!is.null(old_dbias) && sum(abs(old_dbias - dbias))<eps){
      message("Break at iteration ", it)
      break
    }
    old_dbias <- dbias 
  }
  if (it == max_iter){
    message("Did not converged. Stop at iteration ",max_iter)
  }else{
    message("Converged at ",it," iterations")
  }
  return(x)
}


registerDoParallel(cores=1)
ICE_by_chr <- foreach(i=1:5) %dopar% {
  count_table <- matrix(0, ncol=bin_no_chr[i], nrow=bin_no_chr[i])
  #count_table is a part of the big matrix
  
  for(j in 1:nrow(HiC_reads_split[[as.character(i)]])){
    chr_a <- i; chr_b <- HiC_reads_split[[as.character(i)]][j,4]
    if(chr_a!=chr_b){next}
    a <- HiC_reads_split[[as.character(i)]][j,2]; b <- HiC_reads_split[[as.character(i)]][j,5]
    count_table[a,b] <- count_table[a,b]+1
  }
  count_table <- count_table + t(count_table)
#  diag(count_table) <- diag(count_table)/2 #including teh diagnal
  diag(count_table) <- 0 #excluding the diagnal
  print(paste("filling chromosome",i, "--complete"))
  
  #start ICE
  bin_coverage <- apply(count_table, 1, sum)
  bins_removed <- numeric()
  stats <- boxplot(log(bin_coverage+1), plot=F)
  idx <- which(log(bin_coverage+1)>stats$stats[5,1])
  bins_removed <- c(bins_removed, idx)
  message("<<<<<<<<<filter_upper: ",length(idx)," >>>>>>>>>")
  idx <- which(log(bin_coverage+1)<stats$stats[1,1])
  bins_removed <- c(bins_removed, idx)
  message("<<<<<<<<<filter_lower: ",length(idx)," >>>>>>>>>")
  
  count_table[bins_removed,] <- 0; count_table[,bins_removed] <- 0
  normICE <- function(idata, max_iter=100, eps=1e-4){
    message("Start Iterative Correction ...")
    xmat <- IterativeCorNormalization(idata, max_iter=max_iter, eps=eps)
    return(xmat)
  }
  ICE_normalized <- normICE(count_table, max_iter=1500, eps=1e-4)
  #in case rounding is needed (to reduce file size)
  ICE_normalized <- round(ICE_normalized,3)
  
  return(ICE_normalized)
}
registerDoParallel(cores=1)

rm(HiC_reads_split); garbage <- gc(verbose = F)

#write normalized chromosomes:
for(i in 1:length(ICE_by_chr)){
  file_name <- paste("ICE_", sample_id, "_chr",i,"_binsize_",bin/1000,"kb", sep="")
  write.table(ICE_by_chr[[i]], file=file_name, sep="\t", row.names=F, col.names=F, quote=F)
}

