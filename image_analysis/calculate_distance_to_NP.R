



# this R script does following:
# 1) extract nuclear boundary
# 2) calculate the distance of each pixel inside the nucleus from the nuclear periphery
# 3) correlate the confocal signals (Green/Red channels) to the distances obtained from step 2.

# load library
options(stringsAsFactors=F, scipen=999)
library(readxl)
setwd("~/Documents/ZMBP/DNA_methylation_and_NP/FISH/for_distance_calculation")

# read the xlsx file, get informatino of sheets:
input_files <- c("dapi_col_2.xlsx", "green_col_2.xlsx", "red_col_2.xlsx") # make sure to follow this sequence: DAPI, Green, Red
slice_no <- length(excel_sheets(path = input_files[1]))
nucleus_slice <- nucleus_border_slice <- distance_slice <- DAPI_data <- list()

for(s in 1:slice_no){
  #read DAPI image
  DAPI_data[[s]] <- read_xlsx(input_files[1], sheet=s, col_names=F)
  DAPI_data[[s]] <- data.matrix(DAPI_data[[s]])
  
  #fill nucleus
  row_number <- nrow(DAPI_data[[s]]); col_number <- ncol(DAPI_data[[s]])
  nucleus <- surface <- matrix(0, row_number, col_number)
  for(i in 1:row_number){
    for(j in 1:col_number){
      if(sum(DAPI_data[[s]][i,1:j])>=1&sum(DAPI_data[[s]][i,j:col_number])>=1&sum(DAPI_data[[s]][1:i,j])>=1&sum(DAPI_data[[s]][i:row_number,j])>=1){
        nucleus[i,j] <- 1 # label all pixels overlapping with DAPI-stained regions
      }
    }
  }
  
  #find the surface of this nuclei
  #NOTE: here I assume that this nucleus does not have a irregular shape !
  nuclei_pt <- which(nucleus==1,arr.ind=T)
  for(i in 1:nrow(nuclei_pt)){
    if((nucleus[nuclei_pt[i,1]-1,nuclei_pt[i,2]]+nucleus[nuclei_pt[i,1]+1,nuclei_pt[i,2]]+
        nucleus[nuclei_pt[i,1],nuclei_pt[i,2]-1]+nucleus[nuclei_pt[i,1],nuclei_pt[i,2]+1])<4){
      surface[nuclei_pt[i,1],nuclei_pt[i,2]] <- 1
    }
  }
  nucleus_slice[[s]] <- nucleus; nucleus_border_slice[[s]] <- surface
}

par(mfrow=c(3,4), pty="s", mar=c(1,1,1,1))
r1<-1; r2<-nrow(DAPI_data[[1]]); c1<-1; c2<-nrow(DAPI_data[[1]])
for(s in 1:slice_no){
  if(s==1){nucleus_border_slice_total <- nucleus_border_slice[[s]]}else{
    nucleus_border_slice_total <- nucleus_border_slice_total+nucleus_border_slice[[s]]
  }
  image(nucleus_border_slice[[s]][r1:r2, c1:c2], col=colorRampPalette(c("blue","red","red","yellow"))(100))
}
image(nucleus_border_slice_total[r1:r2, c1:c2], col=colorRampPalette(c("blue","red","red","yellow"))(100))


boundaries <- c()
for(s in 1:slice_no){
  boundaries_pixel <- which(nucleus_border_slice[[s]]==1, arr.ind=T)
  boundaries_pixel <- cbind(boundaries_pixel, s*0.24/0.07199219, s) #1 pixel=0.03295 um
  boundaries <- rbind(boundaries, boundaries_pixel)
}

# To check all slices:
# library(rgl)
# plot3d(boundaries[,1:3],xlim =c(1,nrow(DAPI_data[[1]])),ylim =c(1,nrow(DAPI_data[[1]])),zlim = c(-nrow(DAPI_data[[1]])/2,nrow(DAPI_data[[1]])/2), size=2)
# detach(package:rgl, unload=T)


dist_two_points <- function(delta_row, delta_col, delta_z){
  distance <- sqrt(delta_row^2+delta_col^2+delta_z^2)
  return(round(distance,2))
}

#calculate the distances between each pixel inside a nuclear slice and all nuclear borders, then choose the smallest as the distance of this pixel to NP:
#NOTE: the three slices close to glass slide were not considered (due to distortion)

####
####
#transform physical distance into pixel:
h <- 0.24; scale <- 0.07199219 #h: distance between slice (um); scale: (um/pixel)

#define slices:
for(s in 1:(slice_no-0)){
  nucleus_pixel <- which(nucleus_slice[[s]]==1, arr.ind=T)
  distance_slice[[s]] <- nucleus_slice[[s]]; distance_slice[[s]][,] <- 0
  for(i in 1:nrow(nucleus_pixel)){
    row_id <- nucleus_pixel[i,1]; col_id <- nucleus_pixel[i,2]
    distances_current <- apply(boundaries[boundaries[,4]==s,], 1 , FUN=function(x){dist_two_points(x[1]-row_id, x[2]-col_id, (x[4]-s)*h/scale)})
    distances_current_min <- min(distances_current)
    
    #considering other slices, ONLY if the vertical distance is shorter than distances_current_min
    delta_max <- floor(distances_current_min/(h/scale))
    if(delta_max>0){
      more_s <- c(max(1, s-delta_max):(s-1), (s+1):min(slice_no-0, s+delta_max))
      more_s <- more_s[more_s>0]
      distances <- apply(boundaries[boundaries[,4]%in%more_s,], 1 , FUN=function(x){dist_two_points(x[1]-row_id, x[2]-col_id, (x[4]-s)*h/scale)})
      distance_slice[[s]][row_id, col_id] <- min(distances_current_min, distances)
    }else{distance_slice[[s]][row_id, col_id] <- distances_current_min}
  }
  print(s)
}

par(mfrow=c(3,4), pty="s", mar=c(1,1,1,1))
for(s in 1:length(distance_slice)){
  image(distance_slice[[s]][r1:r2, c1:c2], col=colorRampPalette(c("blue","red","yellow"))(100))
}


# read green and red signals
green <- red <- list()
for(s in 1:slice_no){
  green[[s]] <- read_xlsx(input_files[2], sheet=s, col_names=F)
  green[[s]] <- data.matrix(green[[s]])
  red[[s]] <- read_xlsx(input_files[3], sheet=s, col_names=F)
  red[[s]] <- data.matrix(red[[s]])
  for(i in 1:row_number){
    for(j in 1:col_number){
      if(nucleus_slice[[s]][i,j]==0){
        green[[s]][i,j] <-NA # this pixel is outside the nuclear slice
        red[[s]][i,j] <-NA # this pixel is outside the nuclear slice
      }
    }
  }
}

par(mfrow=c(2,2), pty="s", mar=c(4,4,2,2))
# draw the plot
red_vector<- as.vector(unlist(red))
green_vector<- as.vector(unlist(green))
distance_vector<- as.vector(unlist(distance_slice))*scale; distance_vector <- distance_vector[!is.na(green_vector)]
red_vector <- red_vector[!is.na(red_vector)]
green_vector <- green_vector[!is.na(green_vector)]


###: keep the pixels with top 5% values
cutoff <- 0.95

plot(red_vector[red_vector>quantile(red_vector,cutoff)], distance_vector[red_vector>quantile(red_vector,cutoff)],
     xlab="gray value",ylab = "distance(μm)",xlim = c(0,250), ylim=c(0,2), pch=16,cex=0.3,col = "red")
plot(green_vector[green_vector>quantile(green_vector,cutoff)], distance_vector[green_vector>quantile(green_vector,cutoff)],
     xlab="gray value",ylab = "distance(μm)",xlim = c(0,250), ylim=c(0,2), pch=16,cex=0.3,col = "green3")

boxplot(distance_vector[red_vector>quantile(red_vector,cutoff)], distance_vector[green_vector>quantile(green_vector,cutoff)],
        col = c("red","green3"),xlab="probe",ylab="distance")
axis(1,at=1:2,labels = c("red","green"))

wilcox.test(distance_vector[red_vector>quantile(red_vector,cutoff)],distance_vector[green_vector>quantile(green_vector,cutoff)])

# cumulative distribution function:
plot(ecdf(distance_vector[red_vector>quantile(red_vector,cutoff)]), col="red", 
     cex=0.5, ylab="Cumulative density",xlab="Distance", xlim=c(0, 2), do.points = FALSE, verticals = TRUE)
plot(ecdf(distance_vector[green_vector>quantile(green_vector,cutoff)]), col="green3", cex=0.5, add=T,
     do.points = FALSE, verticals = TRUE)
