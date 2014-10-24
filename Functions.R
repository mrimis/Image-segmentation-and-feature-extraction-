##Function for Kmean clustering algorithm implementation on Image data## 

kmeansCluster <- function(img, centers, iter.max=10, nstart=1, algorithm = c("Hartigan-Wong","Lloyd", "Forgy", "MacQueen")){
  clust <- kmeans(as.vector(img), centers, iter.max, nstart, algorithm)
  result <- img
  greatest = -1
  greatest.index = -1
  sorted.clusters = sort(clust$cluster, decreasing=T)
  for (i in 1:centers){
    indices <- which(clust$cluster == i)
    result[indices] <- i
    }
  numbers <- seq(1,centers)
  sizes <- clust$size
  sorted.sizes <- sort(clust$size, decreasing=T)
  mat <- matrix(c(numbers, sizes, sorted.sizes), ncol=3)
  temp <- result
  if (all(mat[,2] %in% mat[,3])){
    for (i in 1:centers){
      if (mat[i,2] != mat[i,3]){
        ind <- which(mat[,3] == mat[i,2], arr.ind=T)
        ind <- ind[1]
        result[which(temp == i)] <- ind
        }
      }
    }
  return(result)
}


# Function for segmenting cell images and extracting Shape,size and texture related features
library(biOps)
library(EBImage)
library(alphahull)
segment_image <- function(path_in,path_out,green_channel=F,red_channel=F,day_nmbr,Drug_name){
  all_files = list.files(path_in)
  image_ix = grep(pattern=".tif", x=all_files)
  image_names = all_files[image_ix]
  min_size=100
  for(i in image_names){
    image1 <-  readImage(paste(path_in, i, sep=""), type="tiff");image2 <- flip(image1)
    if(green_channel){image<- kmeansCluster(image1[,,2], 10);image <- flip(normalize(image))}
    else if(red_channel){image<- kmeansCluster(image1[,,1], 7);image <- flip(normalize(image))}
    y <- image
    y[y > 0.1] <- 1
    y[y <= 0.1] <-0
    image_label = propagate(image, seeds=watershed(distmap(fillHull(bwlabel(opening(y, makeBrush(7, shape='disc'))))),9,1), mask=opening(image>0.06, makeBrush(5, shape='disc')))
    size <- computeFeatures.shape(image_label, properties=FALSE)
    intensity <- computeFeatures.basic(image_label, image, properties=FALSE)
    spatial <- computeFeatures.moment(image_label, properties=FALSE)
    txtr <- computeFeatures.haralick(image_label,image,properties=FALSE)
    s.area <- size[,1]
    s.perimeter <- size [,2]
    s.radius.mean <- size [,3]
    s.radius.min <- size [,4]
    s.radius.max <- size[,5]
    s.circularity <- (4*pi*size[,1]*1.0)/(size[,2]*size[,2])
    intensity[,1] <- intensity[,1]*256
    intensity_mean <- intensity[,1]
    m.cx <- spatial[,1]
    m.cy <- spatial[,2]
    m.majoraxis <- spatial[,3]
    m.roundness <- spatial[,4]
    m.theta <- spatial[,5]
    texture_features <- txtr[,1:25]
    h.asm.s1 <- txtr[,1];h.con.s1<- txtr[,2];h.cor.s1<- txtr[,3];h.var.s1<- txtr[,4];h.idm.s1<- txtr[,5];h.sav.s1<- txtr[,6];h.sva.s1<- txtr[,7];h.sen.s1<- txtr[,8];h.ent.s1<- txtr[,9];h.dva.s1<- txtr[,10];h.den.s1<- txtr[,11];h.f12.s1<- txtr[,12];h.f13.s1<- txtr[,13];h.asm.s2<- txtr[,14];h.con.s2<- txtr[,15];h.cor.s2<- txtr[,16];h.var.s2<- txtr[,17];h.idm.s2<- txtr[,18];h.sav.s2<- txtr[,19];h.sva.s2<- txtr[,20];h.sen.s2<- txtr[,21];h.ent.s2<- txtr[,22];h.dva.s2<- txtr[,23];h.den.s2<- txtr[,24];h.f12.s2<- txtr[,25]
    h.f13.s2 <- txtr[,26]
    results <- cbind(intensity_mean,m.cx,m.cy,m.majoraxis,m.roundness,m.theta,s.area,s.perimeter,s.radius.mean,s.radius.min,s.radius.max,s.circularity,h.asm.s1,h.con.s1,h.cor.s1,h.var.s1,h.idm.s1,h.sav.s1,h.sva.s1,h.sen.s1,h.ent.s1,h.dva.s1,h.den.s1,h.f12.s1,h.f13.s1,h.asm.s2,h.con.s2,h.cor.s2,h.var.s2,h.idm.s2,h.sav.s2,h.sva.s2,h.sen.s2,h.ent.s2,h.dva.s2,h.den.s2,h.f12.s2,h.f13.s2)
    results_df <- as.data.frame(results)
    if (green_channel){
      image_label_clean <- rmObjects(image_label,which(results[,7] < min_size , ))
      colorMode(image_label_clean)=Grayscale
      if (!max(image_label_clean)==0){
        results_clean <- subset(results_df, results[,7] > min_size , select=1:38)
        results_clean$label <- seq(1:length(results_clean[,1]))
		
		if (!day_nmbr == "Not valid" & !Drug_name == "Not valid"){
		results_clean$day <- rep(day_nmbr,length(results_clean[,1]))
		results_clean$treatment_drug <- rep(Drug_name,length(results_clean[,1]))
		results_clean <- results_clean[c("label","day","treatment_drug","intensity_mean","m.cx","m.cy","m.majoraxis","m.roundness","m.theta","s.area","s.perimeter","s.radius.mean","s.radius.min","s.radius.max","s.circularity","h.asm.s1","h.con.s1","h.cor.s1","h.var.s1","h.idm.s1","h.sav.s1","h.sva.s1","h.sen.s1","h.ent.s1","h.dva.s1","h.den.s1","h.f12.s1","h.f13.s1","h.asm.s2","h.con.s2","h.cor.s2","h.var.s2","h.idm.s2","h.sav.s2","h.sva.s2","h.sen.s2","h.ent.s2","h.dva.s2","h.den.s2","h.f12.s2","h.f13.s2")]
		}
		else if (day_nmbr == "Not valid" & Drug_name == "Not valid") {
		results_clean <- results_clean[c("label","intensity_mean","m.cx","m.cy","m.majoraxis","m.roundness","m.theta","s.area","s.perimeter","s.radius.mean","s.radius.min","s.radius.max","s.circularity","h.asm.s1","h.con.s1","h.cor.s1","h.var.s1","h.idm.s1","h.sav.s1","h.sva.s1","h.sen.s1","h.ent.s1","h.dva.s1","h.den.s1","h.f12.s1","h.f13.s1","h.asm.s2","h.con.s2","h.cor.s2","h.var.s2","h.idm.s2","h.sav.s2","h.sva.s2","h.sen.s2","h.ent.s2","h.dva.s2","h.den.s2","h.f12.s2","h.f13.s2")]
		}
		else if (day_nmbr == "Not valid" & !Drug_name == "Not valid"){
		results_clean$treatment_drug <- rep(Drug_name,length(results_clean[,1]))
		results_clean <- results_clean[c("label","treatment_drug","intensity_mean","m.cx","m.cy","m.majoraxis","m.roundness","m.theta","s.area","s.perimeter","s.radius.mean","s.radius.min","s.radius.max","s.circularity","h.asm.s1","h.con.s1","h.cor.s1","h.var.s1","h.idm.s1","h.sav.s1","h.sva.s1","h.sen.s1","h.ent.s1","h.dva.s1","h.den.s1","h.f12.s1","h.f13.s1","h.asm.s2","h.con.s2","h.cor.s2","h.var.s2","h.idm.s2","h.sav.s2","h.sva.s2","h.sen.s2","h.ent.s2","h.dva.s2","h.den.s2","h.f12.s2","h.f13.s2")]
		}
		else if (!day_nmbr == "Not valid" & Drug_name == "Not valid"){
		results_clean$day <- rep(day_nmbr,length(results_clean[,1]))
		results_clean <- results_clean[c("label","day","intensity_mean","m.cx","m.cy","m.majoraxis","m.roundness","m.theta","s.area","s.perimeter","s.radius.mean","s.radius.min","s.radius.max","s.circularity","h.asm.s1","h.con.s1","h.cor.s1","h.var.s1","h.idm.s1","h.sav.s1","h.sva.s1","h.sen.s1","h.ent.s1","h.dva.s1","h.den.s1","h.f12.s1","h.f13.s1","h.asm.s2","h.con.s2","h.cor.s2","h.var.s2","h.idm.s2","h.sav.s2","h.sva.s2","h.sen.s2","h.ent.s2","h.dva.s2","h.den.s2","h.f12.s2","h.f13.s2")]
		}
        overlayObjects <- paintObjects(image_label_clean, image2, col=c('red'))
        overlayObjects <- flip(overlayObjects)    
        width <- length(image[,1])
        height <- length(image[1,])
        write.table(results_clean, paste(path_out,"green_extracted_features_",i,"_results.txt",sep=""), sep="\t",row.names = FALSE, col.names=TRUE)
        tiff(filename = paste(path_out,"overlay_","green_segmented_image_",i,sep=""), width = width, height = height, pointsize = 12,compression="none")
        par(mar=c(0,0,0,0))
        plot(x = NULL, y = NULL, xlim = c(0,as.numeric(width)), ylim = c(0,as.numeric(height)), pch = '', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xaxs = 'i', yaxs = 'i', bty = 'n') # plot empty figure
        rasterImage(overlayObjects, xleft = 0, ybottom = 0, xright = width, ytop = height) 
        text(results_clean$m.cx,results_clean$m.cy,labels=results_clean$label, col='yellow')
        dev.off()
      }
    }
    else{
      results_df$label <- seq(1:length(results_df[,1]))
      results_df <- results_df[c("label","intensity_mean","m.cx","m.cy","m.majoraxis","m.roundness","m.theta","s.area","s.perimeter","s.radius.mean","s.radius.min","s.radius.max","s.circularity","h.asm.s1","h.con.s1","h.cor.s1","h.var.s1","h.idm.s1","h.sav.s1","h.sva.s1","h.sen.s1","h.ent.s1","h.dva.s1","h.den.s1","h.f12.s1","h.f13.s1","h.asm.s2","h.con.s2","h.cor.s2","h.var.s2","h.idm.s2","h.sav.s2","h.sva.s2","h.sen.s2","h.ent.s2","h.dva.s2","h.den.s2","h.f12.s2","h.f13.s2")]
      colorMode(image_label)=Grayscale
      overlayObjects <- paintObjects(image_label, image2, col=c('red'))
      overlayObjects <- flip(overlayObjects)    
      width <- length(image[,1])
      height <- length(image[1,])
      write.table(results_df, paste(path_out,"red_extracted_features_",i,"_results.csv",sep=""),sep="\t",row.names = FALSE, col.names=TRUE)
      tiff(filename = paste(path_out,"overlay_","red_segmented_image_",i,sep=""), width = width, height = height, pointsize = 12,compression="none")
      par(mar=c(0,0,0,0))
      plot(x = NULL, y = NULL, xlim = c(0,as.numeric(width)), ylim = c(0,as.numeric(height)), pch = '', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xaxs = 'i', yaxs = 'i', bty = 'n') # plot empty figure
      rasterImage(overlayObjects, xleft = 0, ybottom = 0, xright = width, ytop = height) 
      text(results_df$m.cx,results_df$m.cy,labels=results_df$label, col='yellow')
      dev.off()
	  }
  }
}
