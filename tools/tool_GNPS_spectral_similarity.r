######################################
# GNPS spectral similarity calculation # 
######################################
library(clue)
# input:   x,y: data frame containing m/z and intensity in ms2 spectra;
#          x.pre,y.pre: precursor ion m/z
GNPS_score <- function(x, x.pre, y, y.pre, mz.tol=0.02){
        colnames(x) <- c("mz","int")
        colnames(y) <- c("mz","int")
        #delete pre
        remove_index <- which( abs(x$mz -x.pre) <= 0.02) 
        if (length(remove_index) !=0)  {
               x <- x[-remove_index,]
              if(nrow(x) == 0) {return(NA)} 
        }
        
        remove_index <- which( abs(y$mz -y.pre) <= 0.02) 
        if (length(remove_index) !=0) {
                y <- y[-remove_index,]
                if(nrow(y) == 0) {return(NA)}
        }
        
        dm <- y.pre - x.pre
        #square root transforms
        x[,2] <- sqrt(x[,2])
        y[,2] <- sqrt(y[,2])
        #scale intensity values
        x[,2] <- x[,2]/sqrt(sum(x[,2]^2))
        y[,2] <- y[,2]/sqrt(sum(y[,2]^2))
        # create alignment matrix
        alignment <- data.frame(matrix(ncol = 4))
        colnames(alignment) <- c("mz.x", "int.x", "int.y","mz.y")
        h <- 1
        for (m in 1:nrow(x)){
                mz.diff1 <- abs(x[m,1] - y[,1])
                mz.diff2 <- abs(x[m,1] - y[,1] + dm)
                if(min(mz.diff1) <= mz.tol){
                alignment[h,1] <- x[m,1]
                alignment[h,2] <- x[m,2]
                alignment[h,3] <- y[mz.diff1 == min(mz.diff1),2][1]
                alignment[h,4] <- y[mz.diff1 == min(mz.diff1),1][1]
                h <- h + 1
                }
                if(dm != 0){
                        if(min(mz.diff2) <= mz.tol){
                        alignment[h,1] <- x[m,1]
                        alignment[h,2] <- x[m,2]
                        alignment[h,3] <- y[mz.diff2 == min(mz.diff2),2][1]
                        alignment[h,4] <- y[mz.diff2 == min(mz.diff2),1][1]
                        h <- h + 1
                        }
                }
        }
        alignment <- alignment[complete.cases(alignment),]

        # create a matrix for selecting the highest scoring subset of matching peaks using 'clue' library
        # each peak is matched to at most one peak in another spectrum
        if(nrow(alignment)==0){score <- 0 }
        if(nrow(alignment)>0){
                mzfrag.x <- unique(alignment$mz.x)
                mzfrag.y <- unique(alignment$mz.y)
                max.length <- max(length(mzfrag.x),length(mzfrag.y))
                matrix <- matrix(0, nrow = max.length, ncol = max.length) #matrix nrow=ncol
        #name the matrix rows and cols with mz
        if(length(mzfrag.x)<=length(mzfrag.y)){
        rownames(matrix) <- c(mzfrag.x, rep(0,(max.length-length(mzfrag.x))))
        colnames(matrix) <- mzfrag.y
        }
        if(length(mzfrag.x)>length(mzfrag.y)){
        rownames(matrix) <- mzfrag.x
        colnames(matrix) <- c(mzfrag.y,rep(0,(max.length-length(mzfrag.y))))
    }

    #fill the matrix with intensity product
    for(h in 1:nrow(alignment)){
      matrix[mzfrag.x==alignment[h,1],mzfrag.y==alignment[h,4]] <- alignment[h,2] * alignment[h,3]
    }
    if(length(mzfrag.x)<length(mzfrag.y)){matrix[(length(mzfrag.x)+1):max.length,]<-0}
    if(length(mzfrag.x)>length(mzfrag.y)){matrix[,(length(mzfrag.y)+1):max.length]<-0}
    #LSAP problem in 'clue' library
    optimal <- solve_LSAP(matrix, maximum = TRUE)
    #calculate GNPS score
    AB <- 0
    for (n in 1:max.length){AB <- AB + matrix[n,optimal[n]] }
    GNPS <- AB
    score <- as.numeric(GNPS)
  }
  return(score)
}


