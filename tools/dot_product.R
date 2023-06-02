#########  dot product function ########
# dot product #
# exp_fragMz, exp_fragInt, ref_fragMz, ref_fragMzInt: separate by ";"
dotProduct <- function(exp_fragMz, exp_fragInt, ref_fragMz, ref_fragMzInt, ms2_tol =0.02 ) {
        # get the vector of experimental fragMz and fragInt
        exp_fragMz <- as.numeric(strsplit(exp_fragMz, ";")[[1]])
        exp_fragInt <- as.numeric(strsplit(exp_fragInt, ";")[[1]])
        exp_fragInt <- 100*exp_fragInt/max(exp_fragInt)
        # calculate the length of intensity vector
        A <- sqrt(sum(exp_fragInt^2))
        # record how many experimental fragment
        totalFrag <- length(exp_fragMz)
        # find the reference frag_mz and frag_int
        ref_fragMz <- as.numeric(strsplit(ref_fragMz, ";")[[1]])
        ref_fragInt <- as.numeric(strsplit(ref_fragMzInt, ";")[[1]])
        ref_fragInt <- 100*ref_fragInt/max(ref_fragInt)
        # calculate the length of intensity vector
        B <- sqrt(sum(ref_fragInt^2))
        # create the aligned data frame
        int_dataframe <- data.frame(matrix(nrow = length(ref_fragMz),ncol=4))
        colnames(int_dataframe) <- c("ref_fragMz","ref_fragInt", "exp_fragInt","exp_fragMz")
        int_dataframe[,1] <- as.data.frame(ref_fragMz)
        int_dataframe[,2] <- as.data.frame(ref_fragInt)
        count <- 0 # record how many counts found
        for (j in 1:length(ref_fragMz)){
                fragMz_index <- which( exp_fragMz %in% exp_fragMz[abs(exp_fragMz - int_dataframe[j,1]) <= ms2_tol] )
                # no frag found in t_fragment, assign 0
                if ( length(fragMz_index)==0) {
                        int_dataframe[j,3] <- 0
                        next
                }
                # assign the found intensity to dataframe
                frag_index <- fragMz_index[which.max(exp_fragInt[fragMz_index])]
                int_dataframe[j,3] <- exp_fragInt[frag_index[1]]
                int_dataframe[j,4] <- exp_fragMz[frag_index[1]]
                count <- count + 1
                #remove the matched fragment from experimental
                exp_fragMz <- exp_fragMz[-frag_index]
                exp_fragInt <- exp_fragInt[-frag_index]
        }
        # calculate the dp
        dp <- sum(int_dataframe$ref_fragInt * int_dataframe$exp_fragInt)
        mz <- paste0(int_dataframe$ref_fragMz[complete.cases(int_dataframe$exp_fragMz)],
                     collapse = ";")
        score <- 100*dp/(A*B)
        results <- c(score, count, mz)
        # return the number of matched fragments and dot product score
        return(results)
}
