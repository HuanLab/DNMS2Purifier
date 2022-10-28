# DNMS2Purifier
# This R script reads in the mzML or mzXML files of DDA data,
#               outputs csv files and mgf files containing the metabolic features table with purified ms2 spectra
# Tingting Zhao, Oct 3, 2022
# Copyright @ The University of British Columbia

working_directory <- getwd()
#--------------- parameters setting --------------------------------------------
assign_rt_tol <- 15      # retention time tolerance (in seconds) for MS2 assignment
assign_ms1mz_tol <- 0.01 # m/z tolerance for MS2 assignment
ratio_ms2mz_tol <- 0.01  # m/z tolerance when calculating the ratio for each fragment ion
int_threshold <- 0.01    # threshold of relative intensity to remove the low abundant fragments 
target_extraction <- TRUE # TRUE:  Purification will be performed on target metabolic features
                           # FALSE: Purification will be performed on all metabolic features
mim_num_spectra <- 3   # the least number of spectra required for each feature to be qualified for purification

#--------- load library --------------
library(xcms) # load xcms package
library(xgboost) # load xgboost package

#---------------------- function to remove H isotope ---------------------
remove_iso <- function(frag_v,int_v) {
        match_ms2mz_tol <- 0.01
        for (j in 1:(length(frag_v)-1)) {
                t_mz <- frag_v[j]
                iso1_index <- vector()
                iso2_index <- vector()
                iso3_index <- vector()
                iso4_index <- vector()
                # check isotope by mz, record the index of satisfying the mz difference (1 isotop)
                if ( min(abs(frag_v - t_mz -1.003)) <= match_ms2mz_tol) {
                        iso1_mz <- which (frag_v %in% frag_v[abs(frag_v - t_mz - 1.003) <= match_ms2mz_tol])
                        # confirm the intensity of the isotope is lower, otherwise skip to next fragment
                        for (n in 1:length(iso1_mz)) {
                                if (int_v[iso1_mz[n]] < int_v[j]){
                                        iso1_index <- c(iso1_index,iso1_mz[n])
                                }
                        }
                } 
                
                if (length(iso1_index) != 0) {
                        # check if there is a 2nd isotope
                        if (min(abs(frag_v - t_mz - 2.006)) <= match_ms2mz_tol) {
                                iso2_mz <- which( frag_v %in% frag_v[abs(frag_v - t_mz - 2.006) <= match_ms2mz_tol])
                                
                                for (n in 1:length(iso2_mz)) {
                                        if (int_v[iso2_mz[n]] < min(int_v[iso1_index]) ){
                                                iso2_index <- c(iso2_index,iso2_mz[n])
                                        }
                                }
                        }
                }
                
                # check if there is a 3rd isotope
                if (length(iso2_index) !=0 ){
                        if (min(abs(frag_v - t_mz - 3.009)) <= match_ms2mz_tol) {
                                iso3_mz <- which( frag_v %in% frag_v[abs(frag_v - t_mz - 3.009) <= match_ms2mz_tol])
                                
                                for (n in 1:length(iso3_mz)) {
                                        
                                        if (int_v[iso3_mz[n]] < min(int_v[iso2_index]) ){
                                                iso3_index <- c(iso3_index,iso3_mz[n])
                                        }
                                }
                        }
                }
                
                # check if this fragment has 4th isotope
                if (length(iso3_index) !=0 ){
                        if (min(abs(frag_v - t_mz - 4.012)) <= match_ms2mz_tol) {
                                iso4_mz <- which( frag_v %in% frag_v[abs(frag_v - t_mz - 4.012) <= match_ms2mz_tol])
                                
                                for (n in 1:length(iso4_mz)) {
                                        
                                        if (int_v[iso4_mz[n]] < min(int_v[iso3_index]) ){
                                                iso4_index <- c(iso4_index,iso4_mz[n])
                                        }
                                }
                                
                        }
                }
                iso_list <- c(iso1_index, iso2_index,iso3_index, iso4_index)
                if (length(iso_list) != 0) {
                        frag_v <- frag_v[-iso_list]
                        int_v <- int_v[-iso_list] }
                # if the j has close to the size of removed frag_v, stop the loop
                if ( j == length(frag_v) ) {break}
        }
        frag_v <- round(frag_v, digits = 4)
        frag_wo_iso <- paste0(frag_v, collapse = ";")
        int_wo_iso <- paste0(int_v, collapse = ";")
        return  (c(frag_wo_iso, int_wo_iso))
}
#---------------------------- main program ---------------------------------
setwd(working_directory)
#------------- feature extraction -------------------
# Load all the .mzXML files or .mzML
DDAfiles <- list.files(pattern = ".mzXML") 
N <- length(DDAfiles)
if (N == 0) {
        DDAfiles <- list.files(pattern = ".mzML")
        N <- length(DDAfiles)
        message(paste0(N, " of mzML raw data have been loaded :)"))
}
message(paste0(N, " of mzXML raw data have been loaded :)"))

if (!target_extraction) {
        message("Start to extract metabolic features")
        msSet <- xcmsSet(DDAfiles,
                         method = "centWave",
                         ppm=12, 
                         peakwidth= c(10,120),
                         mzdiff= 0.01, 
                         snthresh=6, 
                         integrate=1,
                         prefilter= c(3, 100), 
                         noise=100)
        
        # Multiple files need to be aligned
        msSet <- group(msSet,
                       bw=5, 
                       minfrac=0.5,
                       mzwid=0.015, 
                       minsamp=1, 
                       max=50)
        
        msSet <- retcor(msSet, 
                        method= "obiwarp", 
                        profStep=1)
        
        msSet <- group(msSet, 
                       bw=5, 
                       minfrac=0.5,
                       mzwid=0.015, 
                       minsamp=1, 
                       max=50)
        
        msSet <- fillPeaks(msSet)
        
        xc_mzrt <- data.frame(msSet@groups)[,c(1,4)]
        
        ## extract the peak max intensity in each sample
        xc_int <- groupval(msSet, value= "maxo")
        # extract feature table
        featureTable <- as.data.frame(cbind(0,xc_mzrt,xc_int))
        colnames(featureTable)[1:3] <- c("ID", "mz", "RT")
        featureTable$ID <- 1:nrow(featureTable)
        ft <- featureTable
        message("feature table has been extracted :)")
}else{
        # if target_extraction == TRUE, the user needs to update ft_name with name of feature table
        ft_name <- "target_table.csv" 
        target_ft <- read.csv(ft_name)
        message("target feature table has been loaded :)")
}

#-------- ms2 spectra assignment ---------------------
# read mzXML file mass spectrum, save the precursor_v, and rt_v
for (i in 1:N){
        if (i==1){t0 <- Sys.time()}
        
        msfile <- DDAfiles[[i]]
        msfile <- readMSData(msfile, msLevel. = 2, centroided. = T)
        ## record msfile_i
        assign(paste0("msfile_",i), msfile)
        
        ## record pre_v, rt_v, preInt_v
        pre_v <- vector()
        preInt_v <- vector()
        rt_v <- vector()
        
        for (n in 1:length(msfile)){
                
                pre_v <- c(pre_v,msfile[[n]]@precursorMz)
                preInt_v <- c(preInt_v, msfile[[n]]@precursorIntensity)
                rt_v <- c(rt_v, msfile[[n]]@rt)
        }
        
        assign(paste0("pre_v_",i), pre_v)
        assign(paste0("rt_v_",i), rt_v)
        assign(paste0("preInt_v_",i), preInt_v) 
        # if (i == length(DDAfiles)) {print(Sys.time() - t0)}
}
# Assign ms2 to each feature
message("Start assigning MS2 spectra to each feature :)")
if (!target_extraction) {
        for (i in 1:N) {
                ft[,N+3+i] <- NA
                colnames(ft)[N+3+i] <- paste0("ms2mz_",i)
        }
        for (i in 1:N) {
                ft[,N+N+3+i] <- NA
                colnames(ft)[N+N+3+i] <- paste0("ms2Int_",i)
        }
        
        for (i in 1:nrow(ft)) {
                
                t_pre <- ft$mz[i]
                t_rt <- ft$RT[i]
                ## obtain the ms2 in each file
                for (n in 1:N){
                        ### recall msfile_n, pre_v_n, rt_v_n
                        msfile <- get(paste0("msfile_",n))
                        pre_v <- get(paste0("pre_v_",n))
                        preInt_v <- get(paste0("preInt_v_",n))
                        rt_v <- get(paste0("rt_v_",n))
                        
                        ### find the matched precursor 
                        rt_index <- which (rt_v %in% rt_v[abs(rt_v - t_rt) <= assign_rt_tol])
                        selected_pre <- pre_v[rt_index]
                        pre_index <- which (selected_pre %in% selected_pre[abs(selected_pre - t_pre)<= assign_ms1mz_tol])
                        
                        if (length(pre_index)== 0) {next}
                        selected_int <- preInt_v[rt_index][pre_index]
                        int_index <- which(selected_int %in% max(selected_int))[1]
                        
                        numms2 <- 1:length(msfile)
                        final_index <- numms2[rt_index][pre_index][int_index]
                        ms2mz <- msfile[[final_index]]@mz
                        ms2Int <- msfile[[final_index]]@intensity
                        
                        col <- n+ N +3
                        ms2mz_list <- paste0(ms2mz, collapse = ";")
                        ms2Int_list <- paste0(ms2Int, collapse = ";")
                        ft[i, col] <- ifelse(ms2mz_list == "", NA, ms2mz_list)
                        ft[i, col+N] <- ifelse(ms2Int_list == "", NA, ms2Int_list)
                }
        }
}else{  
        # assign ms2 spectra to the feature table the user provides.
        if (ncol(target_ft) == 3){
                # re-frame the feature table to save the ms2 information 
                for (i in 1:N) {
                        target_ft[,3+i] <- NA
                        colnames(target_ft)[3+i] <- DDAfiles[i]
                }
                
                for (i in 1:N) {
                        target_ft[,N+3+i] <- NA
                        colnames(target_ft)[N+3+i] <- paste0("ms2mz_",i)
                }
                
                for (i in 1:N) {
                        target_ft[,N+N+3+i] <- NA
                        colnames(target_ft)[N+N+3+i] <- paste0("ms2Int_",i)
                }
                
                for (i in 1:nrow(target_ft)) {
                        t_pre <- target_ft$litera_mz[i]
                        t_rt <- as.numeric(target_ft$RT[i])
                        
                        ## extract ms2 information(precursor intensity, ms2 m/z, ms2 intensity) in each file
                        for (n in 1:N){
                                ### recall msfile_n, pre_v_n, rt_v_n
                                msfile <- get(paste0("msfile_",n))
                                pre_v <- get(paste0("pre_v_",n))
                                preInt_v <- get(paste0("preInt_v_",n))
                                rt_v <- get(paste0("rt_v_",n))
                                
                                ### find the matched precursor 
                                rt_index <- which (rt_v %in% rt_v[abs(rt_v - t_rt) <= assign_rt_tol ])
                                if (length(rt_index) == 0) {
                                        next
                                }
                                
                                selected_pre <- pre_v[rt_index]
                                pre_index <- which (selected_pre %in% selected_pre[abs(selected_pre - t_pre)<= assign_ms1mz_tol])
                                if (length(pre_index)== 0) {next}
                                if (length(pre_index) == 1) {
                                        int_index <- 1
                                } else {
                                        selected_int <- preInt_v[rt_index][pre_index]
                                        int_index <- which(selected_int %in% max(selected_int))[1]
                                }
                                
                                numMS2 <- 1:length(msfile)
                                final_index <- numMS2[rt_index][pre_index][int_index]
                                ms2mz <- round(msfile[[final_index]]@mz,4)
                                ms2Int <- msfile[[final_index]]@intensity
                                ms2preInt <- msfile[[final_index]]@precursorIntensity
                                
                                col <- n + N +3
                                ms2mz_list <- paste0(ms2mz, collapse = ";")
                                ms2Int_list <- paste0(ms2Int, collapse = ";")
                                target_ft[i,3 + n] <- ms2preInt
                                target_ft[i, col] <- ifelse(ms2mz_list == "", NA, ms2mz_list)
                                target_ft[i, col+N] <- ifelse(ms2Int_list == "", NA, ms2Int_list)
                        }
                        
                }
                ft <- target_ft
        }else{
                ft <- target_ft
                
                for (i in 1:N) {
                        target_ft[,N+3+i] <- NA
                        colnames(target_ft)[N+3+i] <- paste0("ms2mz_",i)
                }
                
                for (i in 1:N) {
                        target_ft[,N+N+3+i] <- NA
                        colnames(target_ft)[N+N+3+i] <- paste0("ms2Int_",i)
                }
                
                for (i in 1:nrow(ft)) {
                        
                        t_pre <- ft$mz[i]
                        t_rt <- ft$RT[i]
                        ## obtain the ms2 in each file
                        for (n in 1:N){
                                ### recall msfile_n, pre_v_n, rt_v_n
                                msfile <- get(paste0("msfile_",n))
                                pre_v <- get(paste0("pre_v_",n))
                                preInt_v <- get(paste0("preInt_v_",n))
                                rt_v <- get(paste0("rt_v_",n))
                                
                                ### find the matched precursor 
                                rt_index <- which (rt_v %in% rt_v[abs(rt_v - t_rt) <= assign_rt_tol])
                                selected_pre <- pre_v[rt_index]
                                pre_index <- which (selected_pre %in% selected_pre[abs(selected_pre - t_pre)<= assign_ms1mz_tol])
                                
                                if (length(pre_index)== 0) {next}
                                selected_int <- preInt_v[rt_index][pre_index]
                                int_index <- which(selected_int %in% max(selected_int))[1]
                                
                                numms2 <- 1:length(msfile)
                                final_index <- numms2[rt_index][pre_index][int_index]
                                ms2mz <- msfile[[final_index]]@mz
                                ms2Int <- msfile[[final_index]]@intensity
                                
                                col <- n+ N +3
                                ms2mz_list <- paste0(ms2mz, collapse = ";")
                                ms2Int_list <- paste0(ms2Int, collapse = ";")
                                ft[i, col] <- ifelse(ms2mz_list == "", NA, ms2mz_list)
                                ft[i, col+N] <- ifelse(ms2Int_list == "", NA, ms2Int_list)
                        }
                }
        }
}
# export feature table with ms2 spectra assigned
write.csv(ft, "feature_table_ms2_assigned.csv", row.names = F)

# Remove isotope in MS2 spectra 
for (i in 1:nrow(ft)){
        for (n in 1:N){
                if (is.na(ft[i, 3+N + n])) {next}
                ms2_mz  <- as.numeric(strsplit(ft[i, 3+N + n], ";")[[1]])
                ms2_int <- as.numeric(strsplit(ft[i, 3+N+N + n], ";")[[1]])
                after <- remove_iso(ms2_mz, ms2_int)
                ft[i, 3+N + n] <- after[1]
                ft[i, 3+N+N + n] <- after[2]
        }
}
# Remove the low abundant fragments in MS2 spectra 
for (i in 1:nrow(ft)) {
        col <- 3+N
        preInt_v <- vector()
        for (j in 4:col) {
                preInt_v <- c(preInt_v, ft[i,j])
        }
        index <- order(preInt_v, decreasing = T)[1]
        col = 3+N+index
        ms2Mz <- ft[i,col]
        ms2Int <- ft[i,col+N]
        ms2Mz <- as.numeric(strsplit(ms2Mz,";")[[1]])
        ms2Int <- as.numeric(strsplit(ms2Int,";")[[1]])
        maxInt <- max(ms2Int)
        ms2Int <- ms2Int/max(ms2Int)
        # plot(x= ms2Mz , y=ms2Int, type = "h", main = paste0("ms2 for",ft[i,2]))
        if (length(ms2Mz) == 1) {
                next
        }
        remove_index <- vector()
        for (k in 1:length(ms2Int)) {
                if (ms2Int[k] <= int_threshold) {
                        remove_index <- c(remove_index, k)
                }
        }
        
        if ( length(remove_index) == 0) {
                next
        } 
        ms2Mz <- ms2Mz[-remove_index]
        ms2Int <- ms2Int[-remove_index]
        ft[i,col] <- paste0(ms2Mz, collapse = ";")
        ms2Int <- ms2Int * maxInt
        ft[i,col+N] <- paste0(ms2Int,collapse = ";")
}
#------------------------------------ factor calculation for each fragment ----------------------------------------
# Create the table to save ratio RSD, appearance rate and relative intensity
factor_tb <- as.data.frame(matrix(ncol=6+3*N))
colnames(factor_tb)[1:6] <- c("ID" , "RT", "pre_mz","original_pre_int" ,"reference_peak_mz","frag_mz")
for (i in 1:N) {
        colnames(factor_tb)[6+i] <- paste0("reference_peak_int",i)
        colnames(factor_tb)[6+N+i] <- paste0("frag_int",i)
        colnames(factor_tb)[6+2*N+i] <- paste0("ratio",i)
}
individual <- factor_tb # used to save ratio information for one fragment

# calculation process
for (i in 1:nrow(ft)){

        t_ID <- ft[i, 1]
        t_pre <- ft[i, 2]
        t_rt <- ft[i, 3]
        
        # find the file with highest peak intensity/ concentration
        preInt_v <- vector()
        end <- 3+N
        for (j in 4:end) {
                preInt_v <- c(preInt_v, ft[i,j])
        }
        t <- order(preInt_v, decreasing = T)[1]
        
        ## if less than 3 mass spectra, skip this feature####
        ms2_spectra <- vector()
        for ( j in 1:N){
                ms2_spectra <- c(ms2_spectra, ft[i, 3+N+j])
        }
        num_spectra <- ms2_spectra[!is.na(ms2_spectra)]
        
        if(length(num_spectra) < mim_num_spectra) {next}
        # no MS2 spectra information for the highest concentration, then ignore this metabolic feature
        if(is.na(ft[i,3+N + t])){next} 
        
        # record the ID, pre_mz, RT, pre_int_original
        original_pre_int <- ft[i,3+t]
        individual$ID <- t_ID
        individual$pre_mz <- t_pre
        individual$RT <- t_rt
        individual$original_pre_int <- original_pre_int
        
        # fragment information of MS2 spectra with the highest precursor concentration
        frag_v <- as.numeric(strsplit(ft[i,3+N + t], ";")[[1]])
        fragInt_v <- as.numeric(strsplit(ft[i,3+N+N+ t],";")[[1]])
        
        # if only one fragment, skip to next feature
        if (length(frag_v) ==1) {next}
        
        # found the index of most intense fragment ion and intensity
        most_intense_index <- which(fragInt_v %in% max(fragInt_v))
        reference_peak_int <- fragInt_v[most_intense_index][1]
        referencePeak <- frag_v[most_intense_index][1]
        individual$reference_peak_mz <- referencePeak
        individual$reference_peak_int1 <- reference_peak_int  # record pre_int to factor_tb
        # found a list of frag
        t_frag_v <- frag_v[-most_intense_index[1]]
        t_fragInt_v <- fragInt_v[-most_intense_index[1]]
        
        # found the ratio for each frag
        for (n in 1:length(t_frag_v)) {
                t_frag <- t_frag_v[n]
                individual$frag_mz <- t_frag
                individual$frag_int1 <- t_fragInt_v[n]
                ratio1 <- t_fragInt_v[n]/reference_peak_int
                individual$ratio1 <- ratio1
                
                # found the ratio in each file
                for (k in 2:N) {
                        
                        t1 <- order(preInt_v, decreasing = T)[k]
                        #if there is no ms2 spectra, assign NA to ratio, peak intensity
                        if (is.na(ft[i, 3+ N + t1])) {
                                assign(paste0("indi_referencePeak_int", k), NA)
                                assign(paste0("indi_frag_int", k), NA)
                                assign(paste0("ratio",k),NA)
                                next
                        }
                        
                        frag_v1 <- as.numeric(strsplit(ft[i, 3+ N + t1],";")[[1]])
                        fragInt_v1 <- as.numeric(strsplit(ft[i,3+N+N + t1],";")[[1]])
                        
                        # found the referencePeak_index and referencePeak_Int1
                        referencePeak_index1 <- which(frag_v1 %in% frag_v1[abs(frag_v1 - referencePeak) <= ratio_ms2mz_tol])
                        
                        # if no reference peak found, break, skip to next MS2 of this precursor 
                        if (length(referencePeak_index1)==0) {
                                assign(paste0("indi_referencePeak_int", k), NA)
                                assign(paste0("indi_frag_int", k), NA)
                                assign(paste0("ratio",k),NA)
                                next
                        }
                        
                        referencePeak_index2 <- which(fragInt_v1 %in% max(fragInt_v1[referencePeak_index1]))[1]
                        assign(paste0("indi_referencePeak_int", k), fragInt_v1[referencePeak_index2])
                        referencePeak_Int1 <- fragInt_v1[referencePeak_index2]
                        
                        # remove pre in frag_v1  and fragInt_v1
                        frag_v1 <- frag_v1[-referencePeak_index2]
                        fragInt_v1 <- fragInt_v1[-referencePeak_index2]
                        
                        #found the target frag_mz's intensity
                        t_frag_index <- which (frag_v1 %in% frag_v1[abs(frag_v1 - t_frag)<= ratio_ms2mz_tol])
                        ## if not found, ratio equal to zero
                        if (length(t_frag_index) == 0) {
                                assign(paste0("indi_frag_int", k),0)
                                assign(paste0("ratio",k), 0)
                        } else{
                                t_frag_index <- which(fragInt_v1 %in% max(fragInt_v1[t_frag_index]))
                                assign(paste0("indi_frag_int", k), fragInt_v1[t_frag_index[1]])
                                assign(paste0("ratio",k), fragInt_v1[t_frag_index[1]]/referencePeak_Int1)
                        }
                }
                for (l in 2:N) {
                        individual[1,6+l] <- get(paste0("indi_referencePeak_int",l))
                }
                
                for (l in 2:N) {
                        individual[1,6+N+l] <- get(paste0("indi_frag_int",l))
                }
                
                for (l in 1:N){
                        individual[1,6+2*N+ l] <- round(get(paste0("ratio",l)),4)
                }
                
                factor_tb <- rbind(factor_tb, individual)
        }
}
factor_tb <- factor_tb[-1,]

# Summarize the information of fragment ions: average Ratio, ratioSD, RSD, appearance_rate, relative intensity
factor_tb$averageRatio <- NA
factor_tb$ratioSD <- NA
factor_tb$RSD <- NA
factor_tb$appearance_rate <- NA
for (i in 1:nrow(factor_tb)) {
        ratio_v <- vector()
        ratio_noZero <- vector()
        zero <- 0
        for (j in 1:N) {
                if (!is.na(factor_tb[i,6+ 2*N + j])) {
                        ratio_v <- c(ratio_v, factor_tb[i,6+ 2*N + j])
                        if (factor_tb[i,6+ 2*N + j] == 0) {
                                zero <- zero +1
                        } else {
                                ratio_noZero <- c(ratio_noZero,factor_tb[i,6+ 2*N + j] )
                        }
                }
        }
        
        factor_tb$averageRatio[i] <- round(mean(ratio_noZero),2)
        factor_tb$ratioSD[i] <- round(sd(ratio_noZero),2)
        factor_tb$RSD[i] <- round(sd(ratio_noZero)/mean(ratio_noZero),2)
        factor_tb$appearance_rate[i] <- round((length(ratio_noZero))/length(ratio_v),2)
        if(length(ratio_v) == 1) {
                factor_tb$RSD[i] <- "One_MS2"
        }
}
factor_tb_more_MS2 <- factor_tb[factor_tb$RSD != "One_MS2",]
factor_tb_noneNA <- factor_tb_more_MS2[!is.na(factor_tb_more_MS2$RSD),]
factor_tb_noneNA$RSD <- as.numeric(factor_tb_noneNA$RSD)
#write.csv(factor_tb_noneNA, "factor_tb_for_all.csv", row.names = F)
factor_tb_NA <- factor_tb[is.na(factor_tb$RSD),]
#------------------------------------ Predict true or false for each fragment ----------------------------------------
tb <- factor_tb_noneNA
n_RSD <- 3*N +6 +3
n_meanRatio <- 2*N +6 + 1
n_appearance_rate <- 3*N +6 + 4
colnames(tb)[n_RSD]
colnames(tb)[n_meanRatio]
colnames(tb)[n_appearance_rate]
sub_tb <- tb[,c(n_RSD, n_meanRatio, n_appearance_rate)]
colnames(sub_tb) <- c("RSD","ratio1","percent")
sub_tb <- as.matrix(sub_tb)

# load model 
model <- readRDS("XGBoost_model.RDS")
message("XGBoost model has been loaded :) ")
pre <- predict(model,sub_tb)
tb$prediction <- pre
tb
for (i in 1:nrow(tb)) {
        if (tb$prediction[i] >= 0.5) {
                tb$prediction[i] <- 1
        }else {
                tb$prediction[i] <- 0  
        }
}
outputcol <- 1:6
RSD_col <- which(colnames(tb) == "RSD")
outputcol <- c(outputcol, RSD_col)
appearance_rate_col <- which(colnames(tb) == "appearance_rate")
outputcol <- c(outputcol, appearance_rate_col)
ratio_col_num <- which(colnames(tb) == "ratio1")
outputcol <- c(outputcol, ratio_col_num)
prediction_col <- which(colnames(tb) == "prediction")
outputcol <- c(outputcol, prediction_col)
prediction_output <- tb[,outputcol]
colnames(prediction_output) <- c("ID", "RT", "precursor_mz", "precursor_int","reference_peak_mz","fragment_mz",
                         "RSD","appearance_rate","relative_intensity","prediction")
write.csv(prediction_output, "prediction_output.csv", row.names = F)

#------------------------------------ Output the original and purified MS2 spectra----------------------------------------
record <- as.data.frame(matrix(ncol = 7+N))
nameN <- 3+N
colnames(record) <- c("ID", "preMz","RT",colnames(ft)[4:nameN],
                      "original_ms2_mz","original_ms2_int",
                      "purified_ms2_mz","purified_ms2_int")
sub_record <- record
ftID <- unique(tb$ID)
for ( i in 1:length(ftID)) {
        ID <- ftID[i]
        sub_tb <- tb[tb$ID==ID,]
        ft_sub <- ft[ft$ID == ID,]
        # purified
        purified_ms2_mz <- vector()
        purified_ms2_int <- vector()
        purified_ms2_mz <- c(purified_ms2_mz, sub_tb$reference_peak_mz[1])
        purified_ms2_int <- c(purified_ms2_int, sub_tb$reference_peak_int1[1])
        
        for (j in 1:nrow(sub_tb)) {
                if (sub_tb$prediction[j] == 1) {
                        purified_ms2_mz <- c(purified_ms2_mz, sub_tb$frag_mz[j])
                        purified_ms2_int <- c(purified_ms2_int, sub_tb$frag_int1[j])
                }
        }
        sub_record$ID <- ID
        sub_record$RT <- sub_tb$RT[1]
        sub_record$preMz <- sub_tb$pre_mz[1]
        
        for (j in 1:N){
                sub_record[,3+j] <- ft_sub[,3+j]
        }
        sub_record$purified_ms2_mz <- paste0(purified_ms2_mz, collapse = ";")
        sub_record$purified_ms2_int <- paste0(purified_ms2_int, collapse = ";")
        record <- rbind(record, sub_record)
        
} 
record <- record[-1,]
naFragment <- factor_tb_NA
for (i in 1:nrow(naFragment)) {
        
        t_ID <- naFragment$ID[i]
        naMz <- naFragment$frag_mz[i]
        naInt <- naFragment$frag_int1[i]
        
        index <- which(record$ID %in% t_ID)
        
        if (length(index) ==0) {
                next
        }
        Mz <- strsplit(record$purified_ms2_mz[index], ";")[[1]]
        Int <- strsplit(record$purified_ms2_int[index],";")[[1]]
        Mz <- c(Mz,naMz)
        Int <- c(Int, naInt)
        
        record$purified_ms2_mz[index] <- paste0(Mz, collapse = ";")
        record$purified_ms2_int[index] <- paste0(Int, collapse = ";")
        record$pNum_frag[index] <- record$pNum_frag[index] + 1
}
# original MS2 spectra
for (i in 1:nrow(record)) {
        ft_sub <- ft[ft[,1] == record$ID[i],]
        # add the intensity information back 
        int_v <- vector()
        for (n in 1:N) {
                record[i, 3+n] <- ft_sub[1, 3+n]
                int_v <- c(int_v, ft_sub[1, 3+n])
        }
        # add back experimental MS2 spectra
        index <- order(int_v, decreasing = T)[1]
        original_ms2_mz <- ft_sub[1, 3+N+index]
        record$original_ms2_mz[i] <- original_ms2_mz
        original_ms2_int <- ft_sub[1, 3+N+N+index]
        record$original_ms2_int[i] <- original_ms2_int
}
write.csv(record, "feature_table_ms2_purified.csv",row.names = F)
# write MGF files
df <- record
output_name <- "purified_ms2.mgf"
sink(output_name)
for(i in 1:nrow(df)){
        cat("BEGIN IONS")
        cat("\n")
        cat("TITLE=",df$ID[i],sep='')
        cat("\n")
        cat("RTINSECONDS=",df$RT[i],sep='')
        cat("\n")
        cat("PEPMASS=",df$preMz[i],sep='')
        cat("\n")
        cat("CHARGE=1+")
        cat("\n")
        mz <- as.numeric(strsplit(df$purified_ms2_mz[i],';')[[1]])
        int <- as.numeric(strsplit(df$purified_ms2_int[i],';')[[1]])
        for(j in 1:length(mz)){
                cat(mz[j], int[j])
                cat("\n")
        }
        cat("END IONS")
        cat("\n")
        cat("\n")
}
sink()
message(Sys.time())

