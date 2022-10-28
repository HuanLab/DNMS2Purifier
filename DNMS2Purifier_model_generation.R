# DNMS2Purifier model generation
# This R script reads in the mzML or mzXML files of DDA data,
#               outputs csv files and mgf files containing the metabolic features table with purified ms2 spectra
# Tingting Zhao, October 3, 2022
# Copyright @ The University of British Columbia

# ------ load library ----------------------------------------------------#
library(xcms) # load xcms package
library(xgboost)

work_directory <- "D:/2021-11-05-Ratio Study/code/DNMS2Purifier_workflow"
#------- parameters setting ----------------------------------------------
assign_rt_tol <- 15 # retention time tolerance (in seconds) for MS2 assignment
assign_ms1mz_tol <- 0.01 # m/z tolerance for MS2 assignment
ratio_ms2mz_tol <- 0.01 # m/z tolerance when calculating the ratio for each fragment ion
match_ms2mz_tol <- 0.02 # m/z tolerance when labeling fragments as true or false against library spectra
int_threshold <- 0.01 #threshold to remove the fragments whose intensity is lower than this value

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
refname <- "standards_information.csv"   # data path of csv file which contains lists of true fragments for target metabolites
setwd(work_directory)
# Load the information of the target metabolites: "name", "mz", "RT" 
target_ft <- read.csv(refname)
target_ft <- target_ft[,1:3] 
message("feature table of targert chemicals has been loaded :)")
DDAfiles <- list.files(pattern = ".mzXML")
if(length(DDAfiles) == 0) {DDAfiles <- list.files(pattern = ".mzML")}
N <- length(DDAfiles)
message(paste0(N, " of mzXML raw data have been loaded :)"))
# add column to the feature table to save the ms2 spectra information 
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

# read mzXML file mass spectrum, save the precursor_v, and rt_v 
message(" start to read raw LC-MS data :)")
for (i in 1:N) {
        if (i==1){ t0 <- Sys.time()}
        msfile <- DDAfiles[[i]]
        msfile <- readMSData(msfile, msLevel. = 2, centroided. = T)

        pre_v <- vector()
        preInt_v <- vector()
        rt_v <- vector()
        for (n in 1:length(msfile)){
                pre_v <- c(pre_v,msfile[[n]]@precursorMz)
                preInt_v <- c(preInt_v, msfile[[n]]@precursorIntensity)
                rt_v <- c(rt_v, msfile[[n]]@rt)
        }
        # save information to msfile_i, pre_v, rt_v, preInt_v
        assign(paste0("msfile_",i), msfile)
        assign(paste0("pre_v_",i), pre_v)
        assign(paste0("rt_v_",i), rt_v)
        assign(paste0("preInt_v_",i), preInt_v) 
}
# assign ms2 to each feature 
message(" start to assign ms2 spectra to each feature :)")
for (i in 1:nrow(target_ft)) {
        if (i==1){t0 <- Sys.time()}
        t_pre <- target_ft$mz[i]
        t_rt <- as.numeric(target_ft$RT[i]) * 60
        ## extract ms2 information(precursor intensity, ms2 m/z, ms2 intensity) in each file
        for (n in 1:N){
                #retrieve msfile_n, pre_v_n, rt_v_n
                msfile <- get(paste0("msfile_",n))
                pre_v <- get(paste0("pre_v_",n))
                preInt_v <- get(paste0("preInt_v_",n))
                rt_v <- get(paste0("rt_v_",n))
                ### find the matched precursor 
                rt_index <- which (rt_v %in% rt_v[abs(rt_v - t_rt) <= assign_rt_tol ])
                if (length(rt_index) == 0) {next}
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
#write.csv(target_ft,"feature_table_ms2_assigned.csv", row.names = F)
#------------------------------------ Remove isotope in MS2 spectra ------------------------------------------------
ft <- target_ft
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
#------------------------------------ Remove the low abundant fragments in MS2 spectra -----------------------------
for (i in 1:nrow(ft)) {
        col <- 3+N
        preInt_v <- vector()
        for (j in 4:col) { preInt_v <- c(preInt_v, ft[i,j])}
        index <- order(preInt_v, decreasing = T)[1]
        col = 3+N+index
        ms2Mz <- ft[i,col]
        ms2Int <- ft[i,col+N]
        ms2Mz <- as.numeric(strsplit(ms2Mz,";")[[1]])
        ms2Int <- as.numeric(strsplit(ms2Int,";")[[1]])
        maxInt <- max(ms2Int)
        ms2Int <- ms2Int/max(ms2Int)
        if (length(ms2Mz) == 1) {next}
        remove_index <- vector()
        for (k in 1:length(ms2Int)) {if (ms2Int[k] <= int_threshold) {remove_index <- c(remove_index, k)}}
        if ( length(remove_index) == 0) {next} 
        ms2Mz <- ms2Mz[-remove_index]
        ms2Int <- ms2Int[-remove_index]
        ft[i,col] <- paste0(ms2Mz, collapse = ";")
        ms2Int <- ms2Int * maxInt
        ft[i,col+N] <- paste0(ms2Int,collapse = ";")
}
#----------------------- Calculate ratio RSD, appearance rate, relative intensity for each fragment --------------
# create the table to save ratio information
factor_tb <- as.data.frame(matrix(ncol=6+3*N))
colnames(factor_tb)[1:6] <- c("ID" , "RT", "pre_mz","original_pre_int" ,"reference_mz","frag_mz")
for (i in 1:N) {
        colnames(factor_tb)[6+i] <- paste0("reference_peak_int",i)
        colnames(factor_tb)[6+N+i] <- paste0("frag_int",i)
        colnames(factor_tb)[6+2*N+i] <- paste0("ratio",i)
}
individual <- factor_tb
for (i in 1:nrow(ft)){
        t_ID <- ft[i, 1]
        t_pre <- ft[i, 2]
        t_rt <- ft[i, 3]
        # find the file with highest peak intensity/concentration
        preInt_v <- vector()
        end <- 3+N
        for (j in 4:end) {preInt_v <- c(preInt_v, ft[i,j])}
        t <- order(preInt_v, decreasing = T)[1]
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
        # if less than 3 mass spectra, skip this feature####
        ms2_spectra <- vector()
        for ( j in 1:N){
                ms2_spectra <- c(ms2_spectra, ft[i, 3+N+j])
        }
        num_spectra <- ms2_spectra[!is.na(ms2_spectra)]
        if(length(num_spectra) < 3) {next}
        
        # if only one fragment, skip to next feature
        if (length(frag_v) ==1) {next}
        # found the index of most intense fragment ion and intensity
        most_intense_index <- which(fragInt_v %in% max(fragInt_v))
        reference_peak_int <- fragInt_v[most_intense_index][1]
        referencePeak <- frag_v[most_intense_index][1]
        individual$reference_mz <- referencePeak
        individual$reference_peak_int1 <- reference_peak_int  # record pre_int to factor_tb
        t_frag_v <- frag_v[-most_intense_index[1]]
        t_fragInt_v <- fragInt_v[-most_intense_index[1]]
        # calculate the ratio for each frag
        for (n in 1:length(t_frag_v)) {
                t_frag <- t_frag_v[n]
                individual$frag_mz <- t_frag
                individual$frag_int1 <- t_fragInt_v[n]
                ratio1 <- t_fragInt_v[n]/reference_peak_int
                individual$ratio1 <- ratio1
                for (k in 2:N) {
                        t1 <- order(preInt_v, decreasing = T)[k]
                        #if there is no ms2 spectra, assign NA to ratio, peak intensity
                        if (is.na(ft[i, 3+ N + t1])) {
                                assign(paste0("indi_referencePeak_int", k), NA)
                                assign(paste0("indi_frag_int", k), NA)
                                assign(paste0("ratio",k),NA)
                                next}
                        frag_v1 <- as.numeric(strsplit(ft[i, 3+ N + t1],";")[[1]])
                        fragInt_v1 <- as.numeric(strsplit(ft[i,3+N+N + t1],";")[[1]])
                        # found the reference Peak_index and reference Peak_Int1
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
                for (l in 2:N) {individual[1,6+l] <- get(paste0("indi_referencePeak_int",l))}
                for (l in 2:N) {individual[1,6+N+l] <- get(paste0("indi_frag_int",l))}
                for (l in 1:N){individual[1,6+2*N+ l] <- round(get(paste0("ratio",l)),4)}
                factor_tb <- rbind(factor_tb, individual)
        }
}
factor_tb <- factor_tb[-1,]
factor_tb$RSD <- NA
factor_tb$percent <- NA
for (i in 1:nrow(factor_tb)) {
        ratio_v <- vector()
        ratio_noZero <- vector()
        zero <- 0
        for (j in 1:N) {
                if (!is.na(factor_tb[i,6+ 2*N + j])) {
                        ratio_v <- c(ratio_v, factor_tb[i,6+ 2*N + j])
                        if (factor_tb[i,6+ 2*N + j] == 0) {zero <- zero +1} else {
                                ratio_noZero <- c(ratio_noZero,factor_tb[i,6+ 2*N + j] )}
                }
        }
        factor_tb$RSD[i] <- round(sd(ratio_noZero)/mean(ratio_noZero),2)
        factor_tb$percent[i] <- round((length(ratio_noZero))/length(ratio_v),2)
}
factor_tb<- factor_tb[!is.na(factor_tb$RSD),]

#------------------------------------ Label each fragment by referring to the library MS2 spectra -------------------
# load true fragment information of target metabolites 
ref <- read.csv(refname)
# Label each fragment and save labeled RSD information into a data frame
labeled_tb <- as.data.frame(matrix(ncol=8+3*N+2))
colnames(labeled_tb) <- c(colnames(factor_tb), "false_or_true","Metabolite")

for (i in 1:nrow(ref)) {
        t_name <- ref$name[i]
        t_mz <- as.numeric(ref$mz[i])
        t_rt <- as.numeric(ref$RT[i])*60
        ref_frag_list <- as.numeric(strsplit(ref$frag[i], ";")[[1]])
        matched_factor_tb <- factor_tb[abs(factor_tb$pre_mz - t_mz) <=0.01,]
        if (nrow(matched_factor_tb) ==0 ) {next}
        matched_factor_tb <- matched_factor_tb[abs(matched_factor_tb$RT*60 - t_rt)<= assign_rt_tol,]
        if (nrow(matched_factor_tb) ==0 ) {next}
        matched_factor_tb$false_or_true <- NA
        for (n in 1:nrow(matched_factor_tb)) {
                exp_frag <- matched_factor_tb$frag_mz[n]
                if (min(abs(ref_frag_list - exp_frag)) <= match_ms2mz_tol) {
                        matched_factor_tb$false_or_true[n] <- 1
                } else{ matched_factor_tb$false_or_true[n] <- 0 }
        }
        matched_factor_tb$Metabolite <- t_name
        labeled_tb <- rbind(labeled_tb, matched_factor_tb)
        
}
labeled_tb <- labeled_tb[-1,]
#write.csv(labeled_tb, "ground_truth_data.csv", row.names = F)
#-------------- please conduct some manual inspection for the ground truth data -------------
#----------------------------- model training -----------------------------------------------
require(xgboost) # load the package for xgboost machine learning  model
set.seed(222)
data <- labeled_tb
trainData <- subset(data, select = c("RSD", "ratio1","percent","false_or_true"))
# XGBoost
bst <- xgboost(data = as.matrix(trainData[,-4]), 
               label = as.numeric(as.character(trainData$false_or_true)),
               max.depth = 4, 
               eta = 1, 
               nthread = 2,
               nrounds = 1500, 
               objective = "binary:logistic")
saveRDS(bst, "XGBoost_model.RDS")