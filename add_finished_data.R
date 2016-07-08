library(dplyr)
stored_data <- list.files("Post_Summaries_Pass/")
res <- list()
keep_names <- c("Min_Max_Min.", "Min_Max_1st.Qu.", "Min_Max_Median", "Min_Max_Mean",    
                "Min_Max_3rd.Qu.", "Min_Max_Max.", "Switch_at_0.1", "Switch_at_0.2", "Switch_at_0.3")
for (file_name in stored_data){
  full_name <- paste0("Post_Summaries_Pass/",file_name)
  load(full_name)
  pos <- which(stored_data == file_name)
  res[[pos]] <- store_df
}

# add filler NA columns so each data frame has same number
# also add a variable for number of bins for future convenience
n_bins <- sapply(res, function(x) {which(colnames(x) == "Min_Max_Min.")}) - 2
max_bin <- max(n_bins)
for (i in 1:length(res)){
  to_add <- max_bin - n_bins[i]
  if (to_add > 0){
    for (j in (n_bins[i] + 1):max_bin){
      col_name <- paste0("Mean_Bin_", j)
      res[[i]][,col_name] <- NA
    }
  }
  res[[i]][,"n_bins"] <- n_bins[i]
}

temp2 <- gregexpr("[0-9]+", stored_data)
finished_triplets <- as.numeric(unique(unlist(regmatches(stored_data, temp2))))



triplet_data <- read.csv(file = "Triplets_pass_criteria.csv")
Triplet_meta = triplet_data[order(triplet_data[,"WinPr"], decreasing=T),]

named_data <- list()
pos <- 1
for (i in finished_triplets){
  named_info <- Triplet_meta[i,1:5]
  named_data[[i]] <- cbind(named_info, res[[pos]])
  pos <- pos + 1
}

new_data <- do.call(rbind, named_data)
merged_data <- left_join(Triplet_meta, new_data)
write.csv(merged_data, file = "Triplets_pass_w_data.csv", row.names = FALSE)

