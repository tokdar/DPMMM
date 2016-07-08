library(dplyr)
stored_data <- list.files("Post_Summaries_Pass/")
res <- list()
keep_names <- c("Min_Max_Min.", "Min_Max_1st.Qu.", "Min_Max_Median", "Min_Max_Mean",    
                "Min_Max_3rd.Qu.", "Min_Max_Max.", "Switch_at_0.1", "Switch_at_0.2", "Switch_at_0.3")
for (file_name in stored_data){
  full_name <- paste0("Post_Summaries_Pass/",file_name)
  load(full_name)
  pos <- which(stored_data == file_name)
  res[[pos]] <- store_df[1,keep_names]
}
data_matrix <- do.call(rbind, res)

percn9_med <- quantile(data_matrix$Min_Max_Median, c(.9))
high_med <- which(data_matrix$Min_Max_Median > percn9_med)

percn9_mean <- quantile(data_matrix$Min_Max_Mean, c(.9))
high_mean <- which(data_matrix$Min_Max_Mean > percn9_mean)

percn9_2switch <- quantile(data_matrix$Switch_at_0.2, c(.9))
high_switch <- which(data_matrix$Switch_at_0.2 > percn9_2switch)

switch_triplets <- intersect(intersect(high_med, high_mean), high_switch)

percn1_med <- quantile(data_matrix$Min_Max_Median, .1)
low_med <- which(data_matrix$Min_Max_Median < percn1_med)

percn1_mean <- quantile(data_matrix$Min_Max_Mean, .1)
low_mean <- which(data_matrix$Min_Max_Mean < percn1_mean)

percn1_switch <- quantile(data_matrix$Switch_at_0.2, .1)
low_switch <- which(data_matrix$Switch_at_0.2 < percn1_switch)

stable_triplets <- intersect(intersect(low_med, low_mean), low_switch)


temp2 <- gregexpr("[0-9]+", stored_data)
finished_triplets <- as.numeric(unique(unlist(regmatches(stored_data, temp2))))
triplet_data <- read.csv(file = "Triplets_pass_criteria.csv")
Triplet_meta = triplet_data[order(triplet_data[,"WinPr"], decreasing=T),]

switch_rows <- finished_triplets[switch_triplets]
stable_rows <- finished_triplets[stable_triplets]

keep_rows <- c(switch_rows, stable_rows)
Triplet_meta$Switch <- "Switch"
Triplet_meta$Switch[stable_rows] <- "Stable"
keep_cols <- c(1:11, which(colnames(Triplet_meta) == "Switch"))
Extreme_Triplets <- Triplet_meta[keep_rows,keep_cols]

urls <- paste0(Extreme_Triplets$Cell, "_",
               "Site", Extreme_Triplets$Site,
               "_Freq", Extreme_Triplets$AltFreq,
               "_Pos", Extreme_Triplets$AltPos,
               ".png")
base_url <- "https://github.com/azz2/DPMMM/tree/master/summary_images"
urls <- paste0(base_url, "/", urls)

Extreme_Triplets$url <- urls

write.csv(Extreme_Triplets, file = "Switch_Stable_Triplets.csv", row.names = FALSE)
