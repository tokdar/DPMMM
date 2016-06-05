make_summary_data <- function(MCMC_results){
  # first get the means for each trail
  t.T = ncol(MCMC_results$ALPHA[[1]])
  trials = length(MCMC_results$ALPHA)
  alpha_mean = matrix(nrow = trials, ncol = t.T)
  for(i in 1:trials){
    alpha_mean[i,] = apply(MCMC_results$ALPHA[[i]],2,mean)
  }
  # now get summary statustics for MinMax
  MinMax_summary = summary(MCMC_results$MinMax.Pred)
  MinMax_matrix = matrix(MinMax_summary, ncol = length(MinMax_summary),
                         nrow = trials, byrow = TRUE)
  
  # count the number of switches in post pred draws
  widthes = c(.05, .1, .15)
  switch_counts = switch_matrix(MCMC_results$ALPHA_pred, widthes)
  switch_mean = apply(switch_counts, 1, mean)
  switch_matrix = matrix(switch_mean, ncol = length(widthes),
                         nrow = trials, byrow = TRUE)

  mean_names = paste0("Mean_Bin_", 1:t.T)
  colnames(alpha_mean) = mean_names
  MinMax_names = paste0("Min_Max_", names(MinMax_summary))
  colnames(MinMax_matrix) = MinMax_names
  switch_names = paste0("Switch_at_", 2*widthes)
  colnames(switch_matrix) = switch_names
  joint_matrix = cbind(alpha_mean, MinMax_matrix, switch_matrix)
  result_df = data.frame(Trial = 1:trials,
                         joint_matrix)
  
  return(result_df)
}

add_IDs <- function(triplets, triplet_meta, list_result_matrices, save_name){
  # this code assumes that the cell ID is in the FIRST COLUMN
  cell_ids = triplet_meta[triplets, 1]
  cell_id = colnames(triplet_meta)[1]
  data_names = colnames(list_result_matrices[[1]])
  for (i in 1:length(list_result_matrices)){
    list_result_matrices[[i]] = cbind(cell_ids[i], list_result_matrices[[i]])
  }
  all_results = do.call(list_result_matrices, cbind)
  colnames(all_results) = c(cell_id, data_names)
  merged_data = left_join(triplet_meta, all_results)
  write.csv(merged_data, file = save_name)
}



