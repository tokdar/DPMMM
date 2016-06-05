switch_counter <- function(my_seq, band_width){
  # calculate sequence mean
  seq_mean = mean(my_seq)
  band_up = seq_mean + band_width
  band_low = seq_mean - band_width
  # classify all points as high, low, or mid
  high_points = my_seq >= band_up
  low_points = my_seq <= band_low
  mid_points = my_seq < band_up & my_seq > band_low
  # code high points as 1, low points as -1, and mid points as -1
  point_pos = high_points - low_points
  # number of switches found
  switches = 0
  # whether we have observed a high, mid, or low point
  reached_high = FALSE
  #in_mid = FALSE
  reached_low = FALSE
  for (i in 1:40){
    current_pos = point_pos[i]
    if (current_pos == 1){
      reached_high = TRUE
    } else if (current_pos == -1){
      reached_low = TRUE
    }
    if (reached_high & reached_low){
      switches = switches + 1
      reached_high = FALSE
      # in_mid = FALSE
      reached_low = FALSE
    }
  }
  return(switches)
}

switch_wrapper <- function(my_seq, widthes){
  switch_counts = sapply(widthes, function(x) {switch_counter(my_seq, x)})
  return(switch_counts)
}

switch_matrix <- function(seq_mat, widthes){
  # this function returns a matrix
  # number of rows is the length of widthes
  # number of columns is the number of rows of seq_mat
  count_mat = apply(seq_mat, 1, function(x) {switch_wrapper(x, widthes)})
  return(count_mat)
}

make_switch_df <- function(list_matrices, widthes, n_to_sample, N){
  sampled_rows = sample(1:N, n_to_sample)
  list_mat = lapply(list_matrices, function(x) {switch_matrix(x[sampled_rows,], widthes)})
  single_mat = do.call(cbind, list_mat)
  final_df = data.frame(radius = widthes, mean = apply(single_mat, 1, mean))
  return(final_df)
}