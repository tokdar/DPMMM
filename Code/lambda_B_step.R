lambda_B_step <-function(B,B_star,alpha, r_B, s_B){
  t.T = dim(alpha)[2]
  r_B_tilde = r_B + apply(B_star - B, 2, sum)
  s_B_tilde = s_B + apply(alpha, 2, sum)
  lambda_B = rgamma(t.T, shape = r_B_tilde, rate = s_B_tilde)
  return(lambda_B)
}