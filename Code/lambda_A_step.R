lambda_A_step <-function(A,A_star,alpha, r_A, s_A){
  t.T = dim(alpha)[2]
  r_A_tilde = r_A + apply(A_star - A, 2, sum)
  s_A_tilde = s_A + apply(1-alpha, 2, sum)
  lambda_A = rgamma(t.T, shape = r_A_tilde, rate = s_A_tilde)
  return(lambda_A)
}