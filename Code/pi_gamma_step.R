pi_gamma_step<-function(gamma,K,alpha_gamma){
  n = sapply(1:K, function(k) sum(gamma==k))
  pi_gamma = rdirichlet(1,n+alpha_gamma)
  return(pi_gamma)
}