m_gamma_step<-function(gamma,eta_bar,K,sigma2_gamma,m_0, sigma2_0){
  n = sapply(1:K, function(k) sum(gamma==k))
  sigma2_tilde = (n/sigma2_gamma + 1/sigma2_0)^(-1)
  m_tilde = sigma2_tilde*(sapply(1:K, function(k) sum(eta_bar*(gamma==k) ) )/sigma2_gamma + m_0/sigma2_0) 
  m_gamma = rnorm(K, m_tilde, sqrt(sigma2_tilde))
  return(m_gamma)
}