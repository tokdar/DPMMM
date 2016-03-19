sigma2_gamma_step<-function(eta_bar, gamma, m_gamma, r_gamma, s_gamma){
  n = sapply(1:K, function(k) sum(gamma==k))
  K = length(m_gamma)
  r.tilde = r_gamma + n/2
  s.tilde = rep(NA,K)
  for(k in 1:K) s.tilde[k] = s_gamma + .5*sum((gamma==k)*(eta_bar - m_gamma[k])^2 )
  sigma2_gamma = rinvgamma(K,shape=r.tilde,scale=s.tilde)
  return(sigma2_gamma)
}