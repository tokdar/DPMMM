gamma_step<-function(eta_bar,m_gamma,sigma2_gamma,pi_gamma){
  nRep = length(eta_bar)
  gamma = rep(NA,nRep)
  for(i in 1:nRep){
    prob.gamma = pi_gamma*exp(-1/(2*sigma2_gamma)*(eta_bar[i]-m_gamma)^2)
    prob.gamma = prob.gamma/sum(prob.gamma)
    gamma[i] = sample(1:K, size = 1, prob = prob.gamma)
  }
  return(gamma)
}