phi_step<-function(eta,eta_bar,zeta,sigma2, sigma2_phi, m_phi){
  nRep = dim(eta)[1]
  t.T = dim(eta)[2]
  tmp1 = eta[,1:(t.T-1)] - matrix(rep(eta_bar,(t.T-1)), nrow = nRep, ncol = (t.T-1), byrow=F)
  eta_sum2 = sum( apply( tmp1^2,1,sum)*zeta )
  
  tmp2 = eta[,2:t.T] - matrix(rep(eta_bar,(t.T-1)), nrow = nRep, ncol = (t.T-1), byrow=F)
  eta_sum =  sum(apply( tmp1*tmp2,1,sum)*zeta)
  
  sigma2_phi_tilde = (eta_sum2/sigma2 + 1/sigma2_phi)^(-1)
  m_phi_tilde = sigma2_phi_tilde*(eta_sum/sigma2 + m_phi/sigma2_phi)
  phi = rtruncnorm(1,mean = m_phi_tilde,sd = sqrt(sigma2_phi_tilde), a = 0, b = 1)
  return(phi)
}