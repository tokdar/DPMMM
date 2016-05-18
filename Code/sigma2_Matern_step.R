sigma2_Matern_step<-function(r_0, s_0, zeta, K.Matern.inv, ell_index, eta, eta_bar){
  n.Z = sum(zeta)
  nRep = dim(eta)[1]
  t.T = dim(eta)[2]
  r.tilde = r_0 + n.Z*t.T/2
  ones = matrix(rep(1,t.T), ncol = 1)
  s.tilde = s_0 + .5*sum( sapply(1:nRep, function(i) t(eta[i,] - eta_bar[i]*ones)%*%K.Matern.inv[[ ell_index[i] ]]%*%(eta[i,] - eta_bar[i]*ones) ) )
  sigma2 = rinvgamma(L,shape=r.tilde,scale=s.tilde)
  return(sigma2)
}