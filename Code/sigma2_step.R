sigma2_step<-function(r_0, s_0, zeta, phi, eta, eta_bar){
  n.Z = sum(zeta)
  nRep = dim(eta)[1]
  t.T = dim(eta)[2]
  r.tilde = r_0 + n.Z*t.T/2
  eta.bar.mat = matrix(rep(eta_bar,t.T), nrow = nRep, byrow=F)
  eta.lag = cbind(eta_bar, eta[,1:(t.T-1)])
  eta.diff = eta.bar.mat + phi*(eta.lag - eta.bar.mat) 
  eta.arg = (eta - eta.diff)^2
  s.tilde = s_0 + .5*sum( apply(eta.arg,1,sum)*zeta )
  sigma2 = rinvgamma(1,shape=r.tilde,scale=s.tilde)
  return(sigma2)
}