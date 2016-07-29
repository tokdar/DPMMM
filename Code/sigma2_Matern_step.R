sigma2_Matern_step<-function(r_0, s_0, zeta, K.Matern.chol, ell_index, eta, eta_bar){
  t.T = dim(eta)[2]
  ones = matrix(rep(1,t.T), ncol = 1)
  r.tilde <- rep(r_0, L); s.tilde <- rep(s_0, L); sigma2 <- rep(NA, L)
  for(l in 1:L){
    ix <- which(ell_index == l)
    if(length(ix) > 0){
      r.tilde[l] <- r_0 + length(ix) * t.T/2
      s.tilde[l] <- s_0 + sum(sapply(ix, function(i) sum(c(backsolve(K.Matern.chol[[l]], eta[i,] - eta_bar[i], transpose = TRUE))^2)))/2
    }
  }
  
  #sigma2 = 1/rgamma(L,shape=r.tilde,rate=s.tilde)
  sigma2 = 1/rtruncgamma(L, low = 1/9, up = 1e4, shape = r.tilde, rate = s.tilde)
  if(any(is.nan(sigma2))) print(s.tilde)
  return(sigma2)
}
