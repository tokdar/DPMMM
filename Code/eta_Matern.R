eta_Matern<-function(A,A_star,B,B_star,Omega, zeta, sigma2, K.Matern, ell_prior, gamma, m_gamma, sigma2_gamma){
  nRep = dim(A)[1]
  t.T = dim(A)[2]
  ones = matrix(rep(1,t.T), ncol = 1)
  kappa = (A + B_star - B) - .5*(A_star + B_star) 
  eta = eta_dynamic = matrix(nrow = nRep, ncol = t.T)
  eta_bar = rep(NA,nRep)
  ones = matrix(rep(1,t.T), ncol = 1)
  L = length(K.Matern)
  K = length(m_gamma)
  ell_index = rep(NA,nRep)
  
  C = C.inv = list()
  for(l in 1:L){
    C[[l]] = C.inv[[l]] = list()
    for(k in 1:K){
      C[[l]][[k]] = sigma2*K.Matern[[l]] + sigma2_gamma[k]*ones%*%t(ones)   
      C.inv[[l]][[k]] = solve(C[[l]][[k]])
    }
  }
  
  for(i in 1:nRep){
    k = gamma[i]
    Omega_i = diag(Omega[i,])
    kappa_i = matrix(kappa[i,], ncol = 1)
    p.L = rep(NA,L)
    for(l in 1:L){
      p.L[[l]] = log(ell_prior[l])
      p.L[l] = p.L[[l]] + -.5*log(det(C[[l]][[k]]) ) - .5*log(det(Omega_i + C.inv[[l]][[k]]) )
      p.L[[l]] = p.L[[l]] -as.numeric( .5*m_gamma[k]^2*t(ones)%*%C.inv[[l]][[k]]%*%ones )
      arg = kappa_i + C.inv[[l]][[k]]%*%ones*m_gamma[k]
      p.L[[l]] = p.L[[l]] + as.numeric( .5*t(arg)%*%solve(Omega_i + C.inv[[l]][[k]])%*%arg )
    }
    
    a = max(p.L)
    log.sum = a + log( sum(exp(p.L - a)) ) 
    p.L = p.L - log.sum
    p.L = exp(p.L)
    
    ell_index[i] = sample(1:L, 1, prob = p.L)
    ##sample ell_index 
    
    C.tilde = solve( Omega_i + C.inv[[ ell_index[i] ]][[k]] ) 
    m.tilde = C.tilde%*% (kappa_i + C.inv[[ ell_index[i] ]][[k]]%*%ones*m_gamma[k] )
    
    eta_dynamic[i,] = mvrnorm(1,m.tilde,C.tilde)
    C.l.inv = solve(sigma2*K.Matern[[ ell_index[i] ]])
    v.2 = (1/sigma2_gamma[k] + t(ones)%*%C.l.inv%*%ones )^(-1)
    m.2 = v.2*(m_gamma[k]/sigma2_gamma[k] + t(ones)%*%C.l.inv%*%eta_dynamic[i,])
    eta_bar[i] = rnorm(1,m.2,sqrt(v.2)) 
    eta[i,] = zeta[i]*eta_dynamic[i,] + (1-zeta[i])*eta_bar[i]
  }
  
  ETA = list(eta_bar, eta, eta_dynamic, ell_index )
  return(ETA)
}

