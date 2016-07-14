eta_Matern_mod <-function(A,A_star,B,B_star,Omega, zeta, sigma2, K.Matern, ell_prior, gamma, m_gamma, sigma2_gamma){
  nRep = dim(A)[1]
  t.T = dim(A)[2]
  kappa = (A + B_star - B) - .5*(A_star + B_star) 
  eta = eta_dynamic = matrix(nrow = nRep, ncol = t.T)
  eta_bar = rep(NA,nRep)
  L = length(K.Matern)
  ones = matrix(rep(1,t.T), ncol = 1)
  K = length(m_gamma)
  ell_index = rep(NA,nRep)
  temp_matrices = lapply(1:L, function(l) {sigma2[l]*K.Matern[[l]]})
  C = C.inv = C.chol = C.inv.ones = list()
  for(l in 1:L){
    C[[l]] = C.inv[[l]] = C.chol[[l]] = C.inv.ones[[l]] = list()
    for(k in 1:K){
      C[[l]][[k]] = temp_matrices[[l]] + sigma2_gamma[k]  
      C.chol[[l]][[k]] = chol(C[[l]][[k]])
      C.inv[[l]][[k]] = chol2inv(C.chol[[l]][[k]])
      C.inv.ones[[l]][[k]] = m_gamma[k]*C.inv[[l]][[k]]%*%ones
    }
  }
  
  for(i in 1:nRep){
    k = gamma[i]
    Omega_i = diag(Omega[i,])
    kappa_i = matrix(kappa[i,], ncol = 1)
    p.L = rep(NA,L)
    L_omega_c_inv_sum = lapply(1:L, function(l) {chol(Omega_i + C.inv[[l]][[k]])})
    L_omega_c_inv_det = lapply(1:L, function(l) {sum(log(diag(L_omega_c_inv_sum[[l]])))})
    L_omega_c_inv_inv = lapply(1:L, function(l) {chol2inv(L_omega_c_inv_sum[[l]])})
    for(l in 1:L){
      p.L[[l]] = log(ell_prior[l])
      p.L[l] = p.L[[l]] + - sum(log(diag(C.chol[[l]][[k]])) ) - L_omega_c_inv_det[[l]]
      p.L[[l]] = p.L[[l]] -as.numeric( .5*m_gamma[k]*t(ones)%*%C.inv.ones[[l]][[k]] )
      arg = kappa_i + C.inv.ones[[l]][[k]]
      p.L[[l]] = p.L[[l]] + as.numeric( .5*t(arg)%*%L_omega_c_inv_inv[[l]]%*%arg )
    }
    
    a = max(p.L)
    log.sum = a + log( sum(exp(p.L - a)) ) 
    p.L = p.L - log.sum
    p.L = exp(p.L)
    
    ell_index[i] = sample(1:L, 1, prob = p.L)
    ##sample ell_index 
    
    C.tilde = L_omega_c_inv_inv[[ell_index[[i]]]] 
    m.tilde = C.tilde%*% (kappa_i + C.inv.ones[[ ell_index[i] ]][[k]] )
    
    eta_dynamic[i,] = mvrnorm(1,m.tilde,C.tilde)
    C.l.inv = chol2inv(chol(sigma2[ ell_index[i] ]*K.Matern[[ ell_index[i] ]]))
    t.ones.C.l.inv = t(ones) %*% C.l.inv
    v.2 = (1/sigma2_gamma[k] + t.ones.C.l.inv%*%ones )^(-1)
    m.2 = v.2*(m_gamma[k]/sigma2_gamma[k] + t.ones.C.l.inv%*%eta_dynamic[i,])
    eta_bar[i] = rnorm(1,m.2,sqrt(v.2)) 
    # this line is unncessary
    # zeta is fixed to be 1
    # eta[i,] = zeta[i]*eta_dynamic[i,] + (1-zeta[i])*eta_bar[i]
    eta[i,] = eta_dynamic[i,]
  }
  
  ETA = list(eta_bar, eta, eta_dynamic, ell_index )
  return(ETA)
}
