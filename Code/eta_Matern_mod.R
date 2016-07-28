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
    k <- gamma[i]
    Omega_i <- diag(Omega[i,], t.T)
    #print(Omega[i,])
    Omega.inv_i <- diag(1/Omega[i,], t.T)
    kappa_i <- kappa[i,]
    z_i <- kappa[i,] / Omega[i,]

    chol_C_plus_Omega.inv <- lapply(1:L, function(l) return(chol(C[[l]][[k]] + Omega.inv_i)))
    log_p.L <- sapply(chol_C_plus_Omega.inv, function(R) return(-sum(log(diag(R))) - 0.5*sum(c(backsolve(R, z_i - m_gamma[k], transpose = TRUE))^2)))
    #print(log_p.L)
    #print(z_i - m_gamma[k])
    log_p.L <- log_p.L + log(ell_prior)
  
    p.L <- exp(log_p.L - max(log_p.L))
    p.L <- p.L / sum(p.L)
    ell_index[i] <- sample(L, 1, prob = p.L)

    C.tilde.inv <- chol2inv(chol(C[[ell_index[i]]][[k]])) + Omega_i
    R.tilde <- chol(C.tilde.inv)
    m.tilde <- backsolve(R.tilde, backsolve(R.tilde, kappa[i,] + m_gamma[k]*c(C.inv[[ell_index[i]]][[k]]%*%as.matrix(rep(1,t.T))), transpose = TRUE))
            
    eta_dynamic[i,] <- m.tilde + c(backsolve(R.tilde, rnorm(t.T)))
    
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
