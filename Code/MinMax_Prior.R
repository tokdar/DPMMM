K.SE = K.SE.inv = list()
t.T = 40
for(l in 1:L){
  K.SE[[l]] = matrix(nrow=t.T,ncol=t.T)
  for(t in 1:t.T){
    for(s in 1:t.T){
      r = abs(t-s)
      K.SE[[l]][t,s] = exp(-r^2/(2*ell[l]^2) ) # Squared exponential
    }
  }
  K.SE[[l]] = K.SE[[l]] + .01*diag(40) 
  K.SE.inv[[l]] = solve(K.SE[[l]])
}

eta = matrix(nrow = nSamples,ncol = t.T)
for(i in 1:nSamples){
  pi_gamma = rDirichlet(1,rep(1,K)/K)
  m_gamma = rnorm(K,0,1)
  sigma2_gamma = rinvgamma(K,shape = r_gamma,scale = s_gamma)
  k = sample(1:K,size=1,prob=pi_gamma)
  sigma2 = rinvgamma(1,shape = r_0, scale = s_0)
  
  eta_bar = rnorm(1,m_gamma[k], sqrt(sigma2_gamma[k]))
  
  rho = rDirichlet(1, ell_0 )
  ell_index = sample(1:L, 1, prob = rho)
  C = sigma2*K.SE[[ ell_index ]]
  eta[i,] = eta_bar + mvrnorm(1,rep(0,t.T),C)
}

alpha = 1/(1+exp(-eta))
MinMax.Prior = apply(alpha,1,max) - apply(alpha,1,min)


