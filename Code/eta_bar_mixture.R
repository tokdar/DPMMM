nSamples = N.MC
PI = rDirichlet(nSamples,rep(alpha_gamma,K))

sigma2_gamma = matrix( rinvgamma(K*nSamples,shape = r_gamma,scale = s_gamma), nrow = nSamples)
m_gamma = matrix( rnorm(K*nSamples,m_0, sqrt(sigma2_0)), nrow = nSamples)
gamma = eta_bar = rep(NA,nSamples)

for(i in 1:nSamples){
  gamma[i] = sample(1:K, 1, prob = PI[i,])
  eta_bar[i] = rnorm(1, m_gamma[i,gamma[i]], sqrt(sigma2_gamma[i,gamma[i]]))
}

ETA_BAR_PRIOR = eta_bar

