MCMC.triplet.raw<-function(triplet, trials, spiketimes, ell_0, ETA_BAR_PRIOR, MinMax.Prior){
  
  Triplet_fig_dir = paste(Fig_dir,"Triplet_",triplet,"/",sep="")
  unlink(Triplet_fig_dir, recursive=T)
  dir.create(Triplet_fig_dir)
  # this function translates the raw time data into counts
  # for bins of width 25ms
  # written by Surja, I believe
  # modifing so that we can run with raw local data
  # no need to fetch
  # retX = Pre_Proc(triplet, bin_width = 25, Triplet_meta)
  bin_width = 25
  AB = Bincounts(trials, spiketimes, frq = c(1,2), 
                   pos = c(24, -6), on.reward = TRUE, start.time = 0,
                   end.time = 1000, bw = bin_width, target = c(25, 150, 800), 
                   match.level = FALSE, AB.eqlevel = FALSE, go.by.soff = TRUE, 
                   n.iter = 1e3, plot = FALSE, faux.dual.mix = FALSE, 
                   faux.dual.int = FALSE, faux.alpha = 0.5, faux.dual.swi = FALSE,
                   faux.target = c(100, 100), nAB = "match.realABcount")
  X = AB[[3]]
  
  
  x1 = AB[[1]]
  x2 = AB[[2]]
  x3 = AB[[3]]
  
  t.T <- nrow(X)
  n1 <- ncol(x1)
  n2 <- ncol(x2)
  n3 <- ncol(x3)
  
  get.hyper1 <- smoogam(x1);
  get.hyper2 <- smoogam(x2);
  m.1 <- get.hyper1$mean; am.1 <- n1 * m.1; bm.1 <- rep(n1, t.T); s.1 <- sqrt(am.1)/bm.1
  m.2 <- get.hyper2$mean; am.2 <- n2 * m.2; bm.2 <- rep(n2, t.T); s.2 <- sqrt(am.2)/bm.2
  
  r_A = am.1; s_A = bm.1
  r_B = am.2; s_B = bm.2
  
  # this function returns a list
  # first is counts for each bin
  # take the transpose
  # now each row is a time series
  X = t(X)
  # then the prior parameters for A
  # these are derived from the single source A
  # they only affect the model through the hyperparameters

  
  # get the number of repetitions, which is the number of rows
  nRep = dim(X)[1]
  # number of times (bins) is the number of columns
  t.T = dim(X)[2]
  
  # L is the number of lengths we are considering (usually 5)
  L = length(ell_0)
  
  # sampling initial values?
  m_gamma = rnorm(K,m_0, sqrt(sigma2_0))
  sigma2_gamma = rinvgamma(K,r_gamma,s_gamma)
  sigma2 = rinvgamma(L,r_0,s_0)
  pi_gamma = rDirichlet(1,rep(alpha_gamma,K))
  
  # sample gammas
  gamma = sample(1:K, size = nRep, prob = pi_gamma, replace = TRUE)
  
  # sample initial lambdas
  lambda_A =  rgamma(t.T, shape = r_A, rate = s_A)
  lambda_B =  rgamma(t.T, shape = r_B, rate = s_B)
  
  
  
  # covariance matrices
  K.SE = K.SE.inv = list()
  for(l in 1:L){
    t.mat = matrix(1:t.T, nrow = t.T, ncol = t.T)
    s.mat = t(t.mat)
    K.SE[[l]] = exp(-abs(t.mat - s.mat)^2/(2*ell[l]^2)) + 1e-12*diag(t.T)
    # solve to get inverse
    K.SE.inv[[l]] = solve(K.SE[[l]])
  }
  
  ell_index = rep(1,nRep)
  ell_prior = rep(1/L,L)
  # get eta matrix
  eta = matrix(nrow = nRep, ncol = t.T)
  for(i in 1:nRep){
    C = sigma2[ ell_index[i] ]*K.SE[[ ell_index[i] ]]
    eta[i,] = mvrnorm(1, rep(0,t.T), C)
  }
  eta_dynamic = eta
  # take the average
  eta_bar = apply(eta,1,mean)
  
  
  ## THIS NEEDS TO BE REMOVED
  # depricated feature
  # set to be 1 and no changed
  zeta = rep(1,nRep)
  
  # these will hold the results
  ALPHA = A_POST = list()
  ALPHA_BAR_POST = ZETA_POST = N_switch_alpha_bar = N_mode_switch = MinMax = matrix(nrow = N.MC, ncol = nRep)
  ETA_BAR_POST = MinMax.Pred = rep(NA,N.MC)
  SIGMA2_POST = Ell_Post = matrix(nrow = N.MC, ncol = L)
  
  lambda_A_POST = lambda_B_POST = matrix(nrow = N.MC, ncol = t.T)
  
  ALPHA = lapply(1:nRep, function(x) {matrix(nrow = N.MC, ncol = t.T)})
  A_POST = ALPHA
  ALPHA_pred = matrix(nrow = N.MC, ncol = t.T) 
  ###############################
  #inside MCMC loop
  alpha = 1/(1+exp(-eta))
  mm = 0
  
  for(m in -burnIn:N.MC*thin){
    
    # Block 1
    A = A_step(X,lambda_A,lambda_B,alpha)
    B = X - A
    A_star = A_star_step(A,lambda_A,alpha)
    B_star = B_star_step(B,lambda_B,alpha)
    #Block 2
    lambda_A = lambda_A_step(A,A_star,alpha, r_A, s_A)
    lambda_B = lambda_B_step(B,B_star,alpha, r_B, s_B)
    
    #Block 3
    Omega = Omega_step(A_star, B_star,eta)
    
    #Block 4 
    # no longer updated zeta, it is fixed at 1
    #zeta = Zeta_step_collapsed(A,A_star,B,B_star,Omega,gamma,m_gamma,sigma2_gamma,sigma2, phi, p)
    ETA = eta_Matern_mod(A,A_star,B,B_star,Omega, zeta, sigma2, K.SE, ell_prior, gamma, m_gamma, sigma2_gamma)
    eta_bar = ETA[[1]]
    eta = ETA[[2]]
    eta_dynamic = ETA[[3]]
    ell_index = ETA[[4]]
    
    alpha = 1/(1+exp(-eta))
    
    ell_prior = ell_prior_step(ell_index,ell_0,L)
    #Block 6
    sigma2 = sigma2_Matern_step(r_0, s_0, zeta, K.SE.inv, ell_index, eta, eta_bar)
    
    #Block 7
    gamma = gamma_step(eta_bar,m_gamma,sigma2_gamma,pi_gamma)
    m_gamma = m_gamma_step(gamma,eta_bar,K,sigma2_gamma,m_0, sigma2_0)
    sigma2_gamma = sigma2_gamma_step(eta_bar, gamma, m_gamma, r_gamma, s_gamma)
    pi_gamma = pi_gamma_step(gamma,K,alpha_gamma)
    
    if(m>0 & (m %%thin) == 0){
      mm = mm + 1
      alpha_bar_tmp = 1/(1+exp(-apply(eta,1,mean) ) )
      # is alpha_modes used for anything??
      alpha_modes = mean(alpha_bar_tmp)
      ALPHA_BAR_POST[mm,] = alpha_bar_tmp
      # do not need ZETA_POST
      ZETA_POST[mm,] = zeta
      SIGMA2_POST[mm,] = sigma2
      # sample from posterior predictive for ETA_BAR
      s = sample(1:K,1, prob = pi_gamma)
      ETA_BAR_POST[mm] = rnorm(1, m_gamma[s], sqrt(sigma2_gamma[s]) )
      
      lambda_A_POST[m,] = lambda_A
      lambda_B_POST[m,] = lambda_B
      # draw from posterior predictive
      # for ETA and alpha
      # this is used to calculate MinMax
      Ell_Post[mm,] = ell_prior[1,]
      l.pred = sample(1:L,1,prob=ell_prior[1,])
      eta.pred = ETA_BAR_POST[mm] + mvrnorm(1,rep(0,t.T),sigma2[l.pred]*K.SE[[l.pred]])
      alpha.pred = 1/(1+exp(-eta.pred))
      ALPHA_pred[mm,] = alpha.pred
      MinMax.Pred[mm] = max(alpha.pred) - min(alpha.pred)
      
      for(i in 1:nRep){
        ALPHA[[i]][mm,] = alpha[i,]
        A_POST[[i]][mm,] = A[i,]
        MinMax[mm,] = apply(alpha,1,max) - apply(alpha,1,min)
      }
    }
    
  }
  
  
  
  ## End of MCMC loop
  ##########################################
  
  alpha_mean = alpha_pt025  = alpha_pt975 = matrix(nrow = nRep, ncol = t.T)
  A_mean = A_pt025 = A_pt975 = matrix(nrow = nRep, ncol = t.T)
  for(i in 1:nRep){
    # compute the mean for alpha
    alpha_mean[i,] = apply(ALPHA[[i]],2,mean)
    # compute 2.5% quantile
    alpha_pt025[i,] = apply(ALPHA[[i]],2,quantile,.025)
    # compute 97.5% quantile
    alpha_pt975[i,] = apply(ALPHA[[i]],2,quantile,.975)
    
    # compute mean, quantiles for A_POST
    A_mean[i,] = apply(A_POST[[i]],2,mean)
    A_pt025[i,] = apply(A_POST[[i]],2,quantile,.025)
    A_pt975[i,] = apply(A_POST[[i]],2,quantile,.975)
  }
  # compute median of MinMax
  Median.Min.Max = apply(MinMax,2,median)
  
  data_to_summarize = list(ALPHA = ALPHA, MinMax.Pred = MinMax.Pred,
                           ALPHA_pred = ALPHA_pred)
  
  store_df = make_summary_data(data_to_summarize)
  # save a data file with the results
  
  
  return_list = list(
    triplet = triplet,
    ALPHA_BAR_POST = ALPHA_BAR_POST,
    ALPHA = ALPHA,
    SIGMA2_POST = SIGMA2_POST,
    ETA_BAR_POST = ETA_BAR_POST,
    lambda_A_POST = lambda_A_POST,
    lambda_B_POST = lambda_B_POST,
    Ell_Post = Ell_Post,
    MinMax.Prior = MinMax.Prior,
    MinMax.Pred = MinMax.Pred,
    A_POST = A_POST,
    ALPHA_pred = ALPHA_pred,
    MinMax = MinMax)
  
  save(file = paste(Triplet_dir,"triplet_",triplet,".RData", sep=""),
       alpha_mean, 
       alpha_pt025,
       alpha_pt975,
       A_mean, 
       A_pt025, 
       A_pt975,
       ALPHA_BAR_POST,
       ETA_BAR_POST, 
       ALPHA_pred,
       Median.Min.Max,
       store_df,
       return_list)  
  return(return_list)
}

