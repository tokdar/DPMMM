MCMC.triplet<-function(triplet, ell_0, ETA_BAR_PRIOR, MinMax.Prior){
  
  Triplet_fig_dir = paste(Fig_dir,"Triplet_",triplet,"/",sep="")
  unlink(Triplet_fig_dir, recursive=T)
  dir.create(Triplet_fig_dir)
  # this function translates the raw time data into counts
  # for bins of width 25ms
  # written by Surja, I believe
  retX = Pre_Proc(triplet, bin_width = 25, Triplet_meta)
  # this function returns a list
  # first is counts for each bin
  X = retX[[1]]
  # then the prior parameters for A
  # these are derived from the single source A
  # they only affect the model through the hyperparameters
  r_A = retX[[2]]
  s_A = retX[[3]]
  # hyperparameters for B
  r_B = retX[[4]]
  s_B = retX[[5]]
  # take the transpose
  # now each row is a time series
  X = t(X)
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
  gamma = rep(NA,nRep)
  
  for(i in 1:nRep){
    k = sample(1:K,size=1,prob=pi_gamma)
    gamma[i] = k
  }
  # sample initial lambdas
  lambda_A =  rgamma(t.T, shape = r_A, rate = s_A)
  lambda_B =  rgamma(t.T, shape = r_B, rate = s_B)
  
  
  
  # covariance matrices
  K.SE = K.SE.inv = list()
  for(l in 1:L){
    K.SE[[l]] = matrix(nrow=t.T,ncol=t.T)
    for(t in 1:t.T){
      for(s in 1:t.T){
        r = abs(t-s)
        K.SE[[l]][t,s] = exp(-r^2/(2*ell[l]^2) ) # Squared exponential
      }
    }
    # add on diagonal
    K.SE[[l]] = K.SE[[l]] + .01*diag(t.T) 
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
  for(i in 1:nRep){
    ALPHA[[i]] = A_POST[[i]] = matrix( nrow = N.MC, ncol = t.T)
  }
  
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
    
    ETA = eta_Matern(A,A_star,B,B_star,Omega, zeta, sigma2, K.SE, ell_prior, gamma, m_gamma, sigma2_gamma)
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
      alpha_modes = mean(alpha_bar_tmp)
      ALPHA_BAR_POST[mm,] = alpha_bar_tmp
      ZETA_POST[mm,] = zeta
      SIGMA2_POST[mm,] = sigma2
      s = sample(1:K,1, prob = pi_gamma)
      ETA_BAR_POST[mm] = rnorm(1, m_gamma[s], sqrt(sigma2_gamma[s]) )
      
      lambda_A_POST[m,] = lambda_A
      lambda_B_POST[m,] = lambda_B

      Ell_Post[mm,] = ell_prior[1,]
      l.pred = sample(1:L,1,prob=ell_prior[1,])
      eta.pred = ETA_BAR_POST[mm] + mvrnorm(1,rep(0,t.T),sigma2[l.pred]*K.SE[[l.pred]])
      alpha.pred = 1/(1+exp(-eta.pred))
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
  rr = seq(from = 0.01, to = 2, len = N.MC)
  
  # density of the inverse gamma
  dInvgamma<-function(x,shape, scale){
    ###################################
    # x: value at which we are evaluating
    # shape: shape parameter for density
    # scale: scale parameter for density
    ###################################
    alpha = shape
    beta = scale
    return ( (beta^alpha/gamma(alpha) )*x^(-alpha - 1)*exp(-beta/x) )
  }
  for(l in 1:L){
  df = data.frame(rr, sigma2 = SIGMA2_POST[,l], density = dInvgamma(rr, shape = r_0, scale = s_0[l]) )
  # create a Sigma2 plot and save it
  pdf(paste(Triplet_fig_dir,"Sigma2_",triplet,"_",l,".pdf", sep=""))
  g<-ggplot(df, aes(x=sigma2)) + 
    geom_density(size = 2) + 
    xlab("") + 
    geom_line(aes(rr, density), color = "blue", size = 2)+
    theme(axis.text=element_text(size=20, color="black"),
          axis.title=element_text(size=24,face="bold"), 
          legend.text=element_text(size=20))
  print(g)
  dev.off()
  }
  
  eb_post = 1/(1+exp(-ETA_BAR_POST))
  eb_prior = 1/(1+exp(-ETA_BAR_PRIOR))
  
  df = data.frame(eta_bar = eb_post, prior = eb_prior)
  # create a plot for Eta and save it
  pdf(paste(Triplet_fig_dir,"Eta_Bar_",triplet,".pdf", sep=""))
  g<-ggplot(df, aes(x=eta_bar)) + 
    geom_density(size = 2) + 
    xlab("") + 
    #xlim(0,1)+
    geom_density( aes(prior), color = "blue", size = 2)+
    theme(axis.text=element_text(size=20, color="black"),
          axis.title=element_text(size=24,face="bold"), 
          legend.text=element_text(size=20))
  print(g)
  dev.off()
  
  df = data.frame(replication = 1:nRep, Z.mean = apply(ZETA_POST,2,mean) )
  # create a plot for Zeta and save it
  # this can probably be removed
  pdf(paste(Triplet_fig_dir,"Zeta_Post_",triplet,".pdf",sep="") )
  g<-ggplot(df, aes(x = replication) )+
    geom_point(aes(x=replication,y=Z.mean), size=4)+
    theme(axis.text=element_text(size=20, color="black"),axis.title=element_text(size=24,face="bold"), legend.text=element_text(size=20))
  print(g)
  dev.off()
  
  lw.1 = qgamma(.025, shape = r_A, rate = s_A)
  up.1 = qgamma(.975, shape = r_A, rate = s_A)
  m.1 = r_A/s_A
  
  lw.2 = apply(lambda_A_POST,2,quantile,.025)
  up.2 = apply(lambda_A_POST,2,quantile,.975)
  m.2 = apply(lambda_A_POST,2,mean)
  df = data.frame(Time = 1:t.T, m.1, up.1, lw.1, m.2, up.2, lw.2)
  # create a plot for Lambda_A
  pdf(paste(Triplet_fig_dir,"Lambda_A_Post_",triplet,".pdf",sep="") )
  p1 <- ggplot(df, aes(Time, m.1))+
    geom_point(color = "blue")+
    geom_line(data=df, color = "blue")+
    geom_ribbon(data=df,aes(ymin=lw.1,ymax=up.1),alpha=0.3, color = "blue", fill = "blue")+ 
    geom_line(aes(Time, m.2), color = "red" )+
    geom_point(aes(Time, m.2), color = "red" )+
    geom_ribbon(data=df, aes(ymin=lw.2,ymax=up.2), color = "red", alpha = .3, fill = "red")+
    ylab("Counts")+
    ylim(0,15)+
    theme(axis.text=element_text(size=20, color="black"),
          axis.title=element_text(size=24,face="bold"), 
          legend.text=element_text(size=20))
  print(p1)   
  dev.off()
  
  
  lw.1 = qgamma(.025, shape = r_B, rate = s_B)
  up.1 = qgamma(.975, shape = r_B, rate = s_B)
  m.1 = r_B/s_B
  
  lw.2 = apply(lambda_B_POST,2,quantile,.025)
  up.2 = apply(lambda_B_POST,2,quantile,.975)
  m.2 = apply(lambda_B_POST,2,mean)
  df = data.frame(Time = 1:t.T, m.1, up.1, lw.1, m.2, up.2, lw.2)
  # create a plot for Lambda B
  pdf(paste(Triplet_fig_dir,"Lambda_B_Post_",triplet,".pdf",sep="") )
  p1 <- ggplot(df, aes(Time, m.1))+
    geom_point(color = "blue")+
    geom_line(data=df, color = "blue")+
    geom_ribbon(data=df,aes(ymin=lw.1,ymax=up.1),alpha=0.3, color = "blue", fill = "blue")+ 
    geom_line(aes(Time, m.2), color = "red" )+
    geom_point(aes(Time, m.2), color = "red" )+
    geom_ribbon(data=df, aes(ymin=lw.2,ymax=up.2), color = "red", alpha = .3, fill = "red")+
    ylab("Counts")+
    ylim(0,15)+
    theme(axis.text=element_text(size=20, color="black"),
          axis.title=element_text(size=24,face="bold"),
          legend.text=element_text(size=20))
  print(p1)   
  dev.off()
  
  colnames(Ell_Post) = as.character(ell)
  df = stack(as.data.frame(Ell_Post))
  names(df) = c("Probability", "lengthScale")
  df$lengthScale <- factor(df$lengthScale,levels = ell,ordered = TRUE)
  # create a plot for Ell
  pdf(paste(Triplet_fig_dir,"Post_Predictive_Ell_",triplet,".pdf",sep="") )
  g<-ggplot(df, aes(lengthScale, Probability)) + geom_boxplot()+
    theme(axis.text=element_text(size=20, color="black"),axis.title=element_text(size=24,face="bold"), legend.text=element_text(size=20))
  print(g)
  dev.off()
  
  df = data.frame(MinMax = MinMax.Pred, Prior = MinMax.Prior)
  # create a plot for MinMax
  pdf(paste(Triplet_fig_dir,"MinMax_Pred_",triplet,".pdf", sep=""))
  g<-ggplot(df, aes(x=MinMax)) + 
    geom_density(size = 2) + 
    xlab("") + 
    geom_density( aes(Prior), color = "blue", size = 2)+
    xlim(0,1)+
    theme(axis.text=element_text(size=20, color="black"),
          axis.title=element_text(size=24,face="bold"), 
          legend.text=element_text(size=20))
  print(g)
  dev.off()
  
  for(i in 1:nRep){
    
    df = data.frame(Time = 1:t.T, 
                    mean = alpha_mean[i,], 
                    CI.025 = alpha_pt025[i,], 
                    CI.975 = alpha_pt975[i,],
                    sample1 = ALPHA[[i]][1,],
                    sample2 = ALPHA[[i]][50,], 
                    sample3 = ALPHA[[i]][100,])
    l2 = sqrt(sum( (alpha_mean[i,] - mean(alpha_mean[i,]) )^2 ) )
    l2 = round(l2,2)
    # alpha_triplit plot
    pdf(paste(Triplet_fig_dir,"Alpha_triplet_",triplet,"_",i,".pdf", sep=""))
    g <- ggplot(df, aes(Time, mean))+
      geom_point(color = "black")+
      geom_line(data=df, color = "black")+
      geom_line(data=df,aes(Time,sample1), color = "gray70")+
      geom_line(data=df,aes(Time,sample2), color = "gray70")+
      geom_line(data=df,aes(Time,sample3), color = "gray70")+
      geom_ribbon(data=df,aes(ymin=CI.025,ymax=CI.975),alpha=0.1, color = "black", fill = "black")+ 
      ylim(0,1)+
      xlab("Time")+
      annotate( "text", x=5, y = .9, size = 8, label = paste("L2 = " ,l2, sep=""))+
      theme(axis.text=element_text(size=20, color="black"),
            axis.title=element_text(size=24,face="bold"), 
            legend.text=element_text(size=20))
    print(g)
    dev.off()
    
    df = data.frame(MinMax = MinMax[,i], Prior = MinMax.Prior)
    # make MinMax plot
    pdf(paste(Triplet_fig_dir,"MinMax_",
              triplet,"_",i,".pdf",
              sep=""))
    g<-ggplot(df, aes(x=MinMax)) + 
      geom_density(size=2) + 
      xlab("") + 
      geom_density( aes(Prior), color = "blue", size = 2)+
      xlim(0,1) + 
      theme(axis.text=element_text(size=20, color="black"),
            axis.title=element_text(size=24,face="bold"), 
            legend.text=element_text(size=20))
    print(g)
    dev.off()
    
  }
  # save a data file with the results
  save(file = paste(Triplet_dir,"triplet_",triplet,".RData", sep=""),
       alpha_mean, 
       alpha_pt025,
       alpha_pt975,
       A_mean, 
       A_pt025, 
       A_pt975,
       ALPHA_BAR_POST,
       ETA_BAR_POST, 
       Median.Min.Max)
}

