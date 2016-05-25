MCMC.plot = function(MCMC.results, all.plots = F){
  triplet = MCMC.results$triplet
  ALPHA_BAR_POST = MCMC.results$ALPHA_BAR_POST
  ALPHA = MCMC.results$ALPHA
  SIGMA2_POST = MCMC.results$SIGMA2_POST
  ETA_BAR_POST = MCMC.results$ETA_BAR_POST
  lambda_A_POST = MCMC.results$lambda_A_POST
  lambda_B_POST = MCMC.results$lambda_B_POST
  Ell_Post = MCMC.results$Ell_Post
  A_POST = MCMC.results$A_POST
  MinMax.Pred = MCMC.results$MinMax.Pred
  MinMax.Prior = MCMC.results$MinMax.Prior
  nRep = dim(MinMax)[2]
  
  Triplet_fig_dir = paste(Fig_dir,"Triplet_",triplet,"/",sep="")
  unlink(Triplet_fig_dir, recursive=T)
  dir.create(Triplet_fig_dir)
  
  if (all.plots){
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
      df = data.frame(rr, 
                      sigma2 = SIGMA2_POST[,l], 
                      density = dInvgamma(rr, shape = r_0, scale = s_0[l]) )
      
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
    
    
    ######
    #Lambda plot
    #####
    # calculate posterior quantities
    lw.A = apply(lambda_A_POST,2,quantile,.025)
    up.A = apply(lambda_A_POST,2,quantile,.975)
    m.A = apply(lambda_A_POST,2,mean)
    lw.B = apply(lambda_B_POST,2,quantile,.025)
    up.B = apply(lambda_B_POST,2,quantile,.975)
    m.B = apply(lambda_B_POST,2,mean)
    # put in single data frame
    full_df = data.frame(Time = rep(1:40, 2), 
                         Sound = c(rep("A", 40), rep("B", 40)), 
                         Mean = c(m.A, m.B),
                         low = c(lw.A, lw.B),
                         up = c(up.A, up.B))
    # plot lambda_A and lambda_B, with error bars
    pdf(paste(Triplet_fig_dir,"lambdas_",triplet,".pdf",sep="") )
    p1 = ggplot(full_df, aes(x = Time, y = Mean, group = Sound, color = Sound)) +
      geom_line() +
      geom_ribbon(aes(ymin = low, ymax = up, 
                      group = Sound, fill = Sound), 
                  alpha = .1,
                  colour = NA) +
      ylab("Counts")
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
    
    ######
    # Create MinMax plot for all densities together
    ######
    df = melt(MinMax, varnames = c("Time", "Trail"))
    # make MinMax plot
    pdf(paste(Triplet_fig_dir,"MinMax_all_",triplet,".pdf", sep=""))
    g=ggplot(df) + 
      geom_density(aes(x = value, group = Trail, fill = Trail),
                   alpha = 1/5) +
      xlim(0,1)
    print(g)
    dev.off()
    
    ######
    # Create plot with all alpha_t functions together
    ######
    # process data for ggplot
    trail_names = paste0("Trail", 1:nrow(alpha_mean))
    row.names(alpha_mean) = trail_names
    mean_df = data.frame(alpha_mean, Trail = trail_names)
    colnames(mean_df) = c(1:ncol(alpha_mean), "Trail")
    mean_df = melt(mean_df, id.vars = "Trail", value.name = "mean", variable.name = "Time")
    # make and save plot
    pdf(paste(Triplet_fig_dir,"alpha_joint_",triplet,".pdf", sep=""))  
    g = ggplot(mean_df, aes(x = Time, y = mean, group = Trail)) + geom_line(aes(group = Trail))
    print(g)
    dev.off()
  }
  
}