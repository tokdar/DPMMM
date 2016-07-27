single_draw_plot <- function(n_samps, ALPHA, trails_to_sample, add_true_sin = F){
  bins = ncol(ALPHA[[1]])
  nRep = length(ALPHA)
  sampled_series = data.frame(Trail = 0, Sample = 0, matrix(NA, nrow = 1, ncol = bins))
  N = nrow(ALPHA[[1]])
  if (nRep < trails_to_sample){
    sampled_trails = 1:nRep
  } else {
    sampled_trails = sample(1:nRep, trails_to_sample)
  }
  for (i in sampled_trails){
    rows = sample(1:N, n_samps)
    trail_name = paste0("Trail_",i)
    sample_alpha_is = ALPHA[[i]][rows,]
    named_samples = data.frame(Trail = trail_name, 
                               Sample = paste0(trail_name, 1:n_samps, trail_name),
                               sample_alpha_is)
    sampled_series = rbind(sampled_series, named_samples)
  }
  sampled_series = sampled_series[-1,]
  colnames(sampled_series) <- c("Trail","Sample", 1:bins)
  plot_df = melt(sampled_series, variable.name = "Time")
  plot_df$Time = 25*as.numeric(plot_df$Time) - 12.5
  if (add_true_sin){
    for (i in sampled_trails){
      this_jit = jit[i]
      this_period = sin.period[i]
      time_span = range(plot_df$Time)
      scaled_fac = compressA[i] - compressB[i]
      shift_fac = compressB[i]
      Time = time_span[1]:time_span[2]
      this_type = control.types[i]
      if (this_type == "smoo_mplxd"){
        value = (sin((Time + this_jit)/this_period*2*pi) + 1)/2*scaled_fac + shift_fac
      } else {
        value = rep(1/2, length(Time))
      }
      sample_int = as.integer(-i)
      sin_df = data.frame(Time = Time, value = value, 
                          Sample = sample_int, Trail = paste0("Trail_",i))
      plot_df = rbind(sin_df, plot_df)
    }

  }
  g = ggplot(plot_df, aes(x = Time, 
                          y = value,
                          group = Sample,
                          colour = Trail)) +
      geom_line(data = plot_df[plot_df$Sample > 0,], alpha = .1) +
      geom_line(data = plot_df[plot_df$Sample <= 0,], size = 1.5, alpha = .5) +
      xlab("Time (ms)") +
      ylab("Alpha") +
      theme(legend.position = "none")

  return(g)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

MCMC.plot = function(MCMC.results, all.plots = F, n_samps, widthes, add_true_sin = F){
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
  ALPHA_pred = MCMC.results$ALPHA_pred
  nRep = length(ALPHA)
  N = nrow(ALPHA[[1]])
  bins = ncol(ALPHA[[1]])
  
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
    full_df = data.frame(Time = rep(1:bins, 2), 
                         Sound = c(rep("A", bins), rep("B", bins)), 
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
  
# make single page of plots
  
  # plot of alpha draws
  p1 = single_draw_plot(n_samps, ALPHA, 2, add_true_sin)
  
  # plot of distributions for length scales
  colnames(Ell_Post) = as.character(ell)
  df = stack(as.data.frame(Ell_Post))
  names(df) = c("Probability", "lengthScale")
  df$lengthScale <- factor(df$lengthScale,levels = ell,ordered = TRUE)
  p2 = ggplot(df, aes(lengthScale, Probability)) + 
    geom_violin(scale = "width",aes(fill = lengthScale), alpha = 1/5) +
    theme(legend.position = "none")
  
  # plot of alpha_bar
  alpha_bar_df = data.frame(alpha_bar =1/(1+exp(-ETA_BAR_POST)))
  p3 = ggplot(alpha_bar_df, aes(alpha_bar)) + geom_density() +
      theme(legend.position = "none") 
  
  # plot of posterior predictive for future alpha*
  df = melt(data.frame(MinMax = MinMax.Pred, Prior = MinMax.Prior))
  p4 = ggplot(df, aes(value, 
                 group = variable,
                 colour = variable,
                 fill = variable)) + geom_density(alpha = 1/5)
  
  # plot lambdas
  lw.A = bins*apply(lambda_A_POST,2,quantile,.025)
  up.A = bins*apply(lambda_A_POST,2,quantile,.975)
  m.A = bins*apply(lambda_A_POST,2,mean)
  lw.B = bins*apply(lambda_B_POST,2,quantile,.025)
  up.B = bins*apply(lambda_B_POST,2,quantile,.975)
  m.B = bins*apply(lambda_B_POST,2,mean)
  # put in single data frame
  full_df = data.frame(Time = rep(1:bins, 2), 
                       Sound = c(rep("A", bins), rep("B", bins)), 
                       Mean = c(m.A, m.B),
                       low = c(lw.A, lw.B),
                       up = c(up.A, up.B))
  full_df$Time = 25*full_df$Time
  # plot lambda_A and lambda_B, with error bars
  p5 = ggplot(full_df, aes(x = Time, y = Mean, group = Sound, color = Sound)) +
    geom_line() +
    geom_ribbon(aes(ymin = low, ymax = up, 
                    group = Sound, fill = Sound), 
                alpha = .1,
                colour = NA) +
    ylab("Firing Rate (Hz)") +
    xlab("Time (ms)")
  
  
  switch_counts = switch_matrix(ALPHA_pred, widthes)
  switch_mean = apply(switch_counts, 1, mean)
  switch_975 = apply(switch_counts, 1, function(x) {quantile(x, .975)})
  switch_25 = apply(switch_counts, 1, function(x) {quantile(x, .025)})
  switch_df = data.frame(Width = 2*widthes,
                         mean = switch_mean,
                         low = switch_25,
                         up = switch_975)
  switch_full = melt(data.frame(width = 2*widthes, switch_counts),
                     id.vars = "width")
  
  p6 = ggplot(switch_df, aes(x = Width, y = mean)) +
    geom_line(size = 2) +
    geom_ribbon(aes(ymin = low, ymax = up), 
                alpha = .1,
                fill = "red") +
    geom_jitter(data = switch_full, aes(x = width, y = value, group = width),
                alpha = 1/500, colour = "blue", size = 1/2) +
    ylab("Number of Switches") +
    xlab("Switch Distance")
  pdf(paste(Triplet_fig_dir,"Summary_Triplet",triplet,".pdf", sep=""))
  multiplot(p1, p2, p3, p4, p5, p6, cols = 2)
  dev.off()
}