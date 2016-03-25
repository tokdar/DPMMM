rm(list=ls())
FIRST_USE = T

burnIn = 1e2#500e3
N.MC = 1e2
thin = 1 #300
if(FIRST_USE) install.packages(c("BayesLogit","MCMCpack", "parallel", "ggplot2"), repos = "http://cran.us.r-project.org")

library(BayesLogit)
#library(MCMCpack) ## trying to work a way around MCMCpack
library(bayesm)
library(parallel)
library(ggplot2)
#load synthetic data and hyperparameters

## new function to draw from 
rDirichlet <- function(n, alpha) return(t(replicate(n, rdirichlet(alpha)))[1:n,])

directory = "~/DPMMM-master/"

Code_dir = paste(directory,"Code/",sep="")
Fig_dir = paste(directory,"Figures/",sep="")
Triplet_dir = paste(directory,"Post_Summaries/",sep="")


source(paste(Code_dir,"A_step.R",sep="") )
source(paste(Code_dir,"A_star_step.R",sep="") )
source(paste(Code_dir,"B_star_step.R",sep="") )
source(paste(Code_dir,"lambda_A_step.R", sep="") )
source(paste(Code_dir,"lambda_B_step.R",sep="") )
source(paste(Code_dir,"Omega_step.R",sep="") )
source(paste(Code_dir,"K_Matern.R",sep="") )
source(paste(Code_dir,"eta_Matern.R", sep="") )
source(paste(Code_dir,"ell_prior_step.R",sep="") )
source(paste(Code_dir,"sigma2_Matern_step.R", sep="") )
source(paste(Code_dir,"p_step.R",sep="") )
source(paste(Code_dir,"gamma_step.R",sep="") )
source(paste(Code_dir,"m_gamma_step.R",sep="") )
source(paste(Code_dir,"sigma2_gamma_step.R",sep="") )
source(paste(Code_dir,"pi_gamma_step.R",sep="") )
source(paste(Code_dir,"Bincounts.R",sep="") )
source(paste(Code_dir,"Data_Pre_Proc.R", sep="") )
source(paste(Code_dir,"MCMC_Triplet.R", sep="") )
#parameters for mixture components
K = 5
m_0= 0 
sigma2_0 = 1 #.01
r_gamma = 101
s_gamma = 1
#sampling for sigma2
s_0 = 50
r_0 = s_0+1
#parameters for pi_gamma
alpha_gamma = 1/K

ell = c(1,2,3,4,5,15)
L = length(ell)
ell_0 = c( rep(.5/(L-1),(L-1) ), .5)

Triplet_meta = read.csv("~/Dropbox/Neuro/Data/Trial_ID.csv", stringsAsFactors=F)
Triplet_meta = unique(Triplet_meta)
Triplet_meta = Triplet_meta[order(Triplet_meta[,"SepBF"], decreasing=T),]
triplets = 2

source(paste(Code_dir,"eta_bar_mixture.R",sep="") )
source(paste(Code_dir,"MinMax_Prior.R",sep="") )
#triplets is the index (or row number) of the triplet in the Triplet_Meta dataframe
pt = proc.time()[3]
mclapply(triplets, function(triplet) MCMC.triplet(triplet, ell_0, ETA_BAR_PRIOR, MinMax.Prior), mc.cores = 1)
proc.time()[3] - pt




