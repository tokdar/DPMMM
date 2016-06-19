# clear the environment
rm(list=ls())
# if set to T, will install needed packages
# once packages are installed, no need to set to T
FIRST_USE = F
if(FIRST_USE) install.packages(c("BayesLogit", "parallel", "ggplot2"), 
                               repos = "http://cran.us.r-project.org")

# set burn in
burnIn = 25e3
# set number of iterations
N.MC = 25e3
# set thinning rate
# if set to 1, no thinning
thin = 1
# number of cores to use if running in parallel
nCores = 4



# load needed libraries
library(BayesLogit)
library(parallel)
library(ggplot2)
library(MASS)
library(reshape2)
library(dplyr)
library(compiler)
enableJIT(3)
## rdirichlet and rinvgamma scraped from MCMCpack
rDirichlet <- function(n, alpha){
#######################################
# This function samples from a dirichlet distribution
#
# Args:
#   n: number of samples
#   alpha: parameter vector
#
# Return:
#   A matrix with n rows where each row is a probability vector
#   whose length is equal to the length of alpha
#######################################
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}

rinvgamma <-function(n,shape, scale = 1) return(1/rgamma(n = n, shape = shape, rate = scale))

#directory = "~/DPMMM/"

# this should be modified for increased flexibility
# right now it requires a specific directory path
# so the code is not portable
# it is also inconvinient because this file is not in the root directory

directory = "/Users/azeemzaman/Documents/Research/Neuro/DPMMM/"
Code_dir = paste(directory,"Code/",sep="")
Fig_dir = paste(directory,"Figures_Pass/",sep="")
Triplet_dir = paste(directory,"Post_Summaries_Pass/",sep="")


source(paste(Code_dir,"A_step.R",sep="") )
source(paste(Code_dir,"A_star_step.R",sep="") )
source(paste(Code_dir,"B_star_step.R",sep="") )
source(paste(Code_dir,"lambda_A_step.R", sep="") )
source(paste(Code_dir,"lambda_B_step.R",sep="") )
source(paste(Code_dir,"Omega_step.R",sep="") )
source(paste(Code_dir,"K_Matern.R",sep="") )
source(paste(Code_dir,"eta_Matern.R", sep="") )
source(paste(Code_dir,"eta_Matern_mod.R", sep="") )
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
source(paste(Code_dir,"MCMC_plot.R", sep="") )
source(paste(Code_dir,"Count_Switches.R", sep="") )
source(paste(Code_dir,"Data_Merge.R", sep="") )

#parameters for mixture components
K = 5
m_0= 0 
sigma2_0 = 1 #.01
r_gamma = 101
s_gamma = 1

ell = c(1,2,3,4,5,15)
L = length(ell)
ell_0 = c( rep(.5/(L-1),(L-1) ), .5)

#sampling for sigma2
delta = 2
r_0 = 51
s_0 = (r_0 - 1)*(1-exp(-delta^2/ell^2)) 

#parameters for pi_gamma
alpha_gamma = 1/K

# radius values for switch counts
widthes = seq(from = 0.01, to = .15, length.out = 20)

# read data from Surja's website
# Triplet_meta = read.csv("http://www2.stat.duke.edu/~st118/Jenni/STCodes/ResultsV2/All-HMM-Poi-selected.csv", 
#                        stringsAsFactors=F)
Triplet_meta = read.csv("/Users/azeemzaman/Documents/Research/Neuro/DPMMM/Triplets_pass_criteria.csv")
Triplet_meta = unique(Triplet_meta)
# Triplet_meta = Triplet_meta[order(Triplet_meta[,"SepBF"], decreasing=T),]
Triplet_meta = Triplet_meta[order(Triplet_meta[,"WinPr"], decreasing=T),]
triplets = 1:4

source(paste(Code_dir,"eta_bar_mixture.R",sep="") )
source(paste(Code_dir,"MinMax_Prior.R",sep="") )

#triplets is the index (or row number) of the triplet in the Triplet_Meta dataframe
pt = proc.time()[3]
MCMC.results = mclapply(triplets, function(triplet) {try(MCMC.triplet(triplet, ell_0, ETA_BAR_PRIOR, MinMax.Prior))}, mc.cores = nCores)
proc.time()[3] - pt
mclapply(MCMC.results, function(x) {try(MCMC.plot(x, F, 2, widthes))}, mc.cores = nCores)

# collect error messages
errors = which(sapply(MCMC.results, typeof) == "character")
log.file = file("log.txt")
sink(log.file, append = TRUE)
for (error in errors){
  print(Triplet_meta[error,1:11])
  print(MCMC.results[[error]])
}
sink()


