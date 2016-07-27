# clear the environment
rm(list=ls())
# if set to T, will install needed packages
# once packages are installed, no need to set to T
FIRST_USE = F
if(FIRST_USE) install.packages(c("BayesLogit", "parallel", "ggplot2"), 
                               repos = "http://cran.us.r-project.org")
exper.name <- "MixA400B100Sine50CompA75B25TestSigma2"
expir.desc <- "Run experiments with fixed shift in smooth multiplexing"
#############################
# generate synthetic data
#############################
# Select control types from this list
# control.type = c("Alike", "weak_switch80","weak_switch90", 
#                  "switch", "average","wt_average","avg_mplxd",
#                   "wt_avg_mplxd", "outside", "smoo_mplxd")
control.type = c("average", "smoo_mplxd")
mix_bool = T
#control.types = rep(control.type, 10)
control.types = sample(control.type, size = 20, prob = c(.5, .5), replace = TRUE)
#number of trials to be simulated for each triplet
n_trials_vec=c(20) 
#number of files(triplets) to be generated for each condition
repeats = 8   
#Note: the spike generator is made to accept lambda in a vector form 
# with lambda for each time bin,
# baseline for A = 400, B = 100
lambdaASeed=400  
lambdaBSeed=100 

#sin.period <- rep(400, n_trials_vec)
sin.period <- seq(from = 200, to = 1200, length.out = 20)
jit <- sample(1:sin.period, n_trials_vec)
compressA <- rep(.75, n_trials_vec)
compressB <- rep(.25, n_trials_vec)
out.dir = paste0("Data/", exper.name, "/")
list.file.names = paste0(out.dir, "Synth_names.txt")

source("Code/Make_Synth_Data.R")
if (mix_bool == FALSE){
  control.types = control.type
  make_synth(control.type, n_trials_vec, repeats, lambdaASeed,
             lambdaBSeed, out.dir, list.file.names)
} else {
  
  make_synth_mix(control.types, n_trials_vec, repeats, lambdaASeed,
                 lambdaBSeed, out.dir, list.file.names,
                 bw = 200,
                 start.time = -200,
                 end.time = 1000,
                 jit,
                 sin.period,
                 compressA,
                 compressB)
}

# set burn in
burnIn = 3000
# set number of iterations
N.MC = 2000
# set thinning rate
# if set to 1, no thinning
thin = 1
# number of cores to use if running in parallel
nCores = 8



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
  sm <- rowSums(x)
  return(x/sm)
}

rinvgamma <-function(n,shape, scale = 1) return(1/rgamma(n = n, shape = shape, rate = scale))

rtruncgamma <- function(n, low = 0, up = Inf, ...){
  unifs <- runif(n)
  return(qgamma(pgamma(low, ...) + unifs*(pgamma(up, ...) - pgamma(low, ...)), ...))
}

smoogam <- function(x){
  require(mgcv)
  T <- nrow(x)
  n <- ncol(x)
  if(T > 1){
    x.dd <- data.frame(cts = c(x), time = rep(1:T, n))
    x.gam <- gam(cts ~ s(time, bs = "ad"), data = x.dd, family = poisson(link = "log"))
    x.pred <- predict(x.gam, data.frame(cts = NA, time = 1:T), se.fit = TRUE)
    mu <- x.pred$fit; sig <- x.pred$se.fit
  } else {
    mu <- log(mean(c(x))); sig <- log(sd(c(x)))
  }
  firingrate.mean <- exp(mu)
  firingrate.vari <- expm1(sig^2/2) * firingrate.mean^2
  return(list(a = 1 / expm1(sig^2/2), b = 1/(expm1(sig^2/2) * firingrate.mean), 
              mean = firingrate.mean, vari = firingrate.vari))
}

#directory = "~/DPMMM/"

# this should be modified for increased flexibility
# right now it requires a specific directory path
# so the code is not portable
# it is also inconvinient because this file is not in the root directory

# directory on my local machine
# directory = "/Users/azeemzaman/Documents/Research/Neuro/DPMMM/"
# directory on Saxon
directory = "/home/grad/azz2/Research/DPMMM/"
Code_dir = paste(directory,"Code/",sep="")
Fig_dir = paste(directory,"Figures/", exper.name, "/", sep="")
Triplet_dir = paste(directory,"Post_Summaries/", exper.name, "/", sep="")
for (direc in c(Code_dir, Fig_dir, Triplet_dir)){
  if (!file.exists(direc)){
    dir.create(direc)
  }
}

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
source(paste(Code_dir,"MCMC_Triplet_raw.R", sep="") )
source(paste(Code_dir,"MCMC_plot.R", sep="") )
source(paste(Code_dir,"Count_Switches.R", sep="") )
source(paste(Code_dir,"Data_Merge.R", sep="") )

#parameters for mixture components
K = 5
m_0= 0 
sigma2_0 = 1 #.01
r_gamma = 101
s_gamma = 1

#ell = c(1,2,3,4,5,15)
# using new length scales
ell = c(3, 6, 8, 20, 80)
L = length(ell)
#ell_0 = c( rep(.5/(L-1),(L-1) ), .5)
ell_0 = rep(1/L, L)
#sampling for sigma2
delta = 2e4
#r_0 = 51
#s_0 = (r_0 - 1)*(1-exp(-delta^2/ell^2)) 
r_0 = 0.1
s_0 = rep(r_0, L)

#parameters for pi_gamma
alpha_gamma = 1/K

# radius values for switch counts
widthes = seq(from = 0.01, to = .15, length.out = 20)

# read data from Surja's website
# Triplet_meta = read.csv("http://www2.stat.duke.edu/~st118/Jenni/STCodes/ResultsV2/All-HMM-Poi-selected.csv", 
#                        stringsAsFactors=F)
# this is the file location on my laptop
# Triplet_meta = read.csv("/Users/azeemzaman/Documents/Research/Neuro/DPMMM/Triplets_pass_criteria.csv")
# this is the file address on Saxon
#Triplet_meta = read.csv("/home/grad/azz2/Research/DPMMM/Triplets_pass_criteria.csv")
# file taking into account already run Triplets
#Triplet_meta = read.csv("/home/grad/azz2/Research/DPMMM/Filtered_All_HMM.csv")
#Triplet_meta = unique(Triplet_meta)
#Triplet_meta = Triplet_meta[order(Triplet_meta[,"SepBF"], decreasing=T),]
# Triplet_meta = Triplet_meta[order(Triplet_meta[,"WinPr"], decreasing=T),]
#triplets = 321:417

# read in the names of tiles
root.dir <- paste0("Data/", exper.name)
list.file.loc <- paste0(directory, root.dir, list.file.names)
file.names <- as.character(read.table(list.file.names)$V1)
full.names <- paste0(root.dir, "/", file.names)
triplet.data <- list()
for (i in 1:length(full.names)){
  triplet.data[[i]] <- list()
  triplet.data[[i]][[1]] <- i
  triplet.data[[i]][[2]] <- read.table(paste0(full.names[i], ".txt"))
  colnames(triplet.data[[i]][[2]]) <- c("TRIAL", "TASKID", "A_FREQ", "B_FREQ",
                                        "XA", "XB", "REWARD", "A_LEVEL",
                                        "B_LEVEL", "SOFF")

  triplet.data[[i]][[3]] <- read.table(paste0(full.names[i], "_spiketimes", ".txt"))
  colnames(triplet.data[[i]][[3]]) <- c("TRIAL2", "TIMES")
}

source(paste(Code_dir,"eta_bar_mixture.R",sep="") )
source(paste(Code_dir,"MinMax_Prior.R",sep="") )

# small section of test code
# compile to byte code
burnIn = 10
N.MC = 50
test_run = MCMC.triplet.raw(triplet.data[[1]][[1]], triplet.data[[1]][[2]], 
                            triplet.data[[1]][[3]], ell_0, ETA_BAR_PRIOR, 
                            MinMax.Prior)
MCMC.plot(test_run, F, 10, widthes)

# reset burnin and N.MC
burnIn = 3000
N.MC = 2000
#triplets is the index (or row number) of the triplet in the Triplet_Meta dataframe
pt = proc.time()[3]
MCMC.results = mclapply(triplet.data, function(triplet) {try(MCMC.triplet.raw(triplet[[1]], triplet[[2]], triplet[[3]], ell_0, ETA_BAR_PRIOR, MinMax.Prior))}, mc.cores = nCores)
proc.time()[3] - pt
mclapply(MCMC.results, function(x) {try(MCMC.plot(x, F, 100, widthes, T))}, mc.cores = nCores)

# collect error messages
errors = which(sapply(MCMC.results, typeof) == "character")
for (error in errors){
  #cat("\n", file = "raw_log.txt", append = TRUE)
  #cat(colnames(Triplet_meta)[1:11], file = "raw_log.txt", append = TRUE)
  #cat("\n", file = "raw_log.txt", append = TRUE)
  #cat(paste(sapply(Triplet_meta[error,1:11], as.character), sep = ","), file = "raw_log.txt", append = TRUE)
  cat("\n", file = "raw_log.txt", append = TRUE)
  cat(as.character(MCMC.results[[error]]), file = "raw_log.txt", append = TRUE)
  cat("===============", file = "raw_log.txt", append = TRUE)
  #print(MCMC.results[[error]])
}

# clean up and conversion
image.dir <- paste0("Images/", exper.name)
if (!file.exists(image.dir)){
  dir.create(image.dir)
}
root.dir <- paste0("Data/", exper.name)
list.file.loc <- paste0(directory, root.dir, list.file.names)
file.names <- as.character(read.table(list.file.names)$V1)
file.names <- sapply(file.names, function(x) {strsplit(x, "/")[[1]][1]})
names_ext <- paste0(file.names, ".png")
curr_files <- paste0("Summary_Triplet",1:length(full.names), ".pdf-1.png")
cbind(curr_files, names_ext)
conv_csv <- paste0(image.dir, "/rename.csv")
write.csv(cbind(curr_files, names_ext), file = conv_csv, row.names = F, col.names = F)
clean_file <- paste0(image.dir, "/cleanup.sh")
str1 = paste0("START=1
END=",length(file.names))
str2 = paste0('EXT=".pdf"
for ((i=START;i<=END;i++));
do
file="../../Figures/', exper.name,'/Triplet_$i/Summary_Triplet$i"')
str3 = 'file+=".pdf"	
echo "Copying $file"
cp $file ./
  echo "Converting Summary_Triplet$i$EXT "
file="Summary_Triplet$i$EXT"
pdftoppm -rx 300 -ry 300 -png "$file" "$file"
done'
str4 = "awk -F, \'{print("
str5 = '"mv \\"" $1 "\\" \\"" $2 "\\"")}'
str6 = "' rename.csv | bash -"
full_file = paste0(paste(str1, str2, str3, str4, sep = "\n"), str5, str6)
cat(full_file, file = clean_file)
