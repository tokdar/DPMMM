Pre_Proc<-function(triplet, bin_width, Triplet_meta){


freq = Triplet_meta[triplet,"AltFreq"]
pos = Triplet_meta[triplet,"AltPos"]
cell = Triplet_meta[triplet,"Site"]
fname = Triplet_meta[triplet,"Cell"]
if (cell != 0){
  trials <- read.table(url(paste("http://www2.stat.duke.edu/~st118/Jenni/STCodes/VC/", 
                                 fname,".txt",sep="") ), sep="\t" )
  # seems to be extra column of NAs
  # probably a result of extra spaces
  # we now drop it
  trials = trials[,-dim(trials)[2]]
  # name remaining columns
  colnames(trials) = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL", "SOFF")
  
  url = paste("http://www2.stat.duke.edu/~st118/Jenni/STCodes/VC/", fname,
              "_cell", cell, "_spiketimes.txt",sep="")
  spiketimes <- read.table(url(paste("http://www2.stat.duke.edu/~st118/Jenni/STCodes/VC/",
                                     fname, "_cell", cell, "_spiketimes.txt",sep="") ), sep="\t" )
  
  # data file generates an extra column of NAs
  # we drop this unneeded column here
  spiketimes = spiketimes[,-dim(spiketimes)[2]]
  colnames(spiketimes) = c("TRIAL2", "TIMES")
} else {
  trials <- read.table(url(paste("http://www2.stat.duke.edu/~st118/Jenni/STCodes/JA/",
                                 fname, ".txt", sep = "")), sep = "")
  trials <- cbind(trials, 600)
  # these files do no have an extra column of NAs
  # no need to drop the last column
  #trials = trials[,-dim(trials)[2]]
  colnames(trials) = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL", "SOFF")
  
  url = paste("http://www2.stat.duke.edu/~st118/Jenni/STCodes/JA/", fname,
              "_spiketimes.txt",sep="")
  spiketimes <- read.table(url(paste("http://www2.stat.duke.edu/~st118/Jenni/STCodes/JA/",
                                     fname, "_spiketimes.txt",sep="") ), sep="" )
  
  # there is no uneeded column of NAs
  # no need to drop the last column
  # spiketimes = spiketimes[,-dim(spiketimes)[2]]
  colnames(spiketimes) = c("TRIAL2", "TIMES")
}
AB = Bincounts(trials, spiketimes, frq = c(freq, 742), pos = c(pos, -144/pos), on.reward = TRUE, start.time = 0, end.time = 1000, bw = bin_width, target = c(25, 150, 800), match.level = FALSE, AB.eqlevel = FALSE, go.by.soff = TRUE, n.iter = 1e3, plot = FALSE, faux.dual.mix = FALSE, faux.dual.int = FALSE, faux.alpha = 0.5, faux.dual.swi = FALSE, faux.target = c(100, 100), nAB = "match.realABcount")


X = AB[[3]]


x1 = AB[[1]]
x2 = AB[[2]]
x3 = AB[[3]]

t.T <- nrow(X)
n1 <- ncol(x1)
n2 <- ncol(x2)
n3 <- ncol(x3)

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


## new style -- adaptive smoothing with gam
get.hyper1 <- smoogam(x1);
get.hyper2 <- smoogam(x2);
m.1 <- get.hyper1$mean; am.1 <- n1 * m.1; bm.1 <- rep(n1, t.T); s.1 <- sqrt(am.1)/bm.1
m.2 <- get.hyper2$mean; am.2 <- n2 * m.2; bm.2 <- rep(n2, t.T); s.2 <- sqrt(am.2)/bm.2

r_A = am.1; s_A = bm.1
r_B = am.2; s_B = bm.2

ret_X = list(X,r_A,s_A,r_B,s_B)
return(ret_X)
}