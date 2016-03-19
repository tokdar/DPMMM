## Display raw spike trains from VC data
rm(list=ls())
fname = "YKIC140131Loc_DoubleSound"

display.raw.VC <- function(fname, cell){
  
  infile1 <- paste("~/Dropbox/Neuro/Data/", fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL", "SOFF"))
  
  infile2 <- paste("~/Dropbox/Neuro/Data/", fname, "_cell", cell, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))
  
  attach(trials)
  attach(spiketimes)
  timestamps <- split(TIMES, TRIAL2)
  ntrials <- length(timestamps)  
  trial.id <- as.numeric(names(timestamps)) ## same as unique(TRIAL2)
  
  layout(matrix(c(rep(1,12),rep(2,4)), 4, 4, byrow = FALSE))
  
  plot(0,0,ty = "n", xlim = range(spiketimes), ylim = range(TRIAL), bty = "n", xlab = "time", ylab = "trial")  
  reward.locator <- rep(NA, ntrials)
  for(i in 1:ntrials){
    jj <- trial.id[i]
    jj2 <- which(TRIAL == jj)
    spks <- timestamps[[i]]
    nspk <- length(spks)
    reward.locator[i] <- REWARD[jj2]
    if(REWARD[jj2]) points(spks, rep(jj, nspk), pch = ".", col = c("gray", "red", "orange")[1 + (spks > 0) + (spks > SOFF[jj2])])  
  }
  
  spk.counts <- sapply(timestamps, length)
  plot(0,0,ty = "n", xlim = c(0,max(spk.counts)), ylim = range(TRIAL), bty = "n", xlab = "spike count", ylab = "", axes = FALSE)
  axis(1)
  segments(rep(0, ntrials), trial.id, spk.counts, trial.id, col = c("gray", "magenta")[1 + reward.locator])
  
  title(main = paste(fname, " (cell ", cell, ")", sep = ""), out = TRUE, line = -2)
  detach(trials)
  detach(spiketimes)
}

pdf("~/Dropbox/Neuro/Figures/Trials_Raw.pdf")
display.raw.VC("YKIC140131Loc_DoubleSound", 1)
dev.off()

## display raw spikes from a triplet
display.triplet.VC <- function(fname, cell, Afrq, Apos, on.reward = FALSE){
  
  infile1 <- paste("~/Dropbox/Neuro/Data/", fname, ".txt", sep = "")
  trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL", "SOFF"))
  
  infile2 <- paste("~/Dropbox/Neuro/Data/", fname, "_cell", cell, "_spiketimes.txt", sep = "")
  spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))
  
  attach(trials)
  attach(spiketimes)
  timestamps.all <- split(TIMES, TRIAL2)
  
  frq <- c(Afrq, 742); pos <- c(Apos, -144/Apos)
  ix1 <- TASKID == 8 & A_FREQ == frq[1] & XA == pos[1]
  ix2 <- TASKID == 8 & A_FREQ == frq[2] & XA == pos[2]
  ix3 <- TASKID == 12 & (A_FREQ == frq[1] & B_FREQ == frq[2] & XA == pos[1] & XB == pos[2]) | (A_FREQ == frq[2] & B_FREQ == frq[1] & XA == pos[2] & XB == pos[1])
  
  if(on.reward){
    ix1 <- ix1 & REWARD == 1
    ix2 <- ix2 & REWARD == 1
    ix3 <- ix3 & REWARD == 1
  } 
  groups <- list(trials[ix1, 1], trials[ix2, 1], trials[ix3, 1])
  
  layout(matrix(rep(1:6, c(2,1,2,1,2,1)), 3,3, byrow = TRUE))
  
  group.names <- c(paste("A @", frq[1], "Hz", pos[1], "deg", sep = ""),
                   paste("B @", frq[2], "Hz", pos[2], "deg", sep = ""),
                   "AB")
  for(gg in 1:3){
    gg.sel <- na.omit(match(groups[[gg]], names(timestamps.all)))
    timestamps <- timestamps.all[gg.sel]
    ntrials <- length(timestamps)  
    trial.id <- as.numeric(names(timestamps.all))[gg.sel] ## same as unique(TRIAL2)
    
    plot(0,0,ty = "n", xlim = range(spiketimes), ylim = range(trial.id), bty = "n", xlab = "time", ylab = "trial")  
    reward.locator <- rep(NA, ntrials)
    for(i in 1:ntrials){
      jj <- trial.id[i]
      jj2 <- which(TRIAL == jj)
      spks <- timestamps[[i]]
      nspk <- length(spks)
      reward.locator[i] <- REWARD[jj2]
      if(REWARD[jj2]) points(spks, rep(jj, nspk), pch = ".", col = c("gray", "red", "orange")[1 + (spks > 0) + (spks > SOFF[jj2])])  
    }
    title(main = group.names[gg], font = 1, line = 0)    
    
    spk.counts <- sapply(timestamps, length)
    plot(0,0,ty = "n", xlim = c(0,max(spk.counts)), ylim = range(TRIAL), bty = "n", xlab = "spike count", ylab = "", axes = FALSE)
    axis(1)
    segments(rep(0, ntrials), trial.id, spk.counts, trial.id, col = c("gray", "royalblue")[1 + reward.locator])
  }
  title(main = paste(fname, " (cell ", cell, ")", sep = ""), out = TRUE, line = -2)
  detach(trials)
  detach(spiketimes)
}

pdf("~/Dropbox/Neuro/Figures/Triplet_Raw.pdf")
display.triplet.VC("YKIC140131Loc_DoubleSound", 1, 1100, 6, on.reward = FALSE)
dev.off()

cell = 1
infile1 <- paste("~/Dropbox/Neuro/Data/", fname, ".txt", sep = "")
trials <- read.table(infile1, col.names = c("TRIAL", "TASKID", "A_FREQ", "B_FREQ", "XA", "XB", "REWARD", "A_LEVEL", "B_LEVEL", "SOFF"))

infile2 <- paste("~/Dropbox/Neuro/Data/", fname, "_cell", cell, "_spiketimes.txt", sep = "")
spiketimes <- read.table(infile2, col.names = c("TRIAL2", "TIMES"))

FREQS <- unique(c(trials$A_FREQ, trials$B_FREQ))
alt.freq <- sort(FREQS[FREQS > 0 & FREQS != 742])
alt.pos <- c(-24, -6, 6, 24)


frq = c(1100, 742)
pos = c(24, -6)
on.reward = TRUE
start.time = 0
end.time = 600
bw = 25
target = c(25, 150, 800)
match.level = FALSE
AB.eqlevel = FALSE
go.by.soff = TRUE
n.iter = 1e3
plot = FALSE
faux.dual.mix = FALSE
faux.dual.int = FALSE
faux.alpha = 0.5
faux.dual.swi = FALSE
faux.target = c(100, 100)
nAB = "match.realABcount"
bin.counter <- function(x, b) return(diff(sapply(b, function(a) sum(x <= a))))


Bincounts <- function(trials, spiketimes, frq = c(1100, 742), pos = c(24, -6), on.reward = TRUE, start.time = 0, end.time = 600, bw = 25, target = c(25, 150, 800), match.level = FALSE, AB.eqlevel = FALSE, go.by.soff = FALSE, n.iter = 1e3, plot = FALSE, faux.dual.mix = FALSE, faux.dual.int = FALSE, faux.alpha = 0.5, faux.dual.swi = FALSE, faux.target = c(100, 100), nAB = "match.realABcount", ...){
  
  attach(trials)
  attach(spiketimes)
  
  timestamps <- split(TIMES, TRIAL2)
  ntrials <- length(timestamps)
  trial.id <- as.numeric(names(timestamps)) ## same as unique(TRIAL2)
  
  ix1 <- TASKID == 8 & A_FREQ == frq[1] & XA == pos[1]
  ix2 <- TASKID == 8 & A_FREQ == frq[2] & XA == pos[2]
  ix3 <- TASKID == 12 & (A_FREQ == frq[1] & B_FREQ == frq[2] & XA == pos[1] & XB == pos[2]) | (A_FREQ == frq[2] & B_FREQ == frq[1] & XA == pos[2] & XB == pos[1])
  
  if(on.reward){
    ix1 <- ix1 & REWARD == 1
    ix2 <- ix2 & REWARD == 1
    ix3 <- ix3 & REWARD == 1
  } 
  

  
  if(min(sum(ix1), sum(ix2), sum(ix3)) == 0){
    detach(trials)
    detach(spiketimes)
    stop("Not enought data")
  }
  
  sing1 <- trials[ix1, 1]
  sing2 <- trials[ix2, 1]
  success <- REWARD[ix3]
  
  if(go.by.soff) end.time <- min(SOFF[ix1 | ix2 | ix3])
  
  if(is.nan(end.time)){
    detach(trials)
    detach(spiketimes)    
    stop("SOFF is NaN")
  }
  
  
  
  brk <- seq(start.time, end.time, bw)
  mpt <- (brk[-1] + brk[-length(brk)]) / 2
  
  spike.bincounter <- function(jj, brk){
    jj1 <- match(jj, trial.id)
    spks <- timestamps[[jj1]]
    return(bin.counter(spks, brk))
  }
  
  Abincounts <- matrix(sapply(sing1, spike.bincounter, brk = brk), nrow = length(mpt))
  Bbincounts <- matrix(sapply(sing2, spike.bincounter, brk = brk), nrow = length(mpt))
  nA <- ncol(Abincounts)
  nB <- ncol(Bbincounts)
  min.samp.size <- min(nA, nB)
  
  nAB.target <- ifelse(nAB == "match.realABcount", sum(ix3), nAB)
  duplx <- trials[ix3, 1]
  dualbincounts <- matrix(sapply(duplx, spike.bincounter, brk = brk), nrow = length(mpt))
  if(faux.dual.mix){
    A.subsamp <- sample(nA, min(ceiling(nAB.target*faux.alpha), nA))
    B.subsamp <- sample(nB, min(ceiling(nAB.target*(1-faux.alpha)), nB))
    ABbincounts <- cbind(Abincounts[,A.subsamp], Bbincounts[,B.subsamp])
  } else if (faux.dual.int) {
    ABbincounts <- replicate(nAB.target, thincombo(Abincounts[,sample(nA,1)], Bbincounts[,sample(nB,1)], faux.alpha))
  } else if (faux.dual.swi) {
    #faux.theta <- bw / pmax(bw, faux.target)
    ABbincounts <- replicate(nAB.target, swicombo(Abincounts[,sample(nA,1)], Bbincounts[,sample(nB,1)], bw = bw, p1 = faux.alpha, mean.stay = faux.target, sd.stay = 0.25 * faux.target))    
  } else {
    ABbincounts <- dualbincounts
  }
  
  detach(trials)
  detach(spiketimes)
  
  nAB <- ncol(ABbincounts)
  return(Abincounts, Bbincounts, ABbincounts)
}




