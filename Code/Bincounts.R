

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
  
  
  
  #if(min(sum(ix1), sum(ix2), sum(ix3)) == 0){
  #  detach(trials)
  #  detach(spiketimes)
  #  stop("Not enought data")
  #}
  
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
  RETX = list(Abincounts, Bbincounts, ABbincounts)
  return(RETX)
}


