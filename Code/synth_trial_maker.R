#this script fabricates data for testing of the analysis and puts it
# in the form used by the regular analysis


## Codes needed
source("poi_spike_generation.R")


#make_trial_data
#take in the trial indices and fill out a matrix with the correct values for that type of trial
#Inputs: trial indices, trial.type(A, B, or AB)
#optional: frq or indices, to track features of control (frq = two lambdas, pos can be used to make many triplets, but usually just pick 6 and -24)
make_trial_data <- function(trial.nums, trial.type = "AB", pos=c(0,-2859)){
  for (i in trial.nums){
    if (trial.type=="AB") {
      taskid=12
      this.frq=c(1,2)
      Alevel=50
      Blevel=50
    }
    else if (trial.type=="A") {
      taskid=8
      this.frq=c(1,0)
      Alevel=50
      Blevel=0
    }
    else {
      taskid=8
      this.frq=c(2,0)
      Alevel=0
      Blevel=50
    }
    
    this.trial <- c(i,taskid,this.frq,pos,1,Alevel,Blevel,1000)
    
    if (exists('trial.data')) trial.data<-rbind(trial.data,this.trial)
    else trial.data<-this.trial
  }
  trial.data
  
}

#make_synth_data
#inputs: vector of trial indices(trial.nums), rates(lambdaA and lambdaB, vectors in hz), trial type(control.type) 
#outputs: spike times in the form of a nx2 matrix
#added.spikes function variables: lambdas=vectors of firing rate in hz, weight factor = % of A like trials, mplx.rate = % A like bins within trial (if 0 then no multiplexing),
# ...       num.ctrials = number of spike trains to generate (set to 1 here to fascilitate labeling)
#control.types: Alike -> lambda = A lambda, B like ->lambda = B lambda, mplx -> across trial lambda switches, average -> lambda = A+B/2, avg_mplx ->within trial switching, weak_mplx -> across trial weakly switching


make_synth_data <- function(trial.nums, start.time = -200, end.time = 800, bw = 25,lambdaA=rep(0,(end.time-start.time)/bw),lambdaB=rep(0,(end.time-start.time)/bw),control.type=NULL, ...){
  #spikemaker
  
  if(!is.null(control.type)){
    if (control.type == "smoo_mplxd"){
      ABspikes <- smooth_mplxd_wrap(trial.nums, 
                                    lambdaA = lambdaA, 
                                    lambdaB=lambdaB, 
                                    time.bins = seq(from = start.time/1000,
                                                    to = end.time/1000,
                                                    by = 0.001))
    } else {
    for (ABindex in trial.nums) {
      if (control.type == "Alike") added.spikes<-gen_synth_binwise(lambdaA=lambdaA, lambdaB=lambdaB, weight.factor = 1, mplx.rate = 0, num.ctrials=1, start.time=start.time, end.time=end.time,bw=bw)
      else if (control.type == "Blike") added.spikes<-gen_synth_binwise(lambdaA=lambdaA, lambdaB=lambdaB, weight.factor = 0, mplx.rate = 0, num.ctrials=1, start.time=start.time, end.time=end.time,bw=bw)
      else if (control.type == "outside") added.spikes<-gen_synth_binwise(lambdaA=lambdaA*1.25, lambdaB=lambdaB, weight.factor = 1, mplx.rate = 0, num.ctrials=1, start.time=start.time, end.time=end.time,bw=bw)
      else if (control.type == "switch") added.spikes<-gen_synth_binwise(lambdaA=lambdaA, lambdaB=lambdaB, weight.factor = .5, mplx.rate = 0, num.ctrials=1, start.time=start.time, end.time=end.time,bw=bw)
      else if (control.type == "average") added.spikes<-gen_synth_binwise(lambdaA=(lambdaA+lambdaB)/2, mplx.rate=0, num.ctrials=1,start.time=start.time, end.time=end.time, bw=bw)
      else if (control.type == "wt_average") added.spikes<-gen_synth_binwise(lambdaA=(.75*lambdaA+.25*lambdaB),mplx.rate=0, num.ctrials=1,start.time=start.time, end.time=end.time, bw=bw)
      else if (control.type == "avg_mplxd") added.spikes<-gen_synth_binwise(lambdaA=lambdaA, lambdaB=lambdaB,mplx.rate = .5, num.ctrials=1,start.time=start.time, end.time=end.time, bw = bw)
      else if (control.type == "wt_avg_mplxd") added.spikes<-gen_synth_binwise(lambdaA=lambdaA, lambdaB=lambdaB,mplx.rate = .75, num.ctrials=1,start.time=start.time, end.time=end.time, bw = bw)
      else if (control.type == "weak_switch85") added.spikes<-gen_synth_binwise(lambdaA=lambdaA*.85, lambdaB=lambdaB*1.15, weight.factor = .5, mplx.rate = 0, num.ctrials=1, start.time=start.time, end.time=end.time,bw=bw)
      else if (control.type == "weak_switch90") added.spikes<-gen_synth_binwise(lambdaA=lambdaA*.90, lambdaB=lambdaB*1.1, weight.factor = .5, mplx.rate = 0, num.ctrials=1, start.time=start.time, end.time=end.time,bw=bw)
      else if (control.type == "weak_switch95") added.spikes<-gen_synth_binwise(lambdaA=lambdaA*.95, lambdaB=lambdaB*1.05, weight.factor = .5, mplx.rate = 0, num.ctrials=1, start.time=start.time, end.time=end.time,bw=bw)
      else if (control.type == "weak_switch80") added.spikes<-gen_synth_binwise(lambdaA=lambdaA*.80, lambdaB=lambdaB*1.20, weight.factor = .5, mplx.rate = 0, num.ctrials=1, start.time=start.time, end.time=end.time,bw=bw)
      #else if (control.type == "smoo_mplxd") added.spikes <- smooth_mplxd_wrap(ntrials = 1, lambdaA = lambdaA, lambdaB=lambdaB, time.bins = seq(from = start.time/1000,
                                                                                                                                                            #  to = end.time/1000,
                                                                                                                                                            #  by = 0.001))
      else added.spikes <- 0
      added.spikes[,1]<-ABindex
      if (exists('ABspikes')) ABspikes<-rbind(ABspikes,added.spikes)
      else ABspikes<-added.spikes
    }
    }
  }
  ABspikes
  
}

rPoiProc <- function(lambda, time.bins = seq(0, 1, .001)){
  ## time measured in seconds. lambda is a vector of values of 
  ## a rate function (in Hz) recorded at the mid-points of time.bins
  
  bin.widths <- diff(time.bins)
  nbins <- length(bin.widths)
  
  bins.with.spike <- which(runif(nbins) < lambda*bin.widths)
  return(runif(length(bins.with.spike), time.bins[bins.with.spike], time.bins[1+bins.with.spike]))
}

smooth_mplxd_wrap <- function(ntrials, lambdaA, lambdaB, time.bins = seq(0, 1, 0.001)){
  jit <- 0.001*sin.period/length(ntrials)
  mid.points <- (time.bins[-1] + time.bins[-length(time.bins)])/2
  lambda <- ((sin((mid.points+jit)/sin.period*2*pi*1000)+1)/2*(lambdaA-lambdaB))+lambdaB
  spiketimes1 <- rPoiProc(lambda, time.bins)
  
  all_trials <- cbind(ntrials[1], 1000*spiketimes1)
  if (length(ntrials) > 1){
    for (trial_num in ntrials[-1]){
      jit <- 0.001*(trial_num - ntrials[1])*sin.period/length(ntrials)
      lambda <- ((sin((mid.points+jit)/sin.period*2*pi*1000)+1)/2*(lambdaA-lambdaB))+lambdaB
      new_spikes <- rPoiProc(lambda, time.bins)
      temp_mat <- cbind(trial_num, 1000*new_spikes)
      all_trials <- rbind(all_trials, temp_mat)
    }
  }
  return(all_trials)
}