
## generate spikes in a binwise poisson fashion
gen_synth_binwise <- function(lambdaA,lambdaB=NULL, weight.factor = 1, mplx.rate=0 ,num.ctrials,start.time=0, end.time=1000, bw=(end.time-start.time)){
   #returns c counts after manufacturing Cs using poisson spike generator with input rate
   #pass in a vector of lambdas.
   #weight.factor is trial level, mplx.rate is within trials
   spikes <- matrix(0,num.ctrials,end.time-start.time)

   for (i in 1:num.ctrials){
      if (runif(1,0,1)<=weight.factor){
         lambda<-lambdaA #randomly picking distribution to draw from based on weight factor
         other.lambda <- lambdaB
      } 
      else {
         lambda <- lambdaB
         other.lambda <- lambdaA
      }
      if (mplx.rate !=0){
         for (bin in 1:length(lambda)){
            if (runif(1,0,1)<mplx.rate)lambda[bin] <- other.lambda[bin] #picking bin rate based on mplx rate
         }
      }

      spikes[i,]<- c(sapply(lambda,makeSpikes,start=1, end=bw))
      
   }
   #manipulating format
   synth.times <-which(spikes==1,arr.ind=TRUE)
   synth.times <- synth.times[order(synth.times[,1]),]
   added.spikes <- synth.times+start.time
   
}

gen_smooth_mplxd <- function(lambdaA,lambdaB=NULL ,num.ctrials,start.time=0, end.time=1000, bw=(end.time-start.time)){
  #handling mplx differently by smoothly fluctuating from A to B
  spikes <- matrix(0,num.ctrials,end.time-start.time)
  
  for (i in 1:num.ctrials){

    max_rate = lambdaA
    min_rate = lambdaB
    jit = sample(1:200,1)
    t=seq(start.time:(end.time-1))
    
    # little ugly, but making the "rate" constantly fluctuate between max and min values
    # BW sets the period of these fluctuations, in this case 200
    lambda = ((sin((t+jit)/bw*2*pi)+1)/2*(max_rate-min_rate))+min_rate 
    
    
    spikes[i,]<- c(sapply(lambda,makeSpikes,start=1, end=1)) #1 ms bins for spike generation
    
  }
  #manipulating format
  synth.times <-which(spikes==1,arr.ind=TRUE)
  synth.times <- synth.times[order(synth.times[,1]),]
  added.spikes <- synth.times+start.time
  
}


## poisson spike generator
makeSpikes <- function(rate=50,start=0,end=1000, numtrains = 1){
   times <- seq(start,end)
   spikes <- matrix(0,numtrains,length(times))
   for(train in 1:numtrains){
      vt <- runif(length(times),0,1)
      spikes[train,] <- vt < (rate*.001)
   }
   spikes
}