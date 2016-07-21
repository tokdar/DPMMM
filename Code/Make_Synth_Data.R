make_synth <- function(control.type, n_trials_vec, repeats,
                            lambdaAseed, lambdaBseed, out.dir,
                            list.file.names, bw = 200,
                            start.time = -200, end.time = 1000){
#script to make the various synthetic trials

source("Code/synth_trial_maker.R")
# create directory to store data
if (!file.exists(out.dir)){
  dir.create(out.dir)
}
for (n_trials in n_trials_vec){
  for (this.control in control.type) {
    curr.dir <- paste0(out.dir, this.control, "/")
    if (!file.exists(curr.dir)){
      dir.create(curr.dir)
    }
    for (rep in c(1:repeats)){
      #row 1:a, 2:b, 3:ab
      trial_nums <- matrix(seq(1:(n_trials*3)),3,n_trials,byrow=TRUE)
      
      #make trial matrix
      trial_data <- rbind(make_trial_data(trial.nums=trial_nums[1,], trial.type="A", pos=c(24,-2859)),
                          make_trial_data(trial.nums=trial_nums[2,], trial.type="B", pos=c(-6,-2859)),
                          make_trial_data(trial.nums=trial_nums[3,], trial.type="AB", pos=c(24,-6)))
      
      #make spike matrix
      rm(Aspikes,Bspikes,ABspikes,spike_data)
      lambdaA=rep(lambdaASeed,(end.time-start.time)/bw) # making lambdas constant, could use a different function for variable lambda
      lambdaB=rep(lambdaBSeed,(end.time-start.time)/bw)
      Aspikes <- make_synth_data(trial.nums=trial_nums[1,], start.time = start.time, end.time = end.time, bw = bw,lambdaA=lambdaA,lambdaB=lambdaB,control.type="Alike")
      Bspikes <- make_synth_data(trial.nums=trial_nums[2,], start.time = start.time, end.time = end.time, bw = bw,lambdaA=lambdaA,lambdaB=lambdaB,control.type="Blike")
      ABspikes <- make_synth_data(trial.nums=trial_nums[3,], start.time = start.time, end.time = end.time, bw = bw,lambdaA=lambdaA,lambdaB=lambdaB,control.type=this.control)
      spike_data<-rbind(Aspikes,Bspikes,ABspikes)
      
      #write out files
      
      outfile<- paste(this.control,"/synth-",this.control,"-N",n_trials,"-A", lambdaASeed,"-B",lambdaBSeed,"-rep",rep, sep = "")
      write(t(trial_data),file=paste(out.dir,outfile,".txt",sep=""),ncolumns=10,sep="\t")
      
      spikeout <- paste(out.dir,outfile,"_spiketimes.txt",sep="")
      write(t(spike_data),file=spikeout,ncolumns=2,sep="\t")
      
      if (exists('file.list')) {
        file.list<-rbind(file.list,outfile)
      } else {
        file.list<-outfile
      } 
    }
  }
}

write(file.list,file=list.file.names)
}

make_synth_mix <- function(control.types, n_trials_vec, repeats,
                       lambdaAseed, lambdaBseed, out.dir,
                       list.file.names, bw = 200,
                       start.time = -200, end.time = 1000){
  #script to make the various synthetic trials
  
  source("Code/synth_trial_maker.R")
  # create directory to store data
  if (!file.exists(out.dir)){
    dir.create(out.dir)
  }
  for (n_trials in n_trials_vec){
      for (rep in c(1:repeats)){
        #row 1:a, 2:b, 3:ab
        trial_nums <- matrix(seq(1:(n_trials*3)),3,n_trials,byrow=TRUE)
        trial_data <- rbind(make_trial_data(trial.nums=trial_nums[1,], trial.type="A", pos=c(24,-2859)),
                            make_trial_data(trial.nums=trial_nums[2,], trial.type="B", pos=c(-6,-2859)),
                            make_trial_data(trial.nums=trial_nums[3,], trial.type="AB", pos=c(24,-6)))
        
        lambdaA=rep(lambdaASeed,(end.time-start.time)/bw) # making lambdas constant, could use a different function for variable lambda
        lambdaB=rep(lambdaBSeed,(end.time-start.time)/bw)
        Aspikes <- make_synth_data(trial.nums=trial_nums[1,], start.time = start.time, end.time = end.time, bw = bw,lambdaA=lambdaA,lambdaB=lambdaB,control.type="Alike")
        Bspikes <- make_synth_data(trial.nums=trial_nums[2,], start.time = start.time, end.time = end.time, bw = bw,lambdaA=lambdaA,lambdaB=lambdaB,control.type="Blike")
        outfile<- paste("mixed","-N",n_trials,"-A", lambdaASeed,"-B",lambdaBSeed,"-rep",rep, sep = "")
        # this always generates an warning
        # don't know why it is necessary
        # rm(Aspikes,Bspikes,ABspikes,spike_data)
        L_ABspikes <- lapply(1:n_trials, function(x)
                                    {make_synth_data(trial.nums=trial_nums[3,x], start.time = start.time, 
                                    end.time = end.time, bw = bw,lambdaA=lambdaA,
                                    lambdaB=lambdaB,control.type=control.types[x])})
        ABspikes <- do.call(rbind, L_ABspikes)
        spike_data<-rbind(Aspikes,Bspikes,ABspikes)
        
        #write out files
        

        write(t(trial_data),file=paste(out.dir,outfile,".txt",sep=""),ncolumns=10,sep="\t")
        
        spikeout <- paste(out.dir,outfile,"_spiketimes.txt",sep="")
        write(t(spike_data),file=spikeout,ncolumns=2,sep="\t")
        
        if (exists('file.list')) {
          file.list<-rbind(file.list,outfile)
        } else {
          file.list<-outfile
        } 
      }
  }
  
  write(file.list,file=list.file.names)
}





