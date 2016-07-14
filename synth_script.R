#script to make the various synthetic trials

setwd("C:/Users/jtm47/Documents/code/RCode/Power")
source("synth_trial_maker.R")

#need to specifiy the trial indices for each trial type
#generate the trial info
#then run through Alike and Blike spike generation, and one AB method
#then output both sets of data into .txt files with some clear naming convention
#

control.type = c("Alike", "weak_switch80","weak_switch90", "switch", "average","wt_average","avg_mplxd","wt_avg_mplxd", "outside", "smoo_mplxd")
n_trials_vec=c(5,20,50) #number of trials to be simulated for each triplet
repeats = 100   #number of files(triplets) to be generated for each condition
lambdaASeed=50  #Note: the spike generator is made to accept lambda in a vector form with lambda for each time bin,
lambdaBSeed=20  #these seeds are used to make that vector repeats of a single rate
bw=200       #Note:can not run correctly with only one bin, and total trial time must be evenly divisible by bin size
             #also: bw is only use in "mplxd" controls which have variable rates within trials
start.time=-200   #Note:these are used only for the spike generation function
end.time=1000     

for (n_trials in n_trials_vec){
   for (this.control in control.type) {
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
         write(t(trial_data),file=paste("~/Data/Synth_Data/",outfile,".txt",sep=""),ncolumns=10,sep="\t")
         
         spikeout <- paste("~/Data/Synth_Data/",outfile,"_spiketimes.txt",sep="")
         write(t(spike_data),file=spikeout,ncolumns=2,sep="\t")
         
        if (exists('file.list')) file.list<-rbind(file.list,outfile)
        else file.list<-outfile
         
      }
   }
}

write(file.list,file="Complete_List.txt")





