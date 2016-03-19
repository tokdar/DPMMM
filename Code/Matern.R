#Matern
rm(list=ls())
library(MASS)
Matern<-function(nu,l,r){
  A = 2^(1-nu)/gamma(nu)
  B = (sqrt(2*nu)*r/l)^nu
  C = besselK(sqrt(2*nu)*r/l,nu)
  return(A*B*C)
}

SE<-function(l,r) return( exp(-r^2/(2*l^2) )  )
x = seq(0,10,len = 1000)

l = 5
plot(x, Matern(2,l,x), type="l", lty = 2, col = "red" )
lines(x, Matern(.5,l,x), type="l", lty=1, col = "blue")
lines(x, Matern(7/2,l,x), type="l", lty=2, col = "green")

df = data.frame(x = rep(x,times = 3), Cov = c( Matern(.5,l,x), Matern(2,l,x), Matern(7/2,l,x)))
df$Nu = factor( rep( c("1/2", "2", "7/2"), each = 1000) )
pdf("~/Dropbox/Neuro/Figures/Matern_Kernel.pdf")
g<-ggplot(df)+
  geom_line(aes(x,Cov, color = Nu))+
  xlab("|t-s|")+
  ylab("Covariance")+
  theme(axis.text=element_text(size=20, color="black"),axis.title=element_text(size=24,face="bold"), legend.text=element_text(size=20))
print(g)
dev.off()

t.T = 40
l = 5
nu = 7/2
sigma2 = 1
C = matrix(nrow = t.T, ncol = t.T)
for(t in 1:t.T){
  for(s in 1:t.T){
    r = abs(t-s)
    C[t,s] = sigma2*ifelse(r==0, 1, Matern(nu,l,r) ) 
  }
}

eta_bar = matrix(rep(0,t.T), ncol = 1)
nSamples = 100
eta = matrix(nrow = nSamples, ncol = t.T)
for(i in 1:nSamples) eta[i,] = mvrnorm(1,eta_bar, C)
alpha = 1/(1+exp(-eta))

m.1 = apply(eta,2,mean); plot(m.1, ylim=c(-.2,.2))
v.1 = apply(eta,2,var); plot(v.1, ylim=c(.4,.6))


df = data.frame(Time = 1:t.T, alpha1 = eta[1,], alpha2 = eta[2,], alpha3 = eta[3,] )
pdf("~/Dropbox/Neuro/Figures/Eta_Trajectories.pdf")
ggplot(df)+
  #geom_point(aes(Time,alpha1))+
  #geom_line(aes(Time,alpha1))+
  geom_point(aes(Time,alpha2), color = "black")+
  geom_line(aes(Time,alpha2), color = "black")+
  geom_point(aes(Time,alpha3), color = "black", alpha = .5)+
  geom_line(aes(Time,alpha3), color = "black", alpha = .5)+
  ylab("eta")+
  theme(axis.text=element_text(size=20, color="black"),axis.title=element_text(size=24,face="bold"), legend.text=element_text(size=20))
dev.off()


df = data.frame(Time = 1:t.T, alpha1 = alpha[1,], alpha2 = alpha[2,], alpha3 = alpha[3,] )
pdf("~/Dropbox/Neuro/Figures/Alpha_Trajectories.pdf")
ggplot(df)+
  #geom_point(aes(Time,alpha1))+
  #geom_line(aes(Time,alpha1))+
  geom_point(aes(Time,alpha2), color = "black")+
  geom_line(aes(Time,alpha2), color = "black")+
  geom_point(aes(Time,alpha3), color = "black", alpha = .5)+
  geom_line(aes(Time,alpha3), color = "black", alpha = .5)+
  ylim(0,1)+
  ylab("alpha")+
  theme(axis.text=element_text(size=20, color="black"),axis.title=element_text(size=24,face="bold"), legend.text=element_text(size=20))
dev.off()


#I think the Matern is the way to go.  


SE<-function(l,r) return( exp(-r^2/(2*l^2) )  )
x = seq(0,20,len = 1000)

df = data.frame(x = rep(x,times = 4), l1 = SE(1,x), l3 = SE(3,x), l5 = SE(5,x), l15 = SE(15,x)  )
pdf("~/Dropbox/Neuro/Figures/SE_Kernel.pdf")
g<-ggplot(df)+
  geom_line(aes(x,l1, color = "1"), size=2)+
  geom_line(aes(x,l3, color = "3"), size=2)+
  geom_line(aes(x,l5, color = "5"), size=2)+
  geom_line(aes(x,l15, color = "15"), size=2)+
  xlab("|t-s|")+
  ylab("Covariance")+
  scale_colour_manual(name = "length", breaks = c("1","3","5","15"), values =c("black", "red","blue","green"))+
  theme(axis.text=element_text(size=20, color="black"),axis.title=element_text(size=24,face="bold"), legend.text=element_text(size=20))
print(g)
dev.off()