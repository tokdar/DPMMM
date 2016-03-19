#Figures
rm(list=ls())
library(ggplot2)

N = 1100
A = Aplus = B = rep(0,N)
A[11:210] = A[411:610] = A[811:1010] = 1
Aplus[1:220] = Aplus[391:630] = Aplus[791:1030] = 1

signal = c(0,1)
mode = 1
B[21]=1
for(i in 22:N){
  if( i%%20 == 0 ){
    B[i] = B[i-21]
  }else{ B[i] = B[i-1]}
}
AB = Aplus*A + (1-Aplus)*B


cols = rep("B",N)
cols[10:211] = cols[410:611] = cols[810:1011] = "A"
df = data.frame(time = 1:N, A = A, B = B, AB = AB, Signal = cols)
df$Signal = as.factor(df$Signal)

pdf("~/Dropbox/Neuro/Figures/Mux_A.pdf")
ggplot(df, aes(x=time) ) + 
  geom_line(aes(x=time,y=A), color = "blue", size = 1) + xlim(-100,1100) + xlab("Time (ms)") + 
  theme(axis.text=element_text(size=20, color="black"),axis.title=element_text(size=24,face="bold"), legend.text=element_text(size=20))
dev.off()

pdf("~/Dropbox/Neuro/Figures/Mux_B.pdf")
ggplot(df, aes(x=time) ) + 
  geom_line(aes(x=time,y=B), color = "red", size = 1) + xlim(-100,1100) + xlab("Time (ms)") + 
  theme(axis.text=element_text(size=20, color="black"),axis.title=element_text(size=24,face="bold"), legend.text=element_text(size=20))
dev.off()

pdf("~/Dropbox/Neuro/Figures/Mux_AB.pdf")
ggplot(df, aes(x=time) ) + 
  geom_line(aes(x=time,y=AB, colour = Signal, group = 1), size = 1) + xlim(-100,1100) + 
  scale_color_manual(breaks = c("A", "B"), values=c("blue", "red")) + xlab("Time (ms)") + 
  theme(axis.text=element_text(size=20, color="black"),axis.title=element_text(size=24,face="bold"), legend.text=element_text(size=20))
dev.off()

plot(A, type="l")
plot(B, type="l")
plot(AB, type="l")


