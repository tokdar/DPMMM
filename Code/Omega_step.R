Omega_step<-function(A_star,B_star,eta){
  nRep = dim(eta)[1]
  t.T = dim(eta)[2]
  Omega = matrix(nrow = nRep, ncol = t.T)
  for(i in 1:nRep){
    for(t in 1:t.T){
      N = A_star[i,t] + B_star[i,t]
      Omega[i,t] = rpg.devroye(num=1,n = N, z = eta[i,t])
    }
  }
  return(Omega)
}