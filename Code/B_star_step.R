B_star_step<-function(B,lambda_B,alpha){
  nRep = dim(alpha)[1]
  t.T = dim(alpha)[2]
  B_star = matrix(nrow = nRep, ncol = t.T)
  for(i in 1:nRep){
    for(t in 1:t.T){
      B_star[i,t] = B[i,t] + rpois(1, lambda_B[t]*alpha[i,t] )
    }
  }
  return(B_star)
}