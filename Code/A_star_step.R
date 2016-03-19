A_star_step<-function(A,lambda_A,alpha){
  nRep = dim(alpha)[1]
  t.T = dim(alpha)[2]
  A_star = matrix(nrow = nRep, ncol = t.T)
  for(i in 1:nRep){
    for(t in 1:t.T){
      A_star[i,t] = A[i,t] + rpois(1,lambda_A[t]*(1-alpha[i,t]))
    }
  }
  return(A_star)
}