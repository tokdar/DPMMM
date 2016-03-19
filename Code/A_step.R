A_step<-function(X,lambda_A,lambda_B,alpha){
  nRep = dim(alpha)[1]
  t.T = dim(alpha)[2]
  A = matrix(nrow=nRep, ncol = t.T)
  for(i in 1:nRep){
    for(t in 1:t.T){
      rate = lambda_A[t]*alpha[i,t]/(lambda_A[t]*alpha[i,t] + lambda_B[t]*(1-alpha[i,t]))
      A[i,t] = rbinom(1,X[i,t], rate)
    }
  }
  return(A)
}