p_step<-function(zeta,a,b){
  nRep = length(zeta)
  nZ = sum(zeta)
  a.tilde = a + nZ
  b.tilde = b + (nRep-nZ)
  p = rbeta(1,a.tilde,b.tilde)
  return(p)
}