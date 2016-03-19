K_Matern<-function(nu,l,r){
  A = 2^(1-nu)/gamma(nu)
  B = (sqrt(2*nu)*r/l)^nu
  C = besselK(sqrt(2*nu)*r/l,nu)
  return(A*B*C)
}