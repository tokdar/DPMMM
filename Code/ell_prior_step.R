ell_prior_step<-function(ell_index,ell_0,L){
  nEll = sapply(1:L, function(l) sum(ell_index==l))
  concentration = ell_0 + nEll 
  ell_prior = rDirichlet(1,concentration)
  return(ell_prior)
}
