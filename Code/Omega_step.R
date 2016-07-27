Omega_step<-function(A_star,B_star,eta){
  eta.dim <- dim(eta)
  eta.len <- prod(eta.dim)
  N <- A_star + B_star
  Omega <- sapply(1:eta.len, function(j) rpg.devroye(num = 1, n = N[j], z = eta[j]))
  dim(Omega) <- eta.dim
  Omega[Omega == 0] <- 1e-12
  return(Omega)
}
