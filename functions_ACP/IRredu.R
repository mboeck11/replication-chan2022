#-----------------------------------------------------------------------------#
# This function computes impulse responses                                    #
#                                                                             #
# See: Chan, J.C.C. (2022) Asymmetric conjugate priors for large Bayesian     #
#         VARs, Quantitative Economics, Vol. 13, pp. 1145-1169.               #
#                                                                             #
# Replication of Matlab Codes by J.C.C. Chan by Maximilian Boeck              #
# 07/11/22                                                                    #
#-----------------------------------------------------------------------------#
IRredu <- function(A,L,nstep,nshock){
  # dimensions
  n    = ncol(A)
  plag = nrow(A)/n
  
  temp = gen_compMat(A, n, plag)
  Cm   = temp$Cm
  Jm   = temp$Jm
  
  impresp = array(NA_real_, c(n, nshock, nstep))
  impresp[,,1] <- L[,1:nshock]
  Cmi <- Cm
  for(ihor in 2:nstep){
    impresp[,,ihor] <- t(Jm) %*% Cmi %*% Jm %*% L[,1:nshock]
    Cmi <- Cmi %*% Cm
  }
  
  return(impresp)
}
