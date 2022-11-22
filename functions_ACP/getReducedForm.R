#-----------------------------------------------------------------------------#
# This function obtains the reduced-form parameters from the structural-form  #
# parameters                                                                  #
#                                                                             #
# See: Chan, J.C.C. (2022) Asymmetric conjugate priors for large Bayesian     #
#         VARs, Quantitative Economics, Vol. 13, pp. 1145-1169.               #
#                                                                             #
# Replication of Matlab Codes by J.C.C. Chan by Maximilian Boeck              #
# 08/11/22                                                                    #
#-----------------------------------------------------------------------------#
getReducedForm <- function(store_Alp, store_Beta, store_Sig, plag){
  # dimensions
  n      = ncol(store_Sig)
  nsim   = nrow(store_Sig)
  k_beta = ncol(store_Beta)
  
  A = diag(n)
  store_Btilde   = array(NA_real_, c(nsim, k_beta/n, n))
  store_Sigtilde = array(NA_real_, c(nsim, n, n))
  for(isim in 1:nsim){
    Alp  = t(store_Alp[isim,,drop=FALSE])
    Beta = t(store_Beta[isim,,drop=FALSE])
    Sig  = store_Sig[isim,,drop=TRUE]
    
    # transform the parameters into reduced-form
    count_alp = 1
    for(ii in 2:n){
      for(jj in 1:(ii-1)){
        A[ii,jj] = Alp[count_alp]
        count_alp = count_alp + 1
      }
    }
    A_inverse = solve(A)
    Sigtilde  = A_inverse %*% diag(Sig) %*% t(A_inverse)
    Btilde    = t(A_inverse %*% t(matrix(Beta, ncol=n)))
    
    store_Btilde[isim,,]   = Btilde
    store_Sigtilde[isim,,] = Sigtilde
  }
  
  return(list(store_Btilde=store_Btilde, store_Sigtilde=store_Sigtilde))
}




