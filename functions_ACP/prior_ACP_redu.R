#-----------------------------------------------------------------------------#
# This function first elicits the asymmetric conjugate prior on the reduced-  #
# form parameterization and then constructs the implied prior on the          #
# structural parameterization                                                 #
#                                                                             #
# See: Chan, J.C.C. (2022) Asymmetric conjugate priors for large Bayesian     #
#         VARs, Quantitative Economics, Vol. 13, pp. 1145-1169.               #
#                                                                             #
# Replication of Matlab Codes by J.C.C. Chan by Maximilian Boeck              #
# 07/11/22                                                                    #
#-----------------------------------------------------------------------------#
prior_ACP_redu <- function(n,plag,kappa,sig2,deterministics,idx_ns){
  
  Ki     = n*plag+ifelse(deterministics$cons,1,0)+ifelse(deterministics$trend,1,0)+ifelse(deterministics$quadratictrend,1,0)
  k_beta = n*Ki
  prior_stru = prior_ACP_stru(n,plag,kappa,sig2,deterministics,idx_ns)
  
  prior <- list()
  prior[["alp0"]]  = prior_stru[["alp0"]]
  prior[["beta0"]] = prior_stru[["beta0"]]
  prior[["Valp"]]  = prior_stru[["Valp"]]
  prior[["nu"]]    = prior_stru[["nu"]]
  prior[["S"]]     = prior_stru[["S"]]
  prior[["Vbeta"]] = matrix(0, k_beta/n, n)
  for(nn in 1:n){
    for(jj in 1:Ki){
      if(nn == 1){
        prior[["Vbeta"]][jj,nn] = prior_stru[["Vbeta"]][jj,nn]
      }else{
        prior[["Vbeta"]][jj,nn] = prior_stru[["Vbeta"]][jj,nn] + 
          sum(
            prior_stru[["Vbeta"]][jj,1:(nn-1)] + prior_stru[["beta0"]][jj,1:(nn-1)]^2/sig2[1:(nn-1),1]
          )
      }
    }
  }
  
  return(prior)
}
