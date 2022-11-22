#-----------------------------------------------------------------------------#
# This function obtains the optimal shrinkage hyperparameter values and the   #
# associated marginal likelihood under the asymmetric conjugate prior         #
#                                                                             #
# See: Chan, J.C.C. (2022) Asymmetric conjugate priors for large Bayesian     #
#         VARs, Quantitative Economics, Vol. 13, pp. 1145-1169.               #
#                                                                             #
# Replication of Matlab Codes by J.C.C. Chan by Maximilian Boeck              #
# 07/11/22                                                                    #
#-----------------------------------------------------------------------------#
get_OptKappa <- function(Y,X,plag,kappa,sig2,deterministics,type,idx_ns){
  
  if(type == "stru"){
    f = function(k) -ml_VAR_ACP(Y, X, plag, 
                                prior_ACP_stru(n, plag, c(exp(k[1]), exp(k[2]), kappa[3], kappa[4]), sig2, deterministics, idx_ns))
  }else if(type == "redu"){
    f = function(k) -ml_VAR_ACP(Y, X, plag, 
                                prior_ACP_redu(n, plag, c(exp(k[1]), exp(k[2]), kappa[3], kappa[4]), sig2, deterministics, idx_ns))
  }
  minimizer = optim(par=log(kappa[1:2]), fn=f, method="BFGS")
  
  ml_opt    = -minimizer$value
  kappa_opt = c(exp(minimizer$par), kappa[3:4])
  
  return(list(ml_opt=ml_opt, kappa_opt=kappa_opt))
}
