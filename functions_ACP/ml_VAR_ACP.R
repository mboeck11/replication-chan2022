#-----------------------------------------------------------------------------#
# This function computes the marginal likelihood value under the asymmetric   #
# conjugate prior                                                             #
#                                                                             #
# See: Chan, J.C.C. (2022) Asymmetric conjugate priors for large Bayesian     #
#         VARs, Quantitative Economics, Vol. 13, pp. 1145-1169.               #
#                                                                             #
# Replication of Matlab Codes by J.C.C. Chan by Maximilian Boeck              #
# 07/11/22                                                                    #
#-----------------------------------------------------------------------------#
ml_VAR_ACP <- function(Y,X,plag,prior){
  # dimensions
  bigT = nrow(Y)
  n    = ncol(Y)
  
  lml = -n*bigT/2*log(2*pi)
  count_alp = 0
  for(nn in 1:n){
    yi = Y[,nn,drop=FALSE]
    
    if(nn > 1){
      mi = rbind(prior[["beta0"]][,nn,drop=FALSE],prior[["alp0"]][(count_alp+1):(count_alp+nn-1),1,drop=FALSE])
      Vi = diag(c(prior[["Vbeta"]][,nn],prior[["Valp"]][(count_alp+1):(count_alp+nn-1),1]))
      
      Xi = cbind(X, -Y[,1:(nn-1)])
    }else{
      mi = prior[["beta0"]][,nn]
      Vi = diag(prior[["Vbeta"]][,nn])
      
      Xi = X
    }
    Si  = prior[["S"]][nn,1]
    nui = prior[["nu"]][nn,1]
    
    # iVi = solve(Vi)
    iVi         = diag(1/diag(Vi))
    Kthetai     = iVi + crossprod(Xi)
    CKthetai    = t(chol(Kthetai))
    thetai_hat  = solve(t(CKthetai)) %*% (solve(CKthetai) %*% (iVi %*% mi + crossprod(Xi, yi)))
    Si_hat      = Si + (crossprod(yi) + t(mi)%*%iVi%*%mi - t(thetai_hat)%*%Kthetai%*%thetai_hat)/2
    
    lml = lml - 1/2*(sum(log(diag(Vi))) + 2*sum(log(diag(CKthetai)))) + nui*log(Si) -
      (nui+bigT/2)*log(Si_hat) + lgamma(nui+bigT/2) - lgamma(nui)
    
    count_alp = count_alp + nn - 1
  }
  
  return(lml)
}
