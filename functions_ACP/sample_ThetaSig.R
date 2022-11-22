#-----------------------------------------------------------------------------#
# This function obtains the posterior draws of alpha, beta, and Sigma under   #
# the asymmetric conjugate prior                                              #
#                                                                             #
# See: Chan, J.C.C. (2022) Asymmetric conjugate priors for large Bayesian     #
#         VARs, Quantitative Economics, Vol. 13, pp. 1145-1169.               #
#                                                                             #
# Replication of Matlab Codes by J.C.C. Chan by Maximilian Boeck              #
# 08/11/22                                                                    #
#-----------------------------------------------------------------------------#
sample_ThetaSig <- function(Y,X,plag,prior,deterministics,nsim){
  # dimensions
  bigT    = nrow(Y)
  n       = ncol(Y)
  k_betai = n*plag + ifelse(deterministics$cons,1,0)+ifelse(deterministics$trend,1,0)+ifelse(deterministics$quadratictrend,1,0)
  k_beta  = n*k_betai
  k_alp   = n*(n-1)/2
  
  # container
  Beta = matrix(0, nsim, k_beta)
  Alp  = matrix(0, nsim, k_alp)
  Sig  = matrix(0, nsim, n)
  
  count_alp = 0
  for(nn in 1:n){
    yi = Y[,nn,drop=FALSE]
    if(nn > 1){
      ki = k_betai + nn-1
      mi = rbind(prior[["beta0"]][,nn,drop=FALSE],prior[["alp0"]][(count_alp+1):(count_alp+nn-1),1,drop=FALSE])
      Vi = diag(c(prior[["Vbeta"]][,nn],prior[["Valp"]][(count_alp+1):(count_alp+nn-1),1]))
      
      Xi = cbind(X, -Y[,1:(nn-1)])
    }else{
      ki = k_betai
      mi = prior[["beta0"]][,nn]
      Vi = diag(prior[["Vbeta"]][,nn])
      
      Xi = X
    }
    Si  = prior[["S"]][nn,1]
    nui = prior[["nu"]][nn,1]
    
    # compute the parameters of the posterior distribution
    iVi         = diag(1/diag(Vi))
    Kthetai     = iVi + crossprod(Xi)
    CKthetai    = t(chol(Kthetai))
    thetai_hat  = solve(t(CKthetai)) %*% (solve(CKthetai) %*% (iVi %*% mi + crossprod(Xi, yi)))
    Si_hat      = Si + (crossprod(yi) + t(mi)%*%iVi%*%mi - t(thetai_hat)%*%Kthetai%*%thetai_hat)/2
    
    # sample theta and sig
    Sigi   = 1/rgamma(nsim, shape=nui+bigT/2, scale=1/Si_hat)
    U      = t(matrix(rnorm(nsim, 0, sqrt(Sigi)), ki, nsim))
    Thetai = matrix(thetai_hat, nsim, ki, byrow=TRUE) + U %*% solve(CKthetai)
    
    # save
    Beta[,((nn-1)*k_betai+1):(k_betai*nn)] = Thetai[,1:k_betai]
    if(nn>1) Alp[,(count_alp+1):(count_alp+nn-1)]   = Thetai[,(k_betai+1):(k_betai+nn-1)]
    Sig[,nn] = Sigi
    
    count_alp = count_alp + nn - 1
  }
  
  return(list(Beta=Beta, Alp=Alp, Sig=Sig))
}
