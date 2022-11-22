#-----------------------------------------------------------------------------#
# This function first directly elicits the asymmetric conjugate prior on the  #
# structural parameterization                                                 #
#                                                                             #
# See: Chan, J.C.C. (2022) Asymmetric conjugate priors for large Bayesian     #
#         VARs, Quantitative Economics, Vol. 13, pp. 1145-1169.               #
#                                                                             #
# Replication of Matlab Codes by J.C.C. Chan by Maximilian Boeck              #
# 07/11/22                                                                    #
#-----------------------------------------------------------------------------#
prior_ACP_stru <- function(n,plag,kappa,sig2,deterministics,idx_ns){
  
  # dimensions
  k_beta = n*(n*plag+ifelse(deterministics$cons,1,0)+ifelse(deterministics$trend,1,0)+ifelse(deterministics$quadratictrend,1,0))
  k_alp  = n*(n-1)/2
  prior = list()
  prior[["beta0"]] = matrix(0,k_beta/n,n)
  prior[["alp0"]]  = matrix(0,k_alp,1)
  prior[["Vbeta"]] = matrix(0,k_beta/n,n)
  prior[["Valp"]]  = matrix(0,k_alp,1)
  prior[["nu"]]    = matrix(0,n,1)
  prior[["S"]]     = matrix(0,n,1)
  
  count_alp = 0
  for(nn in 1:n){
    is_ns = any(idx_ns == nn)
    temp = prior_ACPi(n, plag, nn, kappa, sig2, is_ns, deterministics)
    prior[["beta0"]][,nn] = temp$mi_beta
    prior[["Vbeta"]][,nn] = temp$Vi_beta
    if(nn>1) prior[["alp0"]][(count_alp+1):(count_alp+nn-1)] = temp$mi_alp
    if(nn>1) prior[["Valp"]][(count_alp+1):(count_alp+nn-1)] = temp$Vi_alp
    prior[["nu"]][nn,1] = temp$nui
    prior[["S"]][nn,1]  = temp$Si
    count_alp = count_alp + nn - 1
  }
  
  return(prior)
}

prior_ACPi <- function(n, plag, nn, kappa, sig2, is_ns, deterministics){
  ki_beta = n*plag + ifelse(deterministics$cons,1,0)+ifelse(deterministics$trend,1,0)+ifelse(deterministics$quadratictrend,1,0)
  ki_alp  = nn-1
  mi_beta = matrix(0, ki_beta, 1)
  Vi_beta = matrix(0, ki_beta, 1)
  
  # construct Vi_beta
  for(jj in 1:ki_beta){
    if(jj <= (n*plag)){
      l   = ceiling(jj/n)
      idx = jj %% n
      if(idx == 0){
        idx = n
      }
    }else{
      idx = jj
    }
    
    # own lags
    if(idx == nn){
      Vi_beta[jj,1] = kappa[1]/(l^2*sig2[idx,1])
      if(l == 1 && is_ns) mi_beta[jj,1] = 1
    }else if(jj < n*plag+1){
      Vi_beta[jj,1] = kappa[2]/(l^2*sig2[idx,1])
    }else if(jj > n*plag){
      Vi_beta[jj,1] = kappa[4]
    }
  }
  
  # alpha
  if(ki_alp > 0){
    mi_alp  = matrix(0, ki_alp, 1)
    Vi_alp  = matrix(0, ki_alp, 1)
    
    # construct Vi_alp
    for(jj in 1:ki_alp){
      Vi_alp[jj,1] = kappa[3]/sig2[jj]
    }
  }else{
    mi_alp = Vi_alp = NULL
  }
  
  Si  = sig2[nn]/2
  nui = 1 + nn/2
  
  return(list(mi_beta=mi_beta, Vi_beta=Vi_beta,
              mi_alp=mi_alp, Vi_alp=Vi_alp,
              Si=Si, nui=nui))
}











