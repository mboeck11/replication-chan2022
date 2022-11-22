get_resid_ar <- function(Yraw,plag){
  # get dimensions
  n    = ncol(Yraw)
  bigT = nrow(Yraw)-plag
  
  # sample variance of AR(plag) process
  sigma_sq  <- matrix(0,n,1)
  for(nn in 1:n){
    Ylag_nn        = mlag(Yraw[,nn],plag)
    Ylag_nn        = Ylag_nn[(plag+1):nrow(Ylag_nn),,drop=FALSE]
    Y_nn           = Yraw[(plag+1):nrow(Yraw),nn,drop=FALSE]
    Ylag_nn        = cbind(Ylag_nn,1)
    alpha_nn       = solve(crossprod(Ylag_nn))%*%crossprod(Ylag_nn,Y_nn)
    sigma_sq[nn,1] = (1/bigT)*t(Y_nn-Ylag_nn%*%alpha_nn)%*%(Y_nn-Ylag_nn%*%alpha_nn)
  }
  
  return(sigma_sq)
}
