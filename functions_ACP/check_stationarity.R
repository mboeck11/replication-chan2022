check_stationarity <- function(Yraw, plag){
  n = ncol(Yraw)
  
  # sample variance of AR(plag) process
  idx_ns  <- NULL
  for(nn in 1:n){
    Ylag_nn        = mlag(Yraw[,nn],plag)
    Ylag_nn        = Ylag_nn[(plag+1):nrow(Ylag_nn),,drop=FALSE]
    Y_nn           = Yraw[(plag+1):nrow(Yraw),nn,drop=FALSE]
    Ylag_nn        = cbind(Ylag_nn,seq(1,nrow(Y_nn)))
    alpha_nn       = solve(crossprod(Ylag_nn))%*%crossprod(Ylag_nn,Y_nn)
    Cm             = matrix(0,plag,plag)
    Cm[1,]         = alpha_nn[1:plag,1]
    Cm[2:plag,1:(plag-1)] = diag(plag-1)
    eigs           = max(Re(eigen(Cm)$values))
    if(eigs>1) idx_ns = c(idx_ns, nn)
  }
  
  return(idx_ns)
}
