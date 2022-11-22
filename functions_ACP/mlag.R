mlag <- function(X,plag){
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,plag*N)
  for (ii in 1:plag){
    Xlag[(plag+1):Traw,(N*(ii-1)+1):(N*ii)] <- X[(plag+1-ii):(Traw-ii),(1:N)]
  }
  colnames(Xlag) <- paste0(colnames(X),".lag",rep(seq(plag),each=N))
  return(Xlag)
}
