#-----------------------------------------------------------------------------#
# This function implements the estimation of the Bayesian VAR proposed in     #
# Chan (2022) featuring asymmetric conjugate priors                           #
#                                                                             #
# See: Chan, J.C.C. (2022) Asymmetric conjugate priors for large Bayesian     #
#         VARs, Quantitative Economics, Vol. 13, pp. 1145-1169.               #
#                                                                             #
# Replication of Matlab Codes by J.C.C. Chan by Maximilian Boeck              #
# 07/11/22                                                                    #
#-----------------------------------------------------------------------------#
bvar_ACP <- function(Yraw, plag = 1, args = NULL){
  #----------------------------------------INPUTS----------------------------------------------------#
  # prepare arguments
  draws = burnin = 5000; cons=TRUE; trend=FALSE; quadratictrend=FALSE; thin=1; Ex=NULL; save.prior=FALSE; optimal.kappa=FALSE
  if(!is.null(args)){
    for(aa in c("draws","burnin","cons","trend","quadratictrend","thin","Ex","save.prior","optimal.kappa")){
      if(aa%in%names(args)) assign(aa, args[[aa]])
    }
  }
  arglist=list(Yraw=Yraw, plag=plag, draws=draws, burnin=burnin, cons=cons, trend=trend, quadratictrend=quadratictrend,
               thin=thin, Ex=Ex, save.prior=save.prior)
  deterministics=list(cons=cons, trend=trend, quadratictrend=quadratictrend)
  
  #----------------------------------------PACKAGES--------------------------------------------------#
  require(MASS, quietly=TRUE)
  
  #-------------------------------------------START--------------------------------------------------#
  Traw  <- nrow(Yraw)
  n     <- ncol(Yraw)
  K     <- n*plag
  Ylag  <- mlag(Yraw,plag)
  varNames <- colnames(Yraw)
  if(is.null(varNames)) varNames <- paste0("Y",seq(n))
  varNameslags <- NULL
  for(pp in 1:plag) varNameslags <- c(varNameslags,paste(varNames,".lag",pp,sep=""))
  colnames(Ylag) <- varNameslags
  
  texo <- FALSE; Mex <- 0; Exraw <- NULL
  if(!is.null(Ex)){
    Exraw <- Ex; Mex <- ncol(Exraw)
    texo <- TRUE
    ExNames <- colnames(Exraw)
    if(is.null(ExNames)) ExNames <- rep("Tex",Mex)
    varNameslags <- c(varNameslags, ExNames)
  }
  
  X <- cbind(Ylag,Exraw)
  X <- X[(plag+1):nrow(X),,drop=FALSE]
  Y <- Yraw[(plag+1):Traw,,drop=FALSE]
  bigT  <- nrow(X)
  
  if(cons){
    X <- cbind(X,1)
    varNameslags <- c(varNameslags,"cons")
    colnames(X)[ncol(X)] <- "cons"
  }
  if(trend){
    X <- cbind(X,seq(1,bigT))
    varNameslags <- c(varNameslags,"trend")
    colnames(X)[ncol(X)] <- "trend"
  }
  if(quadratictrend){
    X <- cbind(X,seq(1,bigT)^2)
    varNameslags <- c(varNameslags,"quadratictrend")
    colnames(X)[ncol(X)] <- "quadratictrend"
  }
  
  k = ncol(X)
  v = (n*(n-1))/2
  
  #---------------------------------------------------------------------------------------------------------
  # PRELIMINARIES
  #---------------------------------------------------------------------------------------------------------
  # compute AR residuals of AR(plag)
  sig2 <- get_resid_ar(Yraw, plag)
  sig2 <- matrix(c(2.5926e-05,2.9884e-06,0.1371,7.6602e-04,2.3551e-04,0.2599), n, 1) # delete later
  
  # check stationarity
  idx_ns = check_stationarity(Yraw, plag)
  
  if(optimal.kappa){
    temp_OptKappa = get_OptKappa(Y, X, plag, kappa, sig2, deterministics, type="stru", idx_ns)
    kappa = temp_OptKappa$kappa_opt
  }else{
    kappa = kappa     = c(1, 1, 1, 100)
  }
  prior_stru = prior_ACP_stru(n, plag, kappa, sig2, deterministics, idx_ns)
  prior_redu = prior_ACP_redu(n, plag, kappa, sig2, deterministics, idx_ns)
  
  #---------------------------------------------------------------------------------------------------------
  # SAMPLER MISCELLANEOUS
  #---------------------------------------------------------------------------------------------------------
  ntot  <- draws+burnin
  
  # thinning
  count <- 0
  thindraws    <- draws/thin
  thin.draws   <- seq(burnin+1,ntot,by=thin)
  arglist      <- c(arglist, thindraws=thindraws)
  
  #---------------------------------------------------------------------------------------------------------
  # STORAGES
  #---------------------------------------------------------------------------------------------------------
  A_store       <- array(NA_real_, c(thindraws, k, M))
  L_store       <- array(NA_real_, c(thindraws, M, M))
  res_store     <- array(NA_real_, c(thindraws, bigT, M))
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 1: Sample nbatch draw from posterior
    
    
    #----------------------------------------------------------------------------
    # Step 2: different shrinkage prior setups
    #----------------------------------------------------------------------------
    
    #----------------------------------------------------------------------------
    # Step 3: Sample variances
    #----------------------------------------------------------------------------
    
    #----------------------------------------------------------------------------
    # Step 4: store draws
    #----------------------------------------------------------------------------
    if(irep %in% thin.draws){
      count <- count+1
      
      
    }
    if(irep%%50==0) print(paste0("Round: ",irep))
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,varNameslags,varNames)
  if(SV){
    dimnames(SIGMA_store)=list(NULL,NULL,varNames,varNames)
  }else{
    dimnames(SIGMA_store)=list(NULL,varNames,varNames)
  }
  ret <- list(Y=Y, X=X, A=A_store, SIGMA=SIGMA_store, res=res_store, 
              L=L_store, Sv=Sv_store, pars=pars_store,
              Aprior=Aprior_store, Lprior=Lprior_store, lambda2=lambda2_store, tau=tau_store,
              args=arglist)
  return(ret)
}
