#-----------------------------------------------------------------------------#
# This script reproduces the empirical application in Chan (2022):            #
# estimating a 6- or 15-variable VAR identified with sign restrictions        #
#                                                                             #
# See: Chan, J.C.C. (2022) Asymmetric conjugate priors for large Bayesian     #
#         VARs, Quantitative Economics, Vol. 13, pp. 1145-1169.               #
#                                                                             #
# Replication of Matlab Codes by J.C.C. Chan by Maximilian Boeck              #
# 07/11/22                                                                    #
#-----------------------------------------------------------------------------#
setwd("/users/mboeck/dropbox/projects/replications/chan (2022) qe")
rm(list=ls())

# load packages
library(readxl)

# load functions
funs_acp <- list.files("./functions_ACP", full.names=TRUE)
for(fun in funs_acp) if(grepl("R$",fun)) source(fun)

# settings
plag = 5
dataset = 1
if(dataset == 1){
  nsim = 5000 # number of posterior draws (that satisfy all the restrictions)
}else{
  nsim = 1000 # might take a few days to get 1000 draws
}

nbatch = 50000  # number of posterior draws sampled in a batch
horizon = 36    # impulse response horizon

# load data
names = read_xlsx("./database_2019Q4.xlsx", range="B1:P1", col_names=FALSE)
data = read_xlsx("./database_2019Q4.xlsx", range="B3:P142", col_names=FALSE)
data = as.matrix(data)
colnames(data) = as.character(names)
data = ts(data, start=c(1985,1), frequency=4)

if(dataset == 1){
  # 6 variables + 5 shocks
  # shocks: supply, demand, monetary, investment, financial
  var_id = 1:6 # GDP, GDP deflator, interest rate, investment, S&P 500, spread 1
  idx_ns = c(1,2,4,5)
  
  n = length(var_id)
  
  # sign restrictions
  supply   = matrix(c(1,-1,NA,NA,1,NA),n,1)
  demand   = matrix(c(1,1,1,NA,NA,NA),n,1)
  monetary = matrix(c(1,1,-1,NA,NA,NA),n,1)
  invest   = matrix(c(1,1,1,NA,-1,NA),n,1)
  finc     = matrix(c(1,1,1,NA,1,NA),n,1)
  S        = cbind(supply, demand, monetary, invest, finc)
  m        = ncol(S) # number of shocks
  # row inequalities
  Ridx     = c(2,4,5) # column indices to which Rineq applies to
  nR       = length(Ridx)
  Rineq    = matrix(c(-1,1,1,0,0,0,0,0,0,1,-1,-1,0,0,0,0,0,0),nR,n)
}else if(dataset == 2){
  # 15 variables + 5 shocks
  # shocks: supply, demand, monetary, investment, financial
  var_id = 1:15
  idx_ns = c(1,2,4,5,10,11,12,13,15) # index for variables in levels
  
  n = length(var_id)
  
  # sign restrictions
  supply   = matrix(c(1,-1,NA,NA,1,NA,NA,NA,NA,-1,-1,NA,1,NA,1),n,1)
  demand   = matrix(c(1,1,1,NA,NA,NA,NA,NA,NA,1,1,NA,1,1,NA),n,1)
  monetary = matrix(c(1,1,-1,NA,NA,NA,NA,NA,NA,1,1,NA,1,-1,NA),n,1)
  invest   = matrix(c(1,1,1,NA,-1,NA,NA,NA,NA,1,1,NA,1,1,-1),n,1)
  finc     = matrix(c(1,1,1,NA,1,NA,NA,NA,NA,1,1,NA,1,1,1),n,1)
  S        = cbind(supply, demand, monetary, invest, finc)
  m        = ncol(S) # number of shocks
  # row inequalities
  Ridx     = c(2,4,5) # column indices to which Rineq applies to
  nR       = length(Ridx)
  Rineq    = matrix(c(-1,1,1,rep(0,6),1,-1,-1,rep(0,33)), nR, n)
}

Y0 = data[1:8,var_id] # save first eight observations as the initial condition
Y  = data[9:nrow(data), var_id]

# settings
deterministics=list(cons=TRUE, trend=FALSE, quadratictrend=FALSE)

# dimensions
bigT  = nrow(Y)
n     = ncol(Y)
K     = n*plag

# create X
X  <- mlag(rbind(Y0[(8-plag+1):8,],Y),plag)
X  <- X[(plag+1):nrow(X),,drop=FALSE]

if(deterministics$cons)  X <- cbind(X,1)
if(deterministics$trend) X <- cbind(X,seq(1,bigT))
if(deterministics$quadratictrend) X <- cbind(X,seq(1,bigT)^2)

k = ncol(X)
v = (n*(n-1))/2

# set up prior
sig2 <- get_resid_ar(rbind(Y0[5:8,],Y), 4) # AR(4)

# check stationarity
idx_ns = check_stationarity(Y, plag)

kappa <- c(1,1,1,100)
if(dataset == 2){
  temp_OptKappa = get_OptKappa(Y, X, plag, kappa, sig2, deterministics, type="redu", idx_ns)
  kappa = temp_OptKappa$kappa_opt
}else{
  kappa = kappa     = c(1, 1, 1, 100)
}
prior_stru = prior_ACP_stru(n, plag, kappa, sig2, deterministics, idx_ns)
prior_redu = prior_ACP_redu(n, plag, kappa, sig2, deterministics, idx_ns)

start_time  = Sys.time()
count_sat   = 0 # counter for number of draws that satisfy all the conditions
count_total = 0 # counter for total number of draws
store_response = array(NA_real_, c(nsim, n, m, horizon))
cat(paste0("Computing impulse responses from a ", n,"-variable VAR"))
cat("\t to an one-standard deviation financial shock...")
while(count_sat < nsim){
  # sample nbatch draws from the posterior
  temp_draw  = sample_ThetaSig(Y, X, plag, prior_redu, deterministics, nbatch)
  store_Alp  = temp_draw$Alp
  store_Beta = temp_draw$Beta
  store_Sig  = temp_draw$Sig
  
  count_total = count_total + nbatch
  
  # obtain the reduced-form parameters
  temp_store = getReducedForm(store_Alp, store_Beta, store_Sig, plag)
  store_Btilde   = temp_store$store_Btilde
  store_Sigtilde = temp_store$store_Sigtilde
  
  for(isim in 1:nbatch){ # go through the nbatch posterior draws to find those that fulfill the sign-restrictions
    Sigtilde = store_Sigtilde[isim,,]
    msat = 0 # counter for the number of shocks that satisfies the sign-restrictions
    nRsat = 0 # counter for the number of satisfied row inequalities
    
    randMat = matrix(rnorm(n^2),n,n)
    QR      = qr(randMat)
    Q       = qr.Q(QR)
    L0      = t(chol(Sigtilde))
    L       = L0 %*% Q
    for(mm in 1:m){
      idx   = which(!is.na(S[,mm]))
      nidx  = length(idx)
      signL = sign(L[idx,])
      # check if the mm-th column satisfies the sign-restrictions
      if(sum(signL[,mm] == S[idx,mm]) == nidx){
        msat = msat + 1
      }else if(sum(signL[,mm] == -S[idx,mm]) == nidx){
        L[,mm] = -L[,mm] # change the sign of the mm-th column
        msat = msat + 1
      }else{
        break
      }
    } # END-for checking sign-restrictions
    
    # check row inequalities
    for(jj in 1:nR){
      if(Rineq[jj,,drop=FALSE] %*% L[,Ridx[jj],drop=FALSE] < 0){
        nRsat = nRsat + 1
      }else{
        nRsat = 0
        break
      }
    } # END-for checking row inequalities
    
    if(msat == m && nRsat == nR){
      count_sat = count_sat + 1
      response = IRredu(store_Btilde[isim,1:(n*plag),], L, horizon, m)
      store_response[count_sat,,,] = response
    }
  } # END-for sign-restrictions
  
  # give some output to console
  if(count_total %% 1e+6 == 0){
    cat(paste0("Out of ",round(count_total/1e+6)," million posterior draws ", count_sat, " satisfy all the restrictions."))
    cat("\n")
  }
}
cat(paste0(nsim, " posterior draws that satisfy all the restrictions are obtained."))
cat(paste0("The simulation took ", format(Sys.time()-start_time)))

irf <- apply(store_response, c(2,3,4), quantile, c(.05, .10, .16, .50, .84, .90, .95))

# graphic stuff
shock  = 5 # financial shock
ylim_u = c(.008, .003, .4, .03, .015, .2)
ylim_l = c(-.001, -.001, -.2, -.01, -.005, -.4)
varNames = c("GDP", "GDP Deflator", "3-month Tbill", "Investment", "S&P 500", "Spread")

pdf(file="./figure1_R.pdf")
par(mfrow=c(3,2),mar=c(2.5,2.5,2,0.5))
for(nn in 1:n){
  ylim1=range(irf[,nn,shock,], na.rm=TRUE)
  plot.ts(irf[4,nn,shock,], col="black", xlab="", ylab="", ylim=c(ylim_l[nn], ylim_u[nn]), main=varNames[nn], axes=FALSE)
  grid(col="grey60", lty=1, lwd=0.5)
  polygon(c(1:horizon,rev(1:horizon)), c(irf[1,nn,shock,],rev(irf[7,nn,shock,])),
          col = "grey80", border=NA)
  polygon(c(1:horizon,rev(1:horizon)), c(irf[2,nn,shock,],rev(irf[6,nn,shock,])),
          col = "grey50", border=NA)
  polygon(c(1:horizon,rev(1:horizon)), c(irf[3,nn,shock,],rev(irf[5,nn,shock,])),
          col = "grey30", border=NA)
  lines(irf[4,nn,shock,], col="black", lwd=2, lty=1)
  abline(h=0, col="black", lty=1)
  axis(1, at=seq(1,horizon,by=5), labels=seq(0,horizon,by=5),lwd=2,font=2)
  axis(2, font=2, lwd=2)
  box(lwd=2)
}
dev.off()












