
rm(list=ls())

library(Matrix)
library(Rcpp)
library(MHadaptive)
library(coda)
library(MASS)
library(circular)

sourceCpp(file = "fillmats.cpp")
sourceCpp(file = "fillvec.cpp")

source("buildmap.r")
source("genmap.r")
source("p2psims.r")
source("p2pmcmc.r")
source("dep2p.r")

# jmm 01-17-2017


# likelihood function to be called in the MCMC ---------------------------------
like <- function(newp){
  lk=NA
  a1= (newp[1])
  a2= (newp[2])
  a3= (newp[3])
  b1 = (newp[4])
  b2 = (newp[5])
  b3 = (newp[6])
  shape = newp[7]
  
  if(a3 > 0 & a2 > 0 & a1 > 0 & shape > 0){
    
    PT <- dep2p(n=n, distance=distance, radius=radius,a1=a1,a2=a2,a3=a3,b1=b1,b2=b2, b3=b3)
    P = PT[[1]]
    if(!any(!is.finite(P))){
      G = PT[[2]]
      if(!any(!is.finite(G)) & !any(G < 0)){
        uid = unique(id)
        lk = 0 #nlliks = numeric(length(uid))
        for(i in 1:length(uid)){
          tmp = which(id==uid[i])
          if(length(tmp) > 1 ){
            z <- zs[tmp]
            tt <- tts[tmp]
            nsteps = length(z)
            lik = numeric(nsteps-1)
            for(k in 1:(nsteps-1)){
              
              lik[k] = P[z[k],z[k+1]]
            }
            lk = lk + sum(log(lik))
          }
          if(length(tmp)>2){
            for(k in 1:(nsteps-2)){
              lk = lk + dgamma(tt[k+1],scale=G[z[k],z[k+1]]/shape,shape=shape,log=TRUE)
            }
          }
        }
      }
    }
  }
  return(lk)
}
# ------------------------------------------------------------------------------
npars <- 7
npas <- c(10, 30, 50)
kaps <- c(0, 5, 50)
tol = 1 
minr = 0.0002
maxr = 0.25

nreps <- 10 

M  <- matrix(NA, nreps * length(npas) * length(kaps), 2 + npars)
Lo <- matrix(NA, nreps * length(npas) * length(kaps), 2 + npars)
Hi <- matrix(NA, nreps * length(npas) * length(kaps), 2 + npars)
Rh <- matrix(NA, nreps * length(npas) * length(kaps), 2 + npars)
Ne <- matrix(NA, nreps * length(npas) * length(kaps), 2 + npars)
fila <- 0

for(j in 1:nreps){
  for(jj in 1:length(npas)){
    for(i in 1:length(kaps)){
      fila <- fila + 1
      npatches <- npas[jj]
      kappa <- kaps[i]
      
      map <- genmap(npatches=npatches, tol=tol, minr=minr, maxr=maxr, minx = 0, maxx = 10)
      
      #map <- buildmap(npatches=npatches, tol=tol, minr=minr, maxr=maxr,  sc = 3)
      
      tmp <- p2psims(map = map, nind = 30, kappa = kappa, m = 1e-05)
      
      zs  <- tmp$zs
      id  <- tmp$id
      tts <- tmp$tt/10
      
      x <- map$x
      y <- map$y
      radius <- map$r  
      xv <- t(rbind(x,y))
      distance <- as.matrix(dist(xv, upper = TRUE, diag = TRUE))
      n <- npatches
      
      #-------------------------------------------------------------------------------
      
      count <- 0
      mrhat <- 10
      while(mrhat > 1.4 & count < 10){ 
        res <- p2pmcmc(like=like, npars=npars)
        count <- count + 1
        p.s <- res[[1]]
        pos <- res[[2]]
        
        X <- list(as.mcmc(p.s[, 1, ]), as.mcmc(p.s[, 2, ]), as.mcmc(p.s[, 3, ]))
        
        rhat <- gelman.diag(as.mcmc.list(X))
        ef <- effectiveSize(X)
        mrhat <- rhat$mpsrf
        
        if(mrhat < 1.4){
          X <- as.mcmc(rbind(p.s[, 1, ], p.s[, 2, ], p.s[, 3, ]))
          ci <- HPDinterval(X)
          M[fila, ]  <- c(kappa, npatches, colMeans(X))
          Lo[fila, ] <- c(kappa, npatches, ci[,1])
          Hi[fila, ] <- c(kappa, npatches, ci[,2])
          Rh[fila, ] <- c(kappa, npatches, rhat$psrf[,1])
          Ne[fila, ] <- c(kappa, npatches, ef)
          save(M,Lo,Hi,Rh,Ne, file = paste("res_tol", tol, ".RData", sep = "")) 
        }
      }
    }
  }
}