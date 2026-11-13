library(EbayesThresh)
library(copula)


# eg. phi(x,2) = x(1-log(x))
phi <- function(x, k){
  mapply(function(y, j){
    y*sum((-log(y))^(0:(j-1)) / factorial(0:(j-1)))
  }, x, k)
  
  # sapply(x, function(y){
  #   y*sum((-log(y))^(0:(k-1)) / factorial(0:(k-1)))
  # })
}

# note that phi_k(x) = P(gamma(k, 1) > -ln(y))

# inverse of phi_k is related to the gamma function 
# inv_phi = psi
inv_phi <- function(x, k){
  t <- qgamma(1 - x, shape = k, rate = 1)
  return(exp(-t))
}

#phi_{k,n}(x) for finite block size n 
phiFin <- function(x, k, n){ ## note here n is block size not sample size 
  
  mapply(function(y, j, b){
    y*sum( choose(b, 0:(j-1)) * (y^(-1/b) -1)^(0:(j-1)))
  }, x, k, n)
  # sapply(x, function(y){
  #   y*sum( choose(n,0:(k-1)) * (y^(-1/n)-1)^(0:(k-1)))  
  # })
}

#inverse of phi_{k,n}(x) for finite block size n
inv_phiFin <- function(x, k, n){
  t <- 1 - qbeta(1-x, k, n-k+1)
  return(t^n)
}

# generate a sample of size n from a copula for the largest k order statistics
rcop <- function(n, k){
  U1 <- runif(n)
  
  W <- matrix(runif(n * (k - 1)), nrow = n) # so W[,1] is column of W1s, etc. 
  
  U <- matrix(NA, nrow=n, ncol=k) # will store the copula sample 
  U[,1] <- U1
  
  p <- U1 # running product 
  for(j in 2:k){  # index j tracks which copula column we are updating
    p <- p * W[,(j-1)] # update the product with each W_i; j-1 because the Ws are off set by 1 
    U[,j] <- phi(p, j) 
  }
  return(U)
}

# simulate a sample of size n from the copula of the LKOS for a finite block size
rcopFin <- function(n, k, blockSize){
  U1 <- runif(n)
  
  W <- matrix(runif(n * (k-1)), nrow = n)
  
  U <- matrix(NA, nrow=n, ncol=k) # will store the copula sample 
  U[,1] <- U1
  
  p <- U1 # running product 
  for(j in 2:k){  # index j tracks which copula column we are updating
    p <- p * W[,(j-1)]^(blockSize/(blockSize-j+1)) # update the product with each W_i; j-1 because the Ws are off set by 1 
    U[,j] <- phiFin(p, j, blockSize) 
  }
  
  return(U)
}


# takes in a vector and calculates the rosenblatt transform using the limit copula 
# used as helper function
vecRtransLim <- function(u){
  k <- length(u)
  e <- numeric(k)
  
  e[1] <- u[1] # first is unchanged

  if (k==1) return(e)
  
  # compute all the inverses once (should speed things up)
  v <- inv_phi(u, 1:k) 
  bdry <- phi(v[1:(k-1)], 2:k)
  
  e[2:k] <- ifelse(u[2:k] >= bdry,  1, v[2:k] / v[1:(k-1)])
  e
}


# # for now: updated rosenblatt transform for just the largest 2 
# vecRtransLim.TT <- function(u){
#   k <- length(u)
#   e <- numeric(k)
#   
#   e[1] <- u[1]
#   e[2] <- if (u[2] >= u[1]*(1-log(u[1]))) 1 else inv_phi(u[2], 2) / inv_phi(u[1] , 1)
#   e
# }

# for a whole copula sample in matrix
RtransLim <- function(data){
  return(t(apply(data, 1, vecRtransLim)))
}


# adjusting Rosenblatt test to account for finite block size

F2 <- function(y, n){
  return(n*y^(n-1) - (n-1)*y^(n))
}
F2_inverse <- function(x, n){
  f <- function(y) F2(y, n) - x
  uniroot(f, interval = c(0, 1), tol = 1e-10)$root
}



# helper function
vecRtransFin <- function(u, n){
  k <- length(u) 
  e <- numeric(k)
  
  e[1] <- u[1] # first is unchanged
  if (k==1) return(e)
  
  # compute all the inverses once
  v <- inv_phiFin(u, 1:k, n)
  bdry <- phiFin(v[1:(k-1)], 2:k, n)
  
  e[2:k] <- ifelse(u[2:k] >= bdry, 1, (v[2:k] / v[1:(k-1)])^((n-(2:k) + 1)/n) )

  e
}

# # adjusted temp function for just the largest 2 
# vecRtransFin.TT <- function(u, n){
#   k <- length(u) 
#   e <- numeric(k)
#   
#   e[1] <- u[1] # first is unchanged
#   
#   e[2] <- if (u[2] >= phiFin(u[1], k, n)) 1 else (inv_phiFin(u[2], 2, n) / inv_phiFin(u[1] , 1, n))^( (n-j+1)/n ) 
#   e
# }

# takes in a copula sample of the largest 2 order statistics
# and computes the rosenblatt transform using the copula for a finite block size n

RtransFin <- function(data, n){
  return(t(apply(data, 1, vecRtransFin, n=n)))  
}



#setting up bootstrap process 


# input pseudo-observations 
# 
# computes test stat (distance from the indep copula)
# this is test stat S_n^c from Genest copula paper (cramer-von mises)
testStat<-function(data){
  # sum(
  #   apply(data, 1, function(row) { #apply this function to each row of data
  #     (D(row, data) - prod(row))^2
  #   })
  #   )
  
  #C.n should take in pseudo-observations
  #i am feeding it rosenblatt transformed pseudo-observations 
  
  sum(
    apply(data, 1, function(row){ # apply this function to each row of data
      (C.n(row, data) - prod(row))^2
    })
  )
}


# reports p-value for rosenblatt goodness of fit test
# input a copula sample (pobs) of largest and second largest (not Rosenblatt transformed)
# finite = TRUE means you are testing the finite block size copula, n = block size
# finite = FALSE testing the limit copula 
# bootstraps = number of bootstraps 
rosenblattTest <- function(sample, finite = FALSE, n = NULL, bootstraps = 1000){
  
  if (finite && is.null(n)) {
      stop("Please provide n when finite = TRUE.")
  }
  
  # Test statistic for the observed sample
  if (finite) {
    # we do the testStat of the transformed sample since after transforming it should be indep
    # computes cramer-von mises difference to the indep cop 
    testStat.sample <- testStat(RtransFin(sample, n))
  } else {
    testStat.sample <- testStat(RtransLim(sample))
  }

  # bootstrapped test stat
  
  # testStat.bootstrapped <- replicate(bootstraps, {
  #     if (finite){
  #         # simulate from null copula 
  #         hyp_cop <- rcopFin(nrow(sample), k = ncol(sample), blockSize = n)
  #         hyp_cop <- pobs(hyp_cop) # ranks
  #         # no need for estimation bc no theta
  #         transformed <- RtransFin(hyp_cop, n) #rosenblatt
  #     }
  #     else{ # testing against limit copula 
  #       hyp_cop <- rcop(nrow(sample), k = ncol(sample))
  #       hyp_cop <- pobs(hyp_cop)
  #       # no theta
  #       transformed <- RtransLim(hyp_cop) #rosenblatt 
  #     }
  #     testStat(transformed)
  # })
  L <- as.list(rep(nrow(sample),bootstraps)) #list to contain the bootstrap replicates
  if(finite){
    # generate boot replicates from the null
    L <- lapply(L,FUN=function(l){rcopFin(l, k = ncol(sample), blockSize = n)})
    # take their ranks, no need for estimation because no theta
    L <- lapply(L, FUN=function(l){pobs(l)})
    # transform all the nulls 
    M <- lapply(L,FUN=function(sample){testStat(RtransFin(sample, n=n))})
  }
  else{
    L <- lapply(L,FUN=function(l){rcop(l, k = ncol(sample))})
    L <- lapply(L, FUN=function(l){pobs(l)})
    M <- lapply(L,FUN=function(sample){testStat(RtransLim(sample))})
  }
  # M is a list that contains all the test stats of the bootstrap samples 
  testStat.bootstrapped <- unlist(M)
  mean(testStat.bootstrapped > testStat.sample)
}
  
  
# quantile functions for largest and second largest margins 


gevq1 <- function(a, z){
  #mu, sigma, xi
  if (a[3] != 0){
    a[1] + a[2]*( (-log(z))^(-a[3]) - 1)/a[3]
  }
  else{
    -a[2]*ln(-log(z)) + a[1]
  }
}

gevq2 <- function(a, z){
  if (a[3] != 0){
    a[1] + a[2]*( (-log(inv_phi(z,2)))^(-a[3]) - 1)/a[3]
  }
  else{
    -a[2]*ln(-log(inv_phi(z,2))) + a[1]
  }
}

# marginal CDF of second largest 
gevf2 <- function(a, z){
  gevf(a, z) * (1 - log(gevf(a,z)))
}
