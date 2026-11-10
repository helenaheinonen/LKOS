library(EbayesThresh)
library(copula)


# eg. phi(x,2) = x(1-log(x))
phi <- function(x, k){
  sapply(x, function(y){
    y*sum((-log(y))^(0:(k-1)) / factorial(0:(k-1)))
  })
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
  sapply(x, function(y){
    y*sum( choose(n,0:(k-1)) * (y^(-1/n)-1)^(0:(k-1)))  
  })
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
    p <- p * W[,(j-1)]^(n/(n-j+1)) # update the product with each W_i; j-1 because the Ws are off set by 1 
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
  
  e[2:k] <- 
    inv_phi(u[2:k], 2:k) /
    inv_phi(u[1:(k-1)] , 1:(k-1))
  e
}

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
  
  j<-2:k
  e[j] <- (inv_phiFin(u[j], j, n) / inv_phiFin(u[j-1] , j-1, n))^( (n-j+1)/n ) 
  e
}

# takes in a copula sample of the largest 2 order statistics
# and computes the rosenblatt transform using the copula for a finite block size n

RtransFin <- function(data, n){
  return(t(apply(data, 1, vecRtransFin, n=n)))  
}



#setting up bootstrap process 


# joint empirical CDF evaluated at u
# in the testing process, input data should be pobs 
D <- function(u, data){
  if (length(u) != ncol(data)) {
    stop("Length of u must match number of columns.")
  }
  return(mean(apply(data <= matrix(u, nrow(data), length(u), byrow = TRUE), 1, all)))
}


# input a copula sample 
# computes test stat (distance from the indep copula)
# this is test stat S_n^c from Genest copula paper (cramer-von mises)
testStat<-function(data){
  sum(
    apply(data, 1, function(row) { #apply this function to each row of data
      (D(row, data) - prod(row))^2
    })
    )
}


# reports p-value for rosenblatt goodness of fit test
# input a copula sample of largest and second largest (not Rosenblatt transformed)
# this should be pobs of the data
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
  
  testStat.bootstrapped <- replicate(bootstraps, {
      if (finite){
          # simulate from null copula 
          hyp_cop <- rcopFin(nrow(sample), k = ncol(sample), blockSize = n)
          hyp_cop <- pobs(hyp_cop) # ranks
          # no need for estimation bc no theta
          transformed <- RtransFin(hyp_cop, n) #rosenblatt
      }
      else{ # testing against limit copula 
        hyp_cop <- rcop(nrow(sample), k = ncol(sample))
        hyp_cop <- pobs(hyp_cop)
        # no theta
        transformed <- RtransLim(hyp_cop) #rosenblatt 
      }
      testStat(transformed)
  })
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
