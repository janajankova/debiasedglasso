## De-biased glasso function
## 
## This function allows you to calculate:
## (i) confidence intervals for individual entries Theta[i,j] of a precision matrix Theta
## (ii) p-values (adjusted for multiple testing) for testing a hypothesis of type H_{i,j}: Theta[i,j] = a[i,j]
## @param X: 

debiased.glasso <- function(X, pval = FALSE, scale.X = TRUE, p = ncol(X), n = nrow(X), lambda = sqrt(log(p)/n), 
                            theta.H0 = matrix(0,ncol(X), ncol(X)), alpha = 0.05){
  
  if(scale.X == TRUE)
    X <- scale(X)  
  
  S <- var(X)
  
  theta.glasso <- method.glasso(S, p, n)
  
  theta.desp <- 2 * theta.glasso - theta.glasso %*% S %*% theta.glasso
  
  # estimate of standard deviation
  sigma.hat <- sqrt(theta.glasso^2 +  diag(theta.glasso) %*% t(diag(theta.glasso)))
  
  
  if(pval)
    p.values(theta.desp, theta.H0, sigma.hat, n, alpha)
  
  else
    ci(theta.desp, sigma.hat, n, alpha)
  
}

ci <- function(theta.desp, sigma.hat, n, alpha){
  
  # confidence intervals
  u <- qnorm(1-alpha/2)
  CI.lower.end <- theta.desp - u * sigma.hat / sqrt(n)
  CI.upper.end <- theta.desp + u * sigma.hat / sqrt(n)
  
  list(CI.lower.end,CI.upper.end)
  
}


p.values <- function(theta.desp, theta.H0, sigma.hat, n, alpha = 0.05){
  pval <- 1 - pnorm(abs(theta.desp-theta.H0)/sigma.hat*sqrt(n), mean = 0, sd = 1)
  
  pval.adj <- matrix(p.adjust(pval, method = "BH", n = length(pval)), nrow = nrow(theta.desp))
  selected <- which(pval.adj < alpha, arr.ind=TRUE)
  
  list(pval.adj, sum(pval.adj < alpha), selected)
  
}


method.glasso <- function(S, p, n){
  
  glasso(s = S, rho = 2*sqrt(log(p)/n), penalize.diagonal = FALSE)$wi
  
}

