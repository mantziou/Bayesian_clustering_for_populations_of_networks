rbetat <- function(n, range, alpha, beta) {
  
  # range is a vector of two values
  
  F.a <- pbeta(min(range), alpha, beta)
  F.b <- pbeta(max(range), alpha, beta)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qbeta(u, alpha, beta)
  
}
