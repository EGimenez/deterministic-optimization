source('new_way.r')

optimize_f_ij <- function(S, K, i, j){
  f <- poli_fun(S, K, i, j)
  
  if (i == j){
    x0 <- (K[i, j])^.5
  }else{
    x0 <- (-K[i, j])^.5
  }
  
  
  result <- optim(x0, f, control=list(fnscale=-1, reltol=1e-8))
  
  if(i == j){
    K[i,j] <- (result$par)^2
  }else{
    K[i, j] <- -(result$par)^2
    K[j, i] <- -(result$par)^2
  }
  K
}

exercise_3 <- function(S, pre = 1e-1, max_iter = 100){
  
  n <- dim(S)[1]
  
  eps <- 0.1
  
  while (TRUE){
    K <- diag(n)
    K <- K - eps*matrix(1, nrow = n, ncol = n) 
    
    if (det(K)>0){
      break
      }else{
        eps<-eps/2
      } 
  }

  
  last_value <- Inf
  
  # We make sure that at least each component is maximized once
  for (i in 1:n){
    for (j in i:n){
      aux <- optimize_f_ij(S, K, i, j)
      K <- aux
    }
  }
  
  # Here we could be fancyer
  for (mi in 1:max_iter){
    ij = find_ij(S, K)
    
    i <- ij[1]
    j <- ij[2]
    
    aux <- optimize_f_ij(S, K, i, j)
    K <- aux
  }
  round(K, 2)
}