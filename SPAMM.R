SPAMM <- function(uhats,threshold=FALSE,tval=10^(-12)) {
  Pis <- list()
  Js <- list()
  for (k in c(1:length(uhats))) {
    R <- uhats[[k]]
    J <- NULL
    j = 1
    r <- dim(R)[2]
    while (j <= r) {
      jstar <- which.max(apply(R,1,function(x){
        sum(x^2)
      }))
      vj <- R[jstar,]
      R <- R %*% ( diag(1,r) - vj %*% t(vj) / sum(vj^2))
      J <- union(J,jstar)
      j <- j+1
    }
    
    Pihat <- uhats[[k]] %*% solve(uhats[[k]][J,])
    if(threshold) {
      Pihat <- ifelse(Pihat < tval,0,Pihat)
      vec <- Pihat %*% rep(1,dim(Pihat)[2])
      Pihat <- diag(as.vector(1/vec),length(vec)) %*% Pihat
    } else {
      vec <- abs(Pihat) %*% rep(1,dim(Pihat)[2])
      Pihat <- diag(as.vector(1/vec),length(vec)) %*% Pihat
    }
    Pis[[k]] <- Pihat
    Js[[k]] <- J
  }
  return(list(Pis,Js))
}