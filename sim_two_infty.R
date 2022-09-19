source("misc.R")
source("SPAMM.R")
source("tensor_two_infty.R")
numcores<- detectCores()
cl <- makeCluster(4,type = "FORK")
registerDoParallel(cl,numcores)

print("got here!")
rs <- rep(3,3)
set.seed(9192022)
C <- array(rnorm(prod(rs)), dim=rs) # Core tensor
C <- as.tensor(C)
# adjust minimal separation.
#delta.min = 100000
delta = 10
Sk = k_unfold(C, 1)@data
delta.min = min(svd(Sk)$d)
for (i in 1:3){
  Sk = k_unfold(C, i)@data
  for (k1 in 1:(rs[i]-1)){
    for (k2 in (k1+1):rs[i]){
      delta.min = min(delta.min, min(svd(Sk)$d))
    }
  }
}
C = C * delta / delta.min

r <- 3

#C <- array( ,dim=rep(r,3))
#C <- as.tensor(C)
#svd(k_unfold(C,1)@data )$d
ps <- seq(100,500,50)#,50)#500,50)
ntrials <- 10
sigmas <- seq(1,100,5)#100,10)


finalres_uniform <- foreach(p = ps) %dopar% {
  print(paste0("p=",p))
  toreturn2 <-  foreach(n = c(1:ntrials),.combine=cbind) %dopar% {
    print(paste0("n=",n))
    toreturn <- foreach(j = c(1:length(sigmas)),.combine=rbind) %dopar% {
      #print(j)
      return(runonesim(p,r,C,sigmamax = sigmas[j]))
     }
   return(toreturn)
  }
  return(toreturn2)
}
  
#rowMeans(finalres[[1]])/sigmas
#rowMeans(finalres[[2]])/sigmas
#rowMeans(finalres[[3]])/sigmas


save(finalres_uniform,file = "sim1_9-19.Rdata")
print("first sim done!")

finalres_het1 <- foreach(p = ps) %dopar% {
  print(paste0("p=",p))
  toreturn2 <-  foreach(n = c(1:ntrials),.combine=cbind) %dopar% {
    print(paste0("n=",n))
    toreturn <- foreach(j = c(1:length(sigmas)),.combine=rbind) %dopar% {
      #print(j)
      return(runonesim(p,r,C,sigmamax = sigmas[j],het_amt = .75))
    }
    return(toreturn)
  }
  return(toreturn2)
}

save(finalres_het1,file = "sim2_9-19.Rdata")
print("second sim done!")

finalres_het2 <- foreach(p = ps) %dopar% {
  print(paste0("p=",p))
  toreturn2 <-  foreach(n = c(1:ntrials),.combine=cbind) %dopar% {
    print(paste0("n=",n))
    toreturn <- foreach(j = c(1:length(sigmas)),.combine=rbind) %dopar% {
      #print(j)
      return(runonesim(p,r,C,sigmamax = sigmas[j],het_amt = .5))
    }
    return(toreturn)
  }
  return(toreturn2)
}

save(finalres_het2,file = "sim3_9-19.Rdata")
print("third sim done!")

finalres_het3 <- foreach(p = ps) %dopar% {
  print(paste0("p=",p))
  toreturn2 <-  foreach(n = c(1:ntrials),.combine=cbind) %dopar% {
    print(paste0("n=",n))
    toreturn <- foreach(j = c(1:length(sigmas)),.combine=rbind) %dopar% {
      #print(j)
      return(runonesim(p,r,C,sigmamax = sigmas[j],het_amt = .25))
    }
    return(toreturn)
  }
  return(toreturn2)
}

save(finalres_het3,file = "sim4_9-19.Rdata")
print("finished!")
