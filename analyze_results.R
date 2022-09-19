hetamt1 <- 1
load("sim1_9-15.Rdata")
ps <- seq(100,500,50)
sigmas <- seq(1,100,5)
result1 <- matrix(0,length(sigmas),length(ps))
rownames(result1) <- sigmas
colnames(result1) <- as.character(ps)
for (i in c(1:length(ps)) ) {
  result1[,i] <- as.numeric(rowMeans(finalres_uniform[[i]]))
}
result1


load("sim2_9-15.Rdata")
result2 <- matrix(0,length(sigmas),length(ps))
rownames(result2) <- sigmas
colnames(result2) <- as.character(ps)
for (i in c(1:length(ps)) ) {
  result2[,i] <- as.numeric(rowMeans(finalres_het1[[i]]))
}
result2


load("sim3_9-15.Rdata")
result3 <- matrix(0,length(sigmas),length(ps))
rownames(result3) <- sigmas
colnames(result3) <- as.character(ps)
for (i in c(1:length(ps)) ) {
  result3[,i] <- as.numeric(rowMeans(finalres_het2[[i]]))
}
result3

load("sim4_9-15.Rdata")
result4 <- matrix(0,length(sigmas),length(ps))
rownames(result4) <- sigmas
colnames(result4) <-as.character(ps)
for (i in c(1:length(ps)) ) {
  result4[,i] <- as.numeric(rowMeans(finalres_het2[[i]]))
}
result4

library(reshape)
resultattempt <- data.frame(rbind(result1,result2,result3,result4))
resultattempt$sigmaval <- as.numeric(rep(rownames(result1),4))
resultattempt$hetamt <- as.factor(c(rep(1,length(sigmas)),rep(.75,length(sigmas)),
                                    rep(.5,length(sigmas)),rep(.25,length(sigmas))))

results <- melt(resultattempt,id=c("hetamt","sigmaval"))
results$pfactor <- sub('X','',results$variable)
results$p <- as.numeric(sub('X','',results$variable))
#plot 1: a curve as sigma increases
library(ggplot2)
g <- ggplot(data=results[which(results$hetamt == "1"),],aes(x=sigmaval,y=value,color=as.factor(p)))
g1 <- g + geom_point() + geom_line(linetype="dashed") + xlab("sigma_max") +
  ylab("2 -> infty error") +  guides(color=guide_legend(title="p")) +
  ggtitle("2 -> infty error, averaged over 10 runs")

results$relativeerror <- results$value/results$sigmaval
library(dplyr)

results <-tibble(results)
results_new <- results
results_new$value <- NULL
results_new$sigmaval <- NULL
results_new <- results_new %>% group_by(pfactor,hetamt) %>%
  mutate(meanerror = mean(relativeerror))


g <- ggplot(data=results_new,aes(x=p,y=meanerror,color=hetamt))
g2 <- g + geom_point() + geom_line(linetype="dashed") +
  #ggtitle("2 -> infty error, normalized by sigmamax and then averaged") +
  ylab("Mean relative 2-> infty error")

png("sigma_error.png",units='in',width=5,height=5,res=300)
g1
dev.off()

png("p_relative_error.png",units='in',width=5,height=5,res=300)
g2
dev.off()





