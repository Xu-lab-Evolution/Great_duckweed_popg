library(abc)
library(spatstat)

simulations <- read.table("sim_tot_m3.txt.transformed",header=TRUE)
P.mod01 <- simulations[,1:14]
S.mod01 <- simulations[,27:56] #Original: 56
S.obs <- read.table("obs_intergenic_mean_var.txt.transformed",header=TRUE)
S.obs <- S.obs[,9:38] #Original: 38

rejection <- abc(S.obs,P.mod01,S.mod01,tol=0.1,method="rejection")
regression <- abc(S.obs,P.mod01,S.mod01,tol=0.1,method="loclinear")

par(mfrow=c(3,5))
for (i in 1:length(names(P.mod01))) {
    a <- density(regression$adj.values[,i], weights=regression$weights/sum(regression$weights))
    maximum <- a$x[which(a$y==max(a$y))]
    plot(a,main=names(P.mod01)[i],xlab=paste("Mode:",round(maximum,3)))
    #print confidence intervals
    qs <- quantile.density(a,probs=c(0.05,0.95)) #Function from library spatstat
    print(paste(names(P.mod01)[i],qs[1],qs[2]))
}

cv <- cv4abc(P.mod01,S.mod01,tols=0.05,method="loclinear",nval=100)

