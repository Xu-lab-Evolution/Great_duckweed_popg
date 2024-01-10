library(abc)

numSim <- 20369

sims01 <- read.table("sim_tot_m1.txt",header=TRUE,nrows=numSim,fill=TRUE)
sims02 <- read.table("sim_tot_m2.txt",header=TRUE,nrows=numSim,fill=TRUE)
sims03 <- read.table("sim_tot_m3.txt",header=TRUE,nrows=numSim,fill=TRUE)
obs <- read.table("obs_intergenic_mean_var.txt",header=TRUE)

stats2check <- c("S_pop_AME_mean","S_pop_ASI_mean","S_pop_EUR_mean","S_pop_IND_mean",
                 "TajD_pop_AME_mean","TajD_pop_ASI_mean","TajD_pop_EUR_mean","TajD_pop_IND_mean",
                 "ZnS_pop_AME_mean","ZnS_pop_ASI_mean","ZnS_pop_EUR_mean","ZnS_pop_IND_mean",
                 "daNei_AME_ASI_mean","daNei_AME_EUR_mean","daNei_AME_IND_mean",
                 "daNei_ASI_EUR_mean","daNei_ASI_IND_mean","daNei_EUR_IND_mean",
                 "S_pop_AME_var","S_pop_ASI_var","S_pop_EUR_var","S_pop_IND_var",
                 "TajD_pop_AME_var","TajD_pop_ASI_var","TajD_pop_EUR_var","TajD_pop_IND_var",
                 "ZnS_pop_AME_var","ZnS_pop_ASI_var","ZnS_pop_EUR_var","ZnS_pop_IND_var",
                 "daNei_AME_ASI_var","daNei_AME_EUR_var","daNei_AME_IND_var",
                 "daNei_ASI_EUR_var","daNei_ASI_IND_var","daNei_EUR_IND_var")

simRed01 <- sims01[,stats2check]
simRed02 <- sims02[,stats2check]
simRed03 <- sims03[,stats2check]
obsRed <- obs[,stats2check]

simRed01 <- simRed01[complete.cases(simRed01),]
simRed02 <- simRed02[complete.cases(simRed02),]
simRed03 <- simRed03[complete.cases(simRed02),]

allSimulations <- rbind(simRed01,simRed02,simRed03)
indices <- c(rep("M01",dim(simRed01)[1]),rep("M02",dim(simRed02)[1]),rep("M03",dim(simRed03)[1]))

tolerance <- 0.05

cv.model.choice <- cv4postpr(indices,allSimulations,nval=100,tols=tolerance,method="mnlogistic")
summary(cv.model.choice)

model.choice <- postpr(obsRed,index=indices,sumstat=allSimulations,tol=tolerance,method="mnlogistic")
summary(model.choice)

#####################################

