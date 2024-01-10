args <- commandArgs(trailingOnly=TRUE)
data <- read.table(args[1],header=FALSE)
means <- apply(data,2,mean)
vars <- apply(data,2,var)
mv <- rbind(means,vars)

res <- as.vector(mv)

cat(res,sep="\t")

