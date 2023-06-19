#args <- commandArgs(trailingOnly=TRUE)
args <- "all_obs_stats_header.csv"
data <- read.table(args[1],header=TRUE)
means <- apply(data,2,mean)
vars <- apply(data,2,var)
mv <- rbind(means,vars)

header <- c()
for (i in 1:length(names(data))) {
    header <- c(header,paste0(names(data)[i],"_mean"))
    header <- c(header,paste0(names(data)[i],"_var"))
}

res <- as.vector(mv)
names(res) <- header

write.table(t(res),file="obs_intergenic_mean_var.txt",quote=FALSE,row.names=FALSE,col.names=TRUE)
