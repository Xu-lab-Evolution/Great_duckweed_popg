input_folder <- "input_IND_mogon"
output_folder <- "input_drift_IND"

for (i in 1:20) {
chr <- i
data <- read.table(paste0(input_folder,"/input_Chr",chr,".txt"),header=TRUE)
input <- data[,4:9]
num <- rep(1,length(input[,1]))
input <- cbind(input,num)
write.table(input,file=paste0(output_folder,"/input_drift_Chr",chr,".txt"),row.names=FALSE,quote=FALSE,sep="\t")
}

