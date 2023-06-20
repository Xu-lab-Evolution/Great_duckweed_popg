popC <- read.table("countsAME.frq.count",header=FALSE,skip=1)
popB <- read.table("countsASIA.frq.count",header=FALSE,skip=1)
popA <- read.table("countsIND.frq.count",header=FALSE,skip=1)

for (i in 1:20) {
#Select here the chromosome
chr <- i
print(paste("Processing chromosome",chr))
genetic_positions <- read.table(paste0("genetic_positions_chr",chr,".txt"),header=TRUE)
popC_chr <- subset(popC,popC$V1==chr)
popB_chr <- subset(popB,popB$V1==chr)
popA_chr <- subset(popA,popA$V1==chr)

derived_popC_chr <- apply(popC_chr[,5:6],1,min)
derived_popB_chr <- apply(popB_chr[,5:6],1,min)
derived_popA_chr <- apply(popA_chr[,5:6],1,min)

input <- data.frame(popC_chr$V1,genetic_positions,derived_popA_chr,popA_chr$V4,derived_popB_chr,popB_chr$V4,derived_popC_chr,popC_chr$V4)
#colnames(input) <- c("chr","physpos","genpos","mEUR","nEUR","mASIA","nASIA","mAME","nAME")
colnames(input) <- c("chr","physpos","genpos","mIND","nIND","mASIA","nASIA","mAME","nAME")

input <- input[input$mAME!=0,] #Get rid of any remaining monomorphic sites in the outgroup AME.
write.table(input,file=paste0("input_Chr",chr,".txt"),row.names=FALSE,quote=FALSE,sep="\t")
}

