output_folder <- "outputs_IND"
data <- read.table(paste0(output_folder,"/output_Chr_all.txt"),header=TRUE)
standardized_clr_scores_anc <- numeric(dim(data)[1])
standardized_clr_scores_eur <- numeric(dim(data)[1]) #suffix eur, but could be europe or india.
standardized_clr_scores_asi <- numeric(dim(data)[1])
data <- cbind(data,standardized_clr_scores_anc,standardized_clr_scores_eur,standardized_clr_scores_asi)

for (i in 1:20) { #For each of the 20 chromosomes
	data$standardized_clr_scores_anc[data$Chr==i] <- (data$X3PCLR.Anc[data$Chr==i]-mean(data$X3PCLR.Anc[data$Chr==i]))/sd(data$X3PCLR.Anc[data$Chr==i])
	data$standardized_clr_scores_eur[data$Chr==i] <- (data$X3PCLR.A[data$Chr==i]-mean(data$X3PCLR.A[data$Chr==i]))/sd(data$X3PCLR.A[data$Chr==i])
	data$standardized_clr_scores_asi[data$Chr==i] <- (data$X3PCLR.B[data$Chr==i]-mean(data$X3PCLR.B[data$Chr==i]))/sd(data$X3PCLR.B[data$Chr==i])
}

#Ancestral population EUR (or IND) ASIA
top01_anc <- data[order(data$standardized_clr_scores_anc,decreasing=TRUE),]
ntop01_anc <- ceiling(length(top01_anc$Chr)/100) #Top 1%
top01_anc_final <- top01_anc[1:ntop01_anc,c(1,4,5,8,18)]

#EUR (or IND)
top01_eur <- data[order(data$standardized_clr_scores_eur,decreasing=TRUE),]
ntop01_eur <- ceiling(length(top01_eur$Chr)/100)
top01_eur_final <- top01_eur[1:ntop01_eur,c(1,4,5,10,19)]

#ASIA
top01_asi <- data[order(data$standardized_clr_scores_asi,decreasing=TRUE),]
ntop01_asi <- ceiling(length(top01_asi$Chr)/100)
top01_asi_final <- top01_asi[1:ntop01_asi,c(1,4,5,12,20)]

for (i in 1:length(top01_anc_final[,1])) {
	if (as.numeric(top01_anc_final[i,1])<10) {
		top01_anc_final[i,1] <- paste0("ChrS0",top01_anc_final[i,1])
	} else {
		top01_anc_final[i,1] <- paste0("ChrS",top01_anc_final[i,1])
	}
}
for (i in 1:length(top01_eur_final[,1])) {
	if (as.numeric(top01_eur_final[i,1])<10) {
		top01_eur_final[i,1] <- paste0("ChrS0",top01_eur_final[i,1])
	} else {
		top01_eur_final[i,1] <- paste0("ChrS",top01_eur_final[i,1])
	}
}
for (i in 1:length(top01_asi_final[,1])) {
	if (as.numeric(top01_asi_final[i,1])<10) {
		top01_asi_final[i,1] <- paste0("ChrS0",top01_asi_final[i,1])
	} else {
		top01_asi_final[i,1] <- paste0("ChrS",top01_asi_final[i,1])
	}
}

#Write bed files
write.table(top01_anc_final[,c(1,2,3)],file=paste0(output_folder,"/results_ancestral.bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(top01_eur_final[,c(1,2,3)],file=paste0(output_folder,"/results_india.bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(top01_asi_final[,c(1,2,3)],file=paste0(output_folder,"/results_asia.bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

#Write bed files plus standardized CLR scores
write.table(top01_anc_final[,c(1,2,3,5)],file=paste0(output_folder,"/results_ancestral_plusScores.bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(top01_eur_final[,c(1,2,3,5)],file=paste0(output_folder,"/results_india_plusScores.bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(top01_asi_final[,c(1,2,3,5)],file=paste0(output_folder,"/results_asia_plusScores.bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

