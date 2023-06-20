#FUNCTIONS
plot_clr <- function(data,num_column_clr,nameTag) {
	clr_branch <- c()
	colores <- c()
	col1 <- "gray53"
	col2 <- "gray"
	pos <- numeric(numChrom)
	end_chr_indices <- numeric(numChrom)
	for (i in 1:numChrom) {
		data_chr <- subset(data,data$Chr==i)
		branch <- data_chr[,num_column_clr]
		clr_branch <- c(clr_branch,branch)
		len <- length(branch)
		if (i==1) {
			end_chr_indices[i] <- len
			pos[i] <- round(len/2)
		} else {
			end_chr_indices[i] <- end_chr_indices[i-1] + len
			pos[i] <- end_chr_indices[i-1]+round(len/2)
		}
		if (i%%2==1) {
			colores <- c(colores,rep(col1,len))
		} else {
			colores <- c(colores,rep(col2,len))
		}
	}
	thres_1percent <- quantile(clr_branch,probs=0.99)
	colores[which(clr_branch>thres_1percent & colores==col1)] <- "firebrick1"
	colores[which(clr_branch>thres_1percent & colores==col2)] <- "firebrick1" #"chocolate1"
	plot(clr_branch,main=nameTag,col=colores,pch=20,xaxt="n",xlab="Chromosome",ylab="CLR")
	axis(1,at=pos,labels=1:numChrom)
	abline(h=thres_1percent,lty=2,lwd=0.5)
	#write.table(data[which(clr_branch>thres_1percent),1:2],file="ancestral_top1_regions.csv",quote=FALSE,sep=",",row.names=FALSE,col.names=FALSE) #Get list of top 1% regions for the Venn diagram.
}

#MAIN

#PLOTS
#pdf(paste0("figure_top_europe.pdf"),width=12,height=3)
#pdf(paste0("figure_top_asia.pdf"),width=12,height=3)
pdf(paste0("figure_top_india.pdf"),width=12,height=3)

#Plot 3P-CLR scan
#output_folder <- "/home/pablo/sciebo/Muenster/Duckweed/GenomeScan/3P-CLR/outputs_EUR"
output_folder <- "/home/pablo/sciebo/Muenster/Duckweed/GenomeScan/3P-CLR/outputs_IND"
data <- read.table(paste0(output_folder,"/output_Chr_all.txt"),header=TRUE)
standardized_clr_scores_pop <- numeric(dim(data)[1]) #suffix pop, could be europe or india, etc.
data <- cbind(data,standardized_clr_scores_pop)
numChrom <- 20
for (i in 1:numChrom) { #For each of the 20 chromosomes
    #data$standardized_clr_scores_pop[data$Chr==i] <- (data$X3PCLR.A[data$Chr==i]-mean(data$X3PCLR.A[data$Chr==i]))/sd(data$X3PCLR.A[data$Chr==i]) #For Europe
    #data$standardized_clr_scores_pop[data$Chr==i] <- (data$X3PCLR.B[data$Chr==i]-mean(data$X3PCLR.B[data$Chr==i]))/sd(data$X3PCLR.B[data$Chr==i]) #For Asia
    data$standardized_clr_scores_pop[data$Chr==i] <- (data$X3PCLR.A[data$Chr==i]-mean(data$X3PCLR.A[data$Chr==i]))/sd(data$X3PCLR.A[data$Chr==i]) #For India (check output_folder is IND)
}
#plot_clr(data,18,"European branch") #Standardized Europe
#plot_clr(data,18,"SE-Asian branch") #Standardized Asia
plot_clr(data,18,"Indian branch") #Standardized India


dev.off()
