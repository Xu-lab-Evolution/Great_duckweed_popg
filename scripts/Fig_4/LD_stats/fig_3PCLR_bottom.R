#FUNCTIONS
plotGenes <- function(ann,chromstart,chromend,startY) {
   transFactor <- 1 #For x-axis to be given in Mega-bases change this value to 1e6
   yheights <- rep(c(0.1,1.2,2.3,3.4,4.5)+startY,ceiling(dim(ann)[1]/4))
   yheights <- yheights[1:dim(ann)[1]]
   xheights <- apply(ann[,2:3],1,mean)
   for (i in 1:length(ann$start)) {
      rect(ann$start[i]/transFactor,yheights[i],ann$stop[i]/transFactor,yheights[i]+1,border=NA,col="lightblue")
      text(xheights[i]/transFactor,yheights[i]+0.5,labels=gsub("SpGA2022_","",ann$name[i]),cex=1.25) #Change here for gene number size.
   }
}

#MAIN
library(Sushi)
options(scipen = 999)

args=commandArgs(trailingOnly=TRUE)

input_folder <- "/home/pablo/Documents/Duckweed/LD"

gtf <- read.delim("/home/pablo/Documents/Duckweed/annotation/SpGA2022.1.putative_function.domain_added_onlyGenes.gff",h=F,sep="\t",skip=1)

#Read candidate locus and window size
chr <- as.numeric(args[1]) #Chromosome
if (chr<10) {
  chrom <- paste0("ChrS0",chr)
} else {
  chrom <- paste0("ChrS",chr)
}

if (args[2]==0) {
   locus <- as.numeric(args[3]) #locus in chromosome
   window <- as.numeric(args[4]) #Window size
   chromstart <- floor(locus-window/2)
   chromend <- ceiling(locus+window/2)
} else if (args[2]==1) {
   locus <- as.numeric(args[3]) #locus in chromosome
   chromstart <- as.numeric(args[3])
   chromend <- as.numeric(args[4])
}

#Retrieve locus region from annotation file
ann <- subset(gtf,gtf$V1==chrom & gtf$V3=="gene" & gtf$V4>=chromstart & gtf$V5<=chromend)
ann <- ann[,c(1,4,5,9)]
a <- unlist(lapply(ann[,4],strsplit,";"))
a <- a[grep("^ID=",a)]
loci_names <- gsub("ID=","",a)
ann <- cbind(ann[,1:3],loci_names)
names(ann) <- c("chrom","start","stop","name")

#PLOTS
pdf(paste0("figure_bottom_chr9.pdf"),width=10,height=10) #Run as: Rscript fig_3PCLR_bottom.R 9 1 5950000 6050000
#pdf(paste0("figure_bottom_chr14.pdf"),width=10,height=10) #Run as: Rscript fig_3PCLR_bottom.R 14 1 770000 870000
mat=matrix(c(1,1,2,2,3,3,4),nrow=7,byrow=T)
layout(mat)
par(mar = c(3, 4, 1, 4) + 0.3)

#Plot statistics
pops <- c("ASIA","IND","EUR")
for (i in 1:length(pops)) {
	locus_name <- paste0("locus_chr_",chr,"_position_",locus,"_",pops[i])

	options(scipen = 0.1)
	#Read Pi and Tajima's D results
	d <- read.table(paste0(input_folder,"/3P-CLR/",locus_name,"/",locus_name,".windowed.pi"),h=T)
	TD <- read.table(paste0(input_folder,"/3P-CLR/",locus_name,"/",locus_name,".Tajima.D"),h=T)

	#Plot Pi and Tajima's D
	plot(d$BIN_END,d$PI,xlab="",ylab=expression(pi),type="l",lwd=2,xaxt='n',cex.lab=2,cex.axis=1.5,ylim=c(0,0.007)) #0.017 for other chromosomes
	means <- read.table(paste0(input_folder,"/stats.",chrom),header=FALSE)
	abline(h=means[2],lty=2,lwd=1)
	par(new = TRUE)
	plot(TD$BIN_START,TD$TajimaD,type="l",lwd=2,axes=FALSE,xlab="",ylab="",col="red",cex=0.75,cex.axis=2,ylim=c(-3,4))
	axis(side = 4, at = pretty(range(TD$TajimaD,na.rm=TRUE),n=3),col="red",col.ticks="red",col.axis="red",cex.axis=2)
	mtext("Taj D", col="red",side = 4, line = 3,cex=1.25,cex.axis=1.5)
	abline(h=means[3],lty=2,lwd=1,col="red")
}

#Plot gene structure
plot(d$BIN_END,d$PI,xlab="",ylab="",ylim=c(0,5.6),col="white",axes=FALSE) #Doing this just to make sure this plot and the previous one are exactly aligned. 
axis(1,at=pretty(d$BIN_END,n=8),labels=round(pretty(d$BIN_END,n=8)/1000000,2),cex.axis=2.00)
plotGenes(ann,chromstart,chromend,0)


dev.off()
