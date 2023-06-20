my_predict <- function(res_df,ppChr) {
	genetic_position <- numeric(length(ppChr$V2))
	for (i in 1:length(ppChr$V2)) {
		upperLim <- res_df$dc_norm[min(which(res_df$data.Loci>ppChr$V2[i]))] #Upper lim genetic
		lowerLim <- res_df$dc_norm[max(which(res_df$data.Loci<ppChr$V2[i]))] #Lower lim genetic
		genetic_position[i] <- mean(c(lowerLim,upperLim))
	}
	return(genetic_position)
}

#Read recombination rates calculated by LDhat for AME (outgroup for 3P-CLR), chr N
data <- read.table("/home/pablo/sciebo/Muenster/Duckweed/Files_for_Demographics/recombination_whole_chromosomes/stat_interval_ChrS20_AME_rates.txtres.txt",h=T)
data <- data[-1,] #Remove first line with the weird value.
Ne <- 400000 #Ne of AME
r <- data$Mean_rho/(2*Ne) #Probability of recombination at a given site.
dc <- 50*log(1/(1-2*r)) #Formula for distance in centimorgans (taken from Wikipedia)

current_chr <- 20
pp <- read.table("countsAME.frq.count",h=F,skip=1)
ppChr <- subset(pp,pp$V1==current_chr)

chrSizes <- read.table("~/Documents/Duckweed/chromSizes.txt",header=FALSE)
size_chr <- chrSizes$V2[current_chr] #Change here to the size of the current chromosome.
pos_last_snp <- ppChr$V2[length(ppChr$V2)] #Position of last SNP of the current chromosome (check in file countsAME.frq.count).
pos_first_snp <- ppChr$V2[1]

dc_cumsum <- cumsum(dc) #Cumulative distance in centimorgans.
dc_norm <- dc_cumsum/max(dc_cumsum) - (1-pos_last_snp/size_chr) #Normalized distance.
dc_norm[dc_norm<0] <- 0

f <- 1-(max(data$Loci)-min(data$Loci))/max(data$Loci) #Relative position of first SNP.
dc_norm <- dc_norm-f
dc_norm[dc_norm<0] <- 0

res_df <- data.frame(data$Loci,dc_norm)
res_df <- rbind(c(1,0),res_df) 
res_df <- rbind(res_df,c(size_chr,1)) 

plot(res_df$data.Loci,res_df$dc_norm)
points(res_df$data.Loci,(res_df$data.Loci-pos_first_snp)/pos_last_snp,col="red")

#y <- dc_norm
#x <- data$Loci
#fit <- lm(y~x)
#fit2 <- lm(y~poly(x,2,raw=TRUE))
#fit3 <- lm(y~poly(x,3,raw=TRUE))
#fit4 <- lm(y~poly(x,4,raw=TRUE))
#fit <- lm(y~poly(x,6,raw=TRUE))
#xx=seq(1,x[length(x)])
#lines(xx, predict(fit, data.frame(x=xx)), col="red")
#lines(xx, predict(fit2, data.frame(x=xx)), col="green")
#lines(xx, predict(fit3, data.frame(x=xx)), col="blue")
#lines(xx, predict(fit4, data.frame(x=xx)), col="purple")
#lines(xx, predict(fit, data.frame(x=xx)), col="blue")

#Get physical position from input file and transform them to genetic positions
#gen_positions <- predict(fit, data.frame(x=ppChr$V2))
gen_positions <- my_predict(res_df,ppChr)
lines(ppChr$V2, gen_positions, col="green")

if (max(gen_positions)>1) {
	gen_positions_norm <- gen_positions/max(gen_positions) #Make sure range is between 0 and 1.
} else {
	gen_positions_norm <- gen_positions
}
gen_positions_norm[gen_positions_norm<0] <- 0
lines(ppChr$V2, gen_positions, col="purple")

write.table(data.frame(ppChr$V2,gen_positions_norm),file=paste0("genetic_positions_chr",current_chr,".txt"),row.names=FALSE)


