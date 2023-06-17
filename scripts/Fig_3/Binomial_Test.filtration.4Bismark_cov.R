#setwd("e:/project/02.duckweed_popg/GWAS_methyl/methyl_level/test_ML_binormial_test")
library("R.utils")
############################some pars that could be modified easily
fdr_threshold <- 0.01
min_cov <- 5
############################functions
#perform binomial test and return P
bin_func <- function(x,samp_conv_r){
  if (x[1] ==0){
    p <- 0
  } else {
    p <- binom.test(x[1],x[1]+x[2],samp_conv_r)$p.value
  }
  return(p)
}

############################main code
#dealing with input pars
args=commandArgs(T)
conversion <- read.table(args[1],header = T)
#conversion <- read.table("convs.20samp.tab",header = T)

cov_file <- basename(args[2])
cov_dir <- dirname(args[2])
#dirname("/lustre/miifs03/scratch/m2_jgu-EvolTroph/ywang/WGBS_Alexandra/X204SC21103110-Z04-F002/CX_report/IN39_CX_report.txt.nuc.gz")
sample_id <- gsub("_.+","",cov_file) #AM87_CX_report.txt.gz
#sample_id <- "AM87"
cov <- read.table(gzfile(args[2]))
#cov <- read.table(gzfile("AM87_CX_report.txt.gz"))
print(paste0("dealing with ",sample_id,", from file ",args[2],". The conversion table was specified as ",args[1]))


#############################calculate the non-conversion rate for certain sample
samp_conv <- conversion[conversion$Sample_name ==sample_id, "BS_conversion_rate"]
samp_conv <- 1-samp_conv/100


#############################perform binomial test
# p <- binom.test(cov$V4,cov$V4+cov$V5,samp_conv)$p.value
# fdr <- p.adjust(p,"fdr")
p_v <- apply(cov[cov$V4>0,c(4,5)],1,bin_func,samp_conv_r=samp_conv)
fdr <- p.adjust(p_v,"fdr")


#############################output two different files 1)original file with filtration annotation 2)filtered file
cov[which(row.names(cov) %in% names(fdr[fdr > fdr_threshold])),"Filter"] <- "Bin_fail"
cov[-which(row.names(cov) %in% names(fdr[fdr > fdr_threshold])),"Filter"] <- "PASS"
cov[which(cov$V4+cov$V5<min_cov),"Filter"] <- "Low_coverage"
cov[which(cov$V4+cov$V5==0),"Filter"] <- "Non_coverage"

print("Output anno file...")
anno.file.name <- paste0(cov_dir,"/",sample_id,"_CX_report.nuc.anno")
write.table(cov,file=anno.file.name, quote = F, sep = "\t",row.names = F, col.names = F)
gzip(anno.file.name, destname=paste0(anno.file.name,".gz"), overwrite=T, remove=T)

filtered_cov <- cov[which(cov$Filter == "PASS" | cov$Filter =="Bin_fail"),]
filtered_cov[which(filtered_cov$Filter =="Bin_fail"),4] <- 0
filtered_cov <- filtered_cov[,c(1:7)]

print("Output filtered file...")
filter.file.name <- paste0(cov_dir,"/",sample_id,"_CX_report.filtered")
write.table(filtered_cov,file=filter.file.name, quote = F, sep = "\t",row.names = F, col.names = F)
gzip(filter.file.name, destname=paste0(filter.file.name,".gz"), overwrite=T, remove=T)
