library(methylKit)
#library(IRanges)
library(regioneR)
#setwd("e:/project/02.duckweed_popg/GWAS_methyl/check_methylome_Data/calc_ML_using_methylKit")

args = commandArgs(trailingOnly=TRUE)
if (length(args)< 2) {
  #stop("At least 4 argument must be supplied (<cx_report_file> <sample_name> <tmp_path> <out_path>).n", call.=FALSE)
  stop("At least 3 argument must be supplied (<cx_report_file> <path_of_annotation> <out_path>).n", call.=FALSE)
}
cx_report=args[1] #cx_report="e:/project/02.duckweed_popg/GWAS_methyl/check_methylome_Data/calc_ML_using_methylKit/test.bis.chr01.txt"
#sample_id=args[2] #sample_id="test"
#tmp_path=args[3] #tmp_path="e:/project/02.duckweed_popg/GWAS_methyl/check_methylome_Data/calc_ML_using_methylKit"
anno=args[2]
out_path=args[3] #out_path="e:/project/02.duckweed_popg/GWAS_methyl/check_methylome_Data/calc_ML_using_methylKit"


sample_id <- basename(cx_report)
sample_id <- strsplit(sample_id,"\\.")[[1]][1]

tmp_path <- out_path
tmp_db <- paste0(tmp_path,"/methylDB")



##################################################funcs

CX_report_to_GR <-function(file,sample,dbdir,context){
  DB=methRead(#"test.bis.chr01.txt",
    file,
    sample.id=sample,
    assembly="Sp7498v2",
    treatment="0",
    context=context,
    dbtype = "tabix",
    pipeline = "bismarkCytosineReport",
    dbdir = dbdir,
    mincov = 4
  )
  GR <- as(DB,"GRanges")
  return(GR)
}

calc_ML <- function(GR){
  ML <- sum(GR$numCs)/(sum(GR$numCs) + sum(GR$numTs))
  return(ML)
}

calc_ML_mC <- function(GR){
  GR <- mC_GR
  mC <- subset(GR,GR$numCs > 0)
  #mC_ML <- mean(mC$numCs/(mC$numCs+mC$numTs))
  mCs <- length(mC$numCs)
  Cs <- length(GR$numCs)
  mC_ML <- sum(mC$numCs)/(sum(mC$numCs) + sum(mC$numTs))
  prop <- mCs/Cs
  return(c(mC_ML,mCs,Cs,prop))
}



bed_to_GR <- function(file){
  df <- read.table(file)[,c(1:3)]
  dd <- toGRanges(df)
  return(dd)
}

ML_df <- data.frame()
mC_prop <- data.frame()

for (con in c("CpG","CHG","CHH")){
  #con <- "CpG"
  mC_GR <- CX_report_to_GR(file=cx_report,sample=sample_id,dbdir=tmp_db, context = con)
  all_ml <- calc_ML(mC_GR)
  mC_ML <- calc_ML_mC(mC_GR)
  ML_df <- rbind.data.frame(ML_df,data.frame(Context=con,Region="genome_wide",wML=all_ml))
  mC_prop <- rbind.data.frame(mC_prop,data.frame(Context=con,Region="genome_wide",mC_mean_ML=mC_ML[1],mCs=mC_ML[2],Cs=mC_ML[3],mC_prop=mC_ML[4]))
  for (r in c("gene","exon","intron","5utr","3utr","up2k","down2k","TE")){
    #r <- "TE"
    #region_bed_path="e:/project/02.duckweed_popg/anno/SpGA2022_anno/SpGA2022.3/for_submission/4_methyl/SpGA2022."
    bed <- paste0(anno,"/SpGA2022.",r,".bed")
    region_GR <- bed_to_GR(bed)
    ovlp <- overlapRegions(mC_GR,region_GR,only.boolean=T)
    ml <- calc_ML(mC_GR[ovlp,])
    ML_df <- rbind.data.frame(ML_df,data.frame(Context=con,Region=r,wML=ml))
  }
}

ML_df$sample <- sample_id
mC_prop$sample <- sample_id

out_file=paste0(out_path,"/",sample_id,".regional_wML.table.tsv")
mC_prop_out_file=paste0(out_path,"/",sample_id,".global_mC_prop.table.tsv")

r_file=paste0(out_path,"/",sample_id,".regional_wML.Rdata")

write.table(ML_df,file= out_file, quote = F, sep = "\t",
            na = "NA", row.names = F,
            col.names = TRUE)

write.table(mC_prop,file= mC_prop_out_file, quote = F, sep = "\t",
            na = "NA", row.names = F,
            col.names = TRUE)

#write.table(ML_df,file= paste0(sample_id,".regional_wML.table.tsv"), quote = F, sep = "\t",na = "NA", row.names = F,col.names = TRUE)
#save.image(r_file)
