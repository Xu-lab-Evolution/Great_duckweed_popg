############################################################population-wide meta plots
############################################################

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

input_dir=args[1] 
output_dir=args[2] 

calc_pop_wide_ml <- function(methOverRegion){
  #methOverRegion <- Gene.methOverRegion.CG
  pop_methOverRegion <- data.frame()
  methOverRegion$pop <- gsub("\\d+","",methOverRegion$sample_name)
  for (p in unique(methOverRegion$pop)){
    for (bin in unique(methOverRegion$bin_num )){
      #p <- "AM"
      #bin <- 1
      tmp_df <- methOverRegion[intersect(which(methOverRegion$pop == p),which(methOverRegion$bin_num == bin)),]
      tmp_pop_df <-data.frame(region=tmp_df$region[1],bin_num =bin, Methylation_level=mean(tmp_df$Methylation_level), pop=p)
      pop_methOverRegion <- rbind.data.frame(pop_methOverRegion, tmp_pop_df)
    }
  }
  return(pop_methOverRegion)
}

#tab<-read.table(meth,head=T);
meta_plot <- function(tab,con,region,out_d){
  #con <- "CHH" 
  #region <- "TE"
  #tab <- TE.methOverRegion.CHH.pop
  
  tab <- tab[order(tab[,2]),]
  
  adjustXaxis <- 10
  
  min <- min(tab$bin_num)  ## by default the lowest value is -19
  max <- abs(max(tab$bin_num))
  
  flank = paste( -(min -1)/adjustXaxis, "kb", sep = " ")
  xlab <- region
  
  ymin <- 0
  ymax <- max(tab$Methylation_level) + min(tab$Methylation_level)
  
  p <- ggplot(tab, aes(x=bin_num, y=Methylation_level, group=pop, col=pop)) +
    geom_line() + 
    xlab(xlab) +
    scale_x_continuous(breaks=c(min, max), labels=c(paste0("-",flank), flank))+
    theme(legend.title=element_blank()) +
    geom_vline(xintercept = c(1, max + min - 1), linetype = "dashed")+
    #expand_limits(y=0) +
    ylab("Methylation level") +
    #ylim(ymin,ymax)+
    theme_classic() +
    #scale_y_continuous(limits = c(0, ymax),n.breaks = 2) +
    scale_y_continuous(limits = c(0, ymax),n.breaks = 3) +
    scale_colour_manual(values=c("#F7931E","#FF0000", "#0000FF","#39B54A")) 
  
  fig <- paste0(out_d,"/",region,".",con,".meta_plot.pdf")
  ggsave(fig, p, height =5, width =8, unit="cm")
}

#################################################Gene meta plots

#setwd("e:/project/02.duckweed_popg/GWAS_methyl/ViewBS/methOver.Gene")


Gene.methOverRegion.CG <- read.table(paste0(input_dir,"/methOver.Gene/","20samp.Gene.CG_MethOverRegion_CG.txt"),header=T)
Gene.methOverRegion.CHG <- read.table(paste0(input_dir,"/methOver.Gene/","20samp.Gene.CHG_MethOverRegion_CHG.txt"),header=T)
Gene.methOverRegion.CHH <- read.table(paste0(input_dir,"/methOver.Gene/","20samp.Gene.CHH_MethOverRegion_CHH.txt"),header=T)

Gene.methOverRegion.CG.pop <- calc_pop_wide_ml(Gene.methOverRegion.CG)
Gene.methOverRegion.CHG.pop <- calc_pop_wide_ml(Gene.methOverRegion.CHG)
Gene.methOverRegion.CHH.pop <- calc_pop_wide_ml(Gene.methOverRegion.CHH)

meta_plot(tab=Gene.methOverRegion.CG.pop,con="CpG",region="Gene",out_d=output_dir)
meta_plot(tab=Gene.methOverRegion.CHG.pop,con="CHG",region="Gene",out_d=output_dir)
meta_plot(tab=Gene.methOverRegion.CHH.pop,con="CHH",region="Gene",out_d=output_dir)


##################################################TE meta plots

#setwd("e:/project/02.duckweed_popg/GWAS_methyl/ViewBS/methOverTE")

TE.methOverRegion.CG <- read.table(paste0(input_dir,"/methOver.TE/","20samp.CG_MethOverRegion_CG.txt"),header=T)
TE.methOverRegion.CHG <- read.table(paste0(input_dir,"/methOver.TE/","20samp.CHG_MethOverRegion_CHG.txt"),header=T)
TE.methOverRegion.CHH <- read.table(paste0(input_dir,"/methOver.TE/","20samp.CHH_MethOverRegion_CHH.txt"),header=T)

TE.methOverRegion.CG.pop <- calc_pop_wide_ml(TE.methOverRegion.CG)
TE.methOverRegion.CHG.pop <- calc_pop_wide_ml(TE.methOverRegion.CHG)
TE.methOverRegion.CHH.pop <- calc_pop_wide_ml(TE.methOverRegion.CHH)

meta_plot(tab=TE.methOverRegion.CG.pop,con="CpG",region="TE",out_d=output_dir)
meta_plot(tab=TE.methOverRegion.CHG.pop,con="CHG",region="TE",out_d=output_dir)
meta_plot(tab=TE.methOverRegion.CHH.pop,con="CHH",region="TE",out_d=output_dir)
