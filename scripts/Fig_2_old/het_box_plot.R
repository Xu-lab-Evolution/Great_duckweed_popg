#!/usr/bin/env Rscript

library(ggplot2)
# library(ggsignif)
# library(ggpubr)
library(stats)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop(paste0("USAGE: Rscript het_box_plot.R family_info het_file.\n"), call.=FALSE)
} 

family_info <- args[1]
het_file <- args[2]

dir <- dirname(het_file)

print(dir)

##########################read in pop info
pop_info <- read.table(family_info,header = T,comment.char = "")

#pop_info <- read.table("../../data//clonal_family/family_info.pop_full_name.txt",header = T,comment.char = "")

simp_pop_info <- pop_info[,c("Sample_ID", "Population", "Colour")]
rownames(simp_pop_info) <- simp_pop_info$Sample_ID

##########################read in het df from VCFtools's result
het_df <- read.table(het_file,header = T)
#het_df <- read.table("/lustre/miifs01/project/m2_jgu-evoltroph/ywang/02.popg_duckweed/228_Sp_popg/data/Fig_2/01.genome_wide_het/output_het.MAF0.05.het",header = T)
colnames(het_df) <- c("INDV","O_HOM","E_HOM","N_SITES","F")

##########################calculate the het rates:
het_df$HET_R  <- (het_df$N_SITES - het_df$O_HOM)/het_df$N_SITES
het_df_for_vis <- data.frame(het_df$INDV, het_df$HET_R)
colnames(het_df_for_vis) <- c("INDV","Het_R")
rownames(het_df_for_vis) <- het_df_for_vis$INDV
#het_df_for_vis <- het_df_for_vis[,-1]
rownames(het_df_for_vis) <- toupper(rownames(het_df_for_vis))

#rownames(het_df_for_vis) %in% pop_info$Sample_ID
het_df_for_vis <- cbind.data.frame(het_df_for_vis,simp_pop_info[rownames(het_df_for_vis),c("Population", "Colour")])
het_df_for_vis$Population <- factor(het_df_for_vis$Population,levels=c("America", "India","Europe","SE_Asia"))

##########################wilcoxon test for pair-wise comprisons:
comr_res <- pairwise.wilcox.test(het_df_for_vis$Het_R, het_df_for_vis$Population, paired=F, p.adjust.method = "fdr")
#write.table(comr_res,paste0(dir,"/het_df.MAF_0_05.wilcox_test.txt"),quote = F)
write.table(comr_res$p.value,paste0(dir,"/het_df.MAF_0_05.wilcox_test.txt"),quote = F)
##########################for visualization
##setting y axis
ymax <- max(het_df_for_vis$Het_R) * 1.1
##setting significance pos
max_values <- aggregate(Het_R ~ Population, data = het_df_for_vis, FUN = max)
sig_y <- max_values$Het_R * 1.05
sig_x <- 1:4
value <- c("a","b","c","d")
sig_df<-data.frame(sig_x,sig_y,value)


ggplot(data = het_df_for_vis, aes(x = Population, y = Het_R)) +
  geom_boxplot(aes(group = Population, fill = Population), size = 0.2,outlier.alpha =0,colour="black")+
  geom_jitter(width=0.15,color="gray",size=0.4)+
  #geom_dotplot(binaxis='y', stackdir='center',fill="gray",color="gray",alpha=0.8,dotsize=0.6)+
  #geom_point(data=Sp_df_for_plot,aes(x=group, y=genome_wide.wML),color="red",shape=19,size=2.5) +
  #geom_text(data=Sp_df_for_plot,aes(x=group, y=genome_wide.wML,label=Sp_df_for_plot$genome_wide.wML),size=2 ,color="red",vjust=1.4, hjust=-0.4) +
  theme_classic() +
  theme(legend.position="none",) +
  ylab("Genome-wide heterozygosity rate")+
  xlab("") +
  ylim(0,ymax)+
  scale_y_continuous(n.breaks=3)+
  scale_fill_manual(values=c("#F7931E", "#39B54A","#0000FF","#FF0000" )) +
  geom_text(data=sig_df, aes(x = sig_x, y = sig_y, label = value),hjust = 0.5, vjust = 0.5, size = 4)

ggsave(paste0(dir,"/MAF0_05.het_R.pdf"),width=7,height=7,units = "cm")

#write.table(het_df_for_vis,"het_df.MAF_0_2.tsv",quote = F, sep = "\t",row.names = F)
write.table(het_df_for_vis,paste0(dir,"/het_df.MAF_0_05.tsv"),quote = F, sep = "\t",row.names = F)
