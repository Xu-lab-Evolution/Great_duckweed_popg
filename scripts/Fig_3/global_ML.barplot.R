library(ggplot2)
#install.packages("tidyverse")
library(tidyverse)
library(ggpubr)
library(rstatix)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

samp_list=args[1] 
regional_wML_table_path=args[2] 


std.error <- function(x) sd(x)/sqrt(length(x))

data_summary <- function(data, varname, groupnames){ ##########used for calculation of mean and sd
  require(plyr)

  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE) ,se = std.error(x[[col]])
      )
  }

  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# dfx <- data.frame(
#   group = c(rep('A', 8), rep('B', 15), rep('C', 6)),
#   sex = sample(c("M", "F"), size = 29, replace = TRUE),
#   age = runif(n = 29, min = 18, max = 54)
# )
#
# ddply(dfx, .(group, sex), summarize,
#       mean = round(mean(age), 2),
#       sd = round(sd(age), 2),
#       se =  round(std.error(age), 2))

#setwd("e:/project/02.duckweed_popg/GWAS_methyl/regional_ml")

samp20 <- read.table(samp_list)[,1]

wML_tab <- data.frame()

for (i in samp20){ #e:/project/02.duckweed_popg/GWAS_methyl/regional_ml/regional_wML_tsv/AM32_CX_report.regional_wML.table.tsv
  #i <- "AM32"
  df <- read.table(paste0(regional_wML_table_path,"/",i,"_CX_report.regional_wML.table.tsv"),header = T)
  wML_tab <- rbind.data.frame(wML_tab, df)
}

wML_tab$sample <- gsub("_CX_report","",wML_tab$sample)
wML_tab$pop <-  gsub("\\d+","",wML_tab$sample)

six_region_wML_tab <- wML_tab[-c(which(wML_tab$Region=="down2k"), which(wML_tab$Region=="5utr"),which(wML_tab$Region=="3utr")),]

six_region_wML_tab$Region <- factor(six_region_wML_tab$Region, levels=c( "genome_wide","up2k","gene", "exon", "intron", "TE"))
#six_region_wML_tab$wML <- as.numeric(six_region_wML_tab$wML)
#unique(population_db$group)

#data_summary(ToothGrowth, varname = "len", groupnames = c("supp", "dose"))

six_region_wML_tab.sum <- data_summary(six_region_wML_tab, varname = "wML",groupnames = c("Region", "pop","Context"))
six_region_wML_tab.sum$Context <- factor(six_region_wML_tab.sum$Context,levels = c("CpG","CHG","CHH"))

# p <- ggplot(six_region_wML_tab.sum, aes(x = Region, y = wML, fill=pop)) +
#   geom_bar(stat = "identity",  position = position_dodge()) +
#   geom_errorbar(aes(ymin = wML, ymax = wML + sd), width = 0.2, position = position_dodge(0.9)) +
#   #geom_jitter(aes(color = group))  +
#   facet_grid(Context ~ .,scales="free") +
#   scale_colour_brewer(type = "seq", palette = "Spectral") +
#   #theme_classic2() +
#   theme_classic() +
#   scale_fill_manual(values=c("#F7931E","#FF0000", "#0000FF","#39B54A"))
#   #theme(legend.position="none")
# #theme_classic()
# #geom_dotplot(binaxis='y', stackdir='center', dotsize =0.4,position=position_dodge(0.9))


# ggsave("global_ML_plot.pdf",plot =p, width = 16,
#        height = 16,
#        units = "cm")

#only genome-wide wML bar plots

p <- ggplot(six_region_wML_tab.sum[six_region_wML_tab.sum$Region=="genome_wide",], aes(x = pop, y = wML, fill=pop)) +
  geom_bar(stat = "identity",  position = position_dodge()) +
  #geom_errorbar(aes(ymin = wML - se , ymax = wML + se), width = 0.2, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = wML - se , ymax = wML + se), position = position_dodge(0.9),width=0) +
  #geom_jitter(aes(color = group))  +
  facet_grid(Context ~ .,scales="free") +
  scale_colour_brewer(type = "seq", palette = "Spectral") +
  #theme_classic2() +
  theme_classic() +
  scale_fill_manual(values=c("#F7931E","#FF0000", "#0000FF","#39B54A"))
#theme(legend.position="none")
#theme_classic()
#geom_dotplot(binaxis='y', stackdir='center', dotsize =0.4,position=position_dodge(0.9))


ggsave("./output/global_ML_plot.genome-wide.pdf",plot =p, width = 9,
       height = 13,
       units = "cm")


###############################################################stats
##############comprison among regions

pwc.cg <-six_region_wML_tab[which(six_region_wML_tab$Context=="CpG"),] %>%
  wilcox_test(wML ~ Region, p.adjust.method = "bonferroni")
pwc.cg

pwc.chg <- six_region_wML_tab[which(six_region_wML_tab$Context=="CHG"),] %>%
  wilcox_test(wML ~ Region, p.adjust.method = "bonferroni")
pwc.chg

#res.chh <- kruskal.test(mean_ML  ~ group, data = six_region_wML_tab[which(six_region_wML_tab$Context=="CHH"),] )
pwc.chh <- six_region_wML_tab[which(six_region_wML_tab$Context=="CHH"),] %>%
  wilcox_test(wML ~ Region, p.adjust.method = "bonferroni")
pwc.chh

pwc.cg <- as.data.frame(pwc.cg)
pwc.cg$con <- "CG"

pwc.chg <- as.data.frame(pwc.chg)
pwc.chg$con <- "CHG"

pwc.chh <- as.data.frame(pwc.chh)
pwc.chh$con <- "CHH"


comp_among_regions <- rbind.data.frame(pwc.cg,pwc.chg,pwc.chh)

write.table(comp_among_regions,file="./output/comp_among_regions.tsv",sep="\t",quote=F,col.names = T,row.names = F)


########comprison among populations regarding different genomic regions


#unique(all_region_ML_db$sample_name)
#colnames(all_region_ML_db) <- c(all_region_ML_db)
#six_region_wML_tab$pop <- factor(six_region_wML_tab$pop)
six_region_wML_tab[which(six_region_wML_tab$sample == "AM4"),]

write.table(six_region_wML_tab,file="./output/six_region_wML_tab.tsv",sep="\t",quote=F,col.names = T,row.names = F)

all_sig_db <- data.frame()
for (i in c("genome_wide","gene" ,"exon" ,"intron","up2k","TE")){
  for (j in c("CpG","CHG","CHH")){
    #i<- "TE"
    #j <- "CpG"
    db <- six_region_wML_tab[which(six_region_wML_tab$Context == j),]
    db <- db[which(db$Region  == i),]
    tmp_res <- db %>%
      wilcox_test(wML ~ pop, p.adjust.method = "bonferroni")

    tmp_res <- as.data.frame(tmp_res)
    tmp_res$con <- j
    tmp_res$region <- i
    all_sig_db <- rbind.data.frame(all_sig_db,tmp_res)
  }
}

write.table(all_sig_db,file="./output/comp_among_pops.tsv",sep="\t",quote=F,col.names = T,row.names = F)
