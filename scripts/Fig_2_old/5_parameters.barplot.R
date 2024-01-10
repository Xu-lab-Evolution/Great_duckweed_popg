library("ggplot2")
library(cowplot)


# tempdir()
# dir.create(tempdir())
#######################################################5 population genomics parameters barplot

#het_df <- read.table("e:/project/02.duckweed_popg/calc_genome_wide_het/het_df.MAF_0_05.tsv",header = T,comment.char = "")
het_df <- read.table("./output/het_df.MAF_0_05.tsv",header = T,comment.char = "")
het_median_df <- aggregate(het_df[,2], by = list(het_df[, 3]), FUN = median)
rownames(het_median_df) <- het_median_df$Group.1

###########data collected from different calculations 
df <- data.frame(pop=factor(c("America","SE_Asia","India","Europe"),level=c("America","SE_Asia","India","Europe")),pi=c(0.00045,0.0011,0.00061,0.00071),selection=c(0.42,0.37,0.41,0.38),LD=c(33,12,61,100),rec=c(3.158,7.139,4.002,3.142))

rownames(df) <- df$pop

df <- cbind.data.frame(df,het_median_df[rownames(df),2])

####polishing df
colnames(df)[length(colnames(df))] <- "het"
df_vis <- df[,c("pop","pi","selection","LD","rec"  ,"het" )]

#df_vis <- read.table("df_vis.tsv",header = T,stringsAsFactors = T)
pop_fill_color <- c('#F7931E','#FF0000','#39B54A','#0000FF')

pi_p <- ggplot(df_vis[,c("pop","pi")], aes(pop,pi,fill=pop)) +
  geom_bar(stat="identity", width=0.7)+
  theme_classic()+
  #scale_y_discrete(expand =c(0,0))+
  #theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  #expand_limits(y=0)+
  scale_fill_manual(values=pop_fill_color) +
  scale_y_continuous(limits = c(0,max(df_vis$pi)),breaks=c(0,max(df_vis$pi)) ,expand=c(0,0))+
  xlab("") +
  ylab("Π") +
  theme(axis.title.y = element_text(vjust= 0.5,angle=0,size=8))+
  theme(axis.text.y = element_text(size=6))+
  #theme(legend.position="none",axis.ticks.x = element_blank(),axis.text = element_blank())
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x = element_blank()) 

LD_p <-  ggplot(df_vis[,c("pop","LD")], aes(pop,LD,fill=pop))+
  geom_bar(stat="identity", width=0.7)+
  theme_classic()+
  #scale_y_discrete(expand =c(0,0))+
  #theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  #expand_limits(y=0)+
  scale_fill_manual(values=pop_fill_color) +
  scale_y_continuous(limits = c(0,max(df_vis$LD)),breaks=c(0,max(df_vis$LD)) ,expand=c(0,0))+
  xlab("") +
  ylab("LD") +
  theme(axis.title.y = element_text(vjust= 0.5,angle=0,size=8))+
  theme(axis.text.y = element_text(size=6))+
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x = element_blank()) 

rec_p <- ggplot(df_vis[,c("pop","rec")], aes(pop,rec,fill=pop))+
  geom_bar(stat="identity", width=0.7)+
  theme_classic()+
  #scale_y_discrete(expand =c(0,0))+
  #theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  #expand_limits(y=0)+
  scale_fill_manual(values=pop_fill_color) +
  scale_y_continuous(limits = c(0,8),breaks=c(0,8) ,expand=c(0,0))+
  xlab("") +
  ylab("r") +
  theme(axis.title.y = element_text(vjust= 0.5,angle=0,size=8))+
  theme(axis.text.y = element_text(size=6))+
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x = element_blank()) 

#summary(df_vis[,c("pop","selection")])
selection_p <- ggplot(df_vis[,c("pop","selection")], aes(pop,selection,fill=pop))+
  geom_bar(stat="identity", width=0.7)+
  theme_classic()+
  #scale_y_discrete(expand =c(0,0))+
  #theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  #expand_limits(y=0)+
  scale_fill_manual(values=pop_fill_color) +
  coord_cartesian(ylim=c(0.35, 0.45)) +
  scale_y_continuous(breaks=c(0.35,0.45) ,expand=c(0,0))+
  xlab("") +
  ylab("ΠN/ΠS") +
  theme(axis.title.y = element_text(vjust= 0.5,angle=0,size=8))+
  theme(axis.text.y = element_text(size=6))+
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x = element_blank()) 

  
  
het_p <- ggplot(df_vis[,c("pop","het")], aes(pop,het,fill=pop))+
  geom_bar(stat="identity", width=0.7)+
  theme_classic()+
  #scale_y_discrete(expand =c(0,0))+
  #theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  #expand_limits(y=0)+
  scale_fill_manual(values=pop_fill_color) +
  scale_y_continuous(limits = c(0,0.23),breaks=c(0,0.23) ,expand=c(0,0))+
  xlab("") +
  ylab("H") +
  theme(axis.title.y = element_text(vjust= 0.5,angle=0,size=8))+
  theme(axis.text.x = element_text(size=8),axis.ticks.x = element_blank())+
  theme(axis.text.y = element_text(size=6))+
  #theme(legend.position="none",axis.text.x = element_blank()) 
  theme(legend.position="none") 


par5_p <- plot_grid(plotlist = list(pi_p,NULL,selection_p,NULL,LD_p,NULL,rec_p,NULL,het_p), 
                    nrow = 9,
                    rel_heights=c(rep(c(1,0.01),4),1),
                    align = "vh",
                    scale=1.1)

# ggsave(par5_p,file="4pop_5pars.barplot.png",
#        width = 6,
#        height = 10,
#        units = "cm")

ggsave(par5_p,file="./output/4pop_5pars.barplot.pdf",
       width = 6,
       height = 10,
       units = "cm")

