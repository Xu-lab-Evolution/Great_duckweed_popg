library("ggplot2")
library(cowplot)

############################################################2 key SVs AF piechat 

pop_col <- factor(c(rep("SE_Asia",2),rep("India",2),rep("Europe",2),rep("America",2)),level=c("America","SE_Asia","India","Europe"))
af_colour <- c( "black", "#d2d2d2") #pick colours for vis

##################agl62_ins 
agl62_ins = data.frame(subject = pop_col,
                       credit = rep(c("f","t"),4),
                       value = c(87.5,100-87.5,47.3,100-47.3,25,100-25,0,100))

agl62_ins$subject <- factor(agl62_ins$subject)
agl62_ins$credit <- factor(agl62_ins$credit) 


agl62_ins_p <- ggplot(data=agl62_ins, aes(x=" ", y=value, group=credit, colour=credit, fill=credit)) +
  geom_bar(width = 0.2, stat = "identity",color="white") +
  coord_polar("y", start=0) + 
  facet_grid(.~ subject) + theme_void() +
  scale_fill_manual(values=c(af_colour)) +
  theme(legend.position="none")+
  theme(text = element_blank())



##################SOC1_del 
SOC1_del = data.frame(subject = pop_col,
                      credit = rep(c("f","t"),4),
                      value = c(0,100,73,100-73,0,100,0,100))

SOC1_del$subject <- factor(SOC1_del$subject)
SOC1_del$credit <- factor(SOC1_del$credit) 


SOC1_del_p <- ggplot(data=SOC1_del, aes(x=" ", y=value, group=credit, colour=credit, fill=credit)) +
  geom_bar(width = 0.2, stat = "identity",color="white") +
  coord_polar("y", start=0) + 
  facet_grid(.~ subject) +theme_void() +
  scale_fill_manual(values=af_colour) +
  theme(legend.position="none")+
  theme(text = element_text(size=4)) 
  

SV2_piechart_p<- plot_grid(plotlist = list(SOC1_del_p, agl62_ins_p), nrow = 2,align = "vh")

ggsave(SV2_piechart_p,file="./output/4pop_2SV.piechart.pdf",
       width = 6,
       height = 3,
       units = "cm")

