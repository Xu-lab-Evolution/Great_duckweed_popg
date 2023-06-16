library(ggplot2)
#install.packages("scatterplot3d")
library("scatterplot3d")

rm(list=ls())
args <- commandArgs(trailingOnly = FALSE)

#setting working dir, which should be the same dir as where the script is located
script.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
script.dir <- file.path(script.dir, "") 

#args[1]
setwd(script.dir)

getwd()

###read in PCA res
pcs <- read.table( "../../../data/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033.pca.eigenvec",stringsAsFactors = F )
percent <- read.table( "../../../data/Fig_1/00.Population_Structure_and_PCA/SP_228.basic_set.snp.HWE001.LD033.pca.eigenval",stringsAsFactors = F)

percent[,1] <- as.numeric(sprintf("%.2f", percent[,1]))
prec.df <- data.frame(x=paste0("PC",1:20), y=percent[,1])
prec.df$x <- factor(prec.df$x , levels = paste0("PC",1:20))

dir.create("./output/", showWarnings = F)

###bar plot of eigenvalue 
eigenval <- ggplot(data=prec.df[1:10,], aes(x=x, y=y)) +
  geom_bar(stat="identity" )+
  geom_col(width = 0.7)+labs(x= "Principal components",y= "Eigenvalue (%)") 
#scale_y_continuous(limits = c(0,20)) 
# theme_minimal() + labs(x= "Principal components",y= "Eigenvalue (%)")

ggsave(eigenval,filename = "./output/SP228_PCA_eigenval.pdf",width = 12,height = 9,units = "cm")


###PCA 3D
rownames(percent) <- paste0("PC",1:20)

group_info <- read.table( file="../../../data/Fig_1/00.Population_Structure_and_PCA/PCA_group_info", sep="\t", header =T, row.names =1,stringsAsFactors = F,fill = T,comment.char = "")

pcs$Group <- group_info[pcs[,2],2]
pcs$Group <- factor(pcs$Group, levels = c("America", "Asia", "Europe","India"))
pcs[rownames(pcs[pcs$Group=="America",]),"color"] <- "#F7931E"
pcs[rownames(pcs[pcs$Group=="Asia",]),"color"] <- "#FF0000"
pcs[rownames(pcs[pcs$Group=="Europe",]),"color"] <- "#0000FF"
pcs[rownames(pcs[pcs$Group=="India",]),"color"] <- "#39B54A"

rownames(pcs) <- pcs[,2]
pcs <- pcs[,c(-1,-2)]
colnames(pcs) <- c(paste0("PC",1:20),"Group","color")

var <- c("PC1","PC2","PC3")

l_pc1=paste0(var[1]," ","(",percent[var[1],1],"%)")
l_pc2=paste0(var[2]," ","(",percent[var[2],1],"%)")
l_pc3=paste0(var[3]," ","(",percent[var[3],1],"%)")

pdf(file="./output/SP228_PCA_3d.pdf")

#scatterplot3d(pcs[,1:3], angle=135, pch=20,type="h",color=pcs$color);
scatterplot3d(pcs[,1:3], angle=135,pch=1,color=pcs$color,xlab=l_pc1, ylab=l_pc2, zlab=l_pc3,tick.marks=F);
legend("right", legend = levels(pcs$Gr),
       col =  c("#F7931E","#FF0000","#0000FF","#39B54A"), pch = 1, inset = 0.1)

dev.off()



