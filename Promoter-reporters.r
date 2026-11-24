#Analyze and plot promoter reporter microscopy data
setwd('C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R')

library(readxl)
AtDsRed <- read_excel("C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R/AtDsRed.xlsx")
AtDsRed <- read_excel("C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R/AtDsRed_7-15-26.xlsx")

AtDsRed = AtDsRed[11:41,]
AtDsRed[19,9] = "BM"
AtDsRed[which(AtDsRed$Stage == "Shed pollen"),9] = "Shed"
AtDsRed[which(AtDsRed$Stage == "tetrads"),9] = "Tetrad"
rownames(AtDsRed) = paste(AtDsRed$Genotype, AtDsRed$Plant, "-", AtDsRed$Flower, sep='')
stages = AtDsRed$Stage
names(stages) = rownames(AtDsRed)


#barplot(height=as.numeric(AtDsRed[5,7:11]), names.arg = c(0,1,2,3,4), main='Tetrad', xlab='Number of Fluorescent Members in Tetrad', ylab='Number of Tetrads')
#barplot(height=as.numeric(AtDsRed[8,7:11]), names.arg = c(0,1,2,3,4), main='UM', xlab='Number of Fluorescent Members in Tetrad', ylab='Number of Tetrads')
#barplot(height=as.numeric(AtDsRed[3,7:11]), names.arg = c(0,1,2,3,4), main='BM', xlab='Number of Fluorescent Members in Tetrad', ylab='Number of Tetrads')
#barplot(height=as.numeric(AtDsRed[7,7:11]), names.arg = c(0,1,2,3,4), main='TP', xlab='Number of Fluorescent Members in Tetrad', ylab='Number of Tetrads')

library(ggplot2)
library(ggpubr)



Tetrad = ggplot(data = data.frame(members=c(0,1,2,3,4), tetrads=as.numeric(AtDsRed[5,7:11])), 
    aes(x=members, y=tetrads)) + geom_bar(stat="identity") + labs(title = "Tetrad", 
    x="Number of Fluorescent Members in Tetrad", y="Number of Tetrads") + theme(panel.background = element_blank(), 
    panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
UM = ggplot(data = data.frame(members=c(0,1,2,3,4), tetrads=as.numeric(AtDsRed[8,7:11])), 
    aes(x=members, y=tetrads)) + geom_bar(stat="identity") + labs(title = "UM", 
    x="Number of Fluorescent Members in Tetrad", y="Number of Tetrads") + theme(panel.background = element_blank(), 
    panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
BM = ggplot(data = data.frame(members=c(0,1,2,3,4), tetrads=as.numeric(AtDsRed[3,7:11])), 
    aes(x=members, y=tetrads)) + geom_bar(stat="identity") + labs(title = "BM", 
    x="Number of Fluorescent Members in Tetrad", y="Number of Tetrads") + theme(panel.background = element_blank(), 
    panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
TP = ggplot(data = data.frame(members=c(0,1,2,3,4), tetrads=as.numeric(AtDsRed[7,7:11])), 
    aes(x=members, y=tetrads)) + geom_bar(stat="identity") + labs(title = "TP", 
    x="Number of Fluorescent Members in Tetrad", y="Number of Tetrads") + theme(panel.background = element_blank(), 
    panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

svg('AtDsRed_qrt1_barplots.svg', width=4, height=7)
ggarrange(Tetrad, UM, BM, TP, nrow=4, align='v')
dev.off()

AtDsRed2 = as.data.frame(AtDsRed[c(5,8,3,7,2,1),c(12,13,5)])
library(gridExtra)
svg('AtDsRed_qrt1_table.svg', width=10, height=3)
grid.arrange(tableGrob(AtDsRed2))
dev.off()



AtDsRed2 = as.data.frame(t(AtDsRed[names(stages),c(10:14)]))
colnames(AtDsRed2) = names(stages)

as.numeric(colSums(AtDsRed[names(which(stages == 'Tetrad')),10:14]))

Tetrad = ggplot(data = data.frame(members=c(0,1,2,3,4), tetrads=as.numeric(colSums(AtDsRed[names(which(stages == 'Tetrad')),10:14]))), 
    aes(x=members, y=tetrads)) + geom_bar(stat="identity") + labs(title = "Tetrad", 
    x="Number of Fluorescent Members in Tetrad", y="Number of Tetrads") + theme(panel.background = element_blank(), 
    panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
UM = ggplot(data = data.frame(members=c(0,1,2,3,4), tetrads=as.numeric(colSums(AtDsRed[names(which(stages == 'UM')),10:14]))), 
    aes(x=members, y=tetrads)) + geom_bar(stat="identity") + labs(title = "UM", 
    x="Number of Fluorescent Members in Tetrad", y="Number of Tetrads") + theme(panel.background = element_blank(), 
    panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
BM = ggplot(data = data.frame(members=c(0,1,2,3,4), tetrads=as.numeric(colSums(AtDsRed[names(which(stages == 'BM')),10:14]))), 
    aes(x=members, y=tetrads)) + geom_bar(stat="identity") + labs(title = "BM", 
    x="Number of Fluorescent Members in Tetrad", y="Number of Tetrads") + theme(panel.background = element_blank(), 
    panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
TP = ggplot(data = data.frame(members=c(0,1,2,3,4), tetrads=as.numeric(colSums(AtDsRed[names(which(stages == 'TP')),10:14]))), 
    aes(x=members, y=tetrads)) + geom_bar(stat="identity") + labs(title = "TP", 
    x="Number of Fluorescent Members in Tetrad", y="Number of Tetrads") + theme(panel.background = element_blank(), 
    panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
Shed = ggplot(data = data.frame(members=c(0,1,2,3,4), tetrads=as.numeric(colSums(AtDsRed[names(which(stages == 'Shed')),10:14]))), 
    aes(x=members, y=tetrads)) + geom_bar(stat="identity") + labs(title = "TP", 
    x="Number of Fluorescent Members in Tetrad", y="Number of Tetrads") + theme(panel.background = element_blank(), 
    panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

svg('AtDsRed_qrt1_barplots.svg', width=4, height=7)
ggarrange(Tetrad, UM, BM, TP, nrow=4, align='v')
dev.off()


