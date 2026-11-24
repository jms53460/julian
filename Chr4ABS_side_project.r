#Figuring out genes density and expression in maize pollen for regions on Ch4 between centromere (110,759,999) and ABS array (183,250,000) and from ABS array to end of chromosome
load('MaizeRemap_8-10_2025.RData')
load("Zm2.RData")
ZmChr4 = Zm_genes[Zm_genes$Chr == 4,]
ZmChr4R1 = ZmChr4[ZmChr4$Position > 110759999 & ZmChr4$Position < 183250000,] #1070 genes between centromere and ABS array
ZmChr4R2 = ZmChr4[ZmChr4$Position >= 183250000,] #1496 genes between ABS array and end of chromosome

Zm2_TP = names(Zm2_stages[which(Zm2_stages == "TP")]) #76
Zm2_TP_shed = names(Zm2_stages[which(Zm2_stages == "TP_shed")]) #23
Zm2_TPM[rownames(ZmChr4R1),Zm2_TP]

R1over10in5pTP = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R1),Zm2_TP]>10)/length(Zm2_TP)  > 0.05)) #325 genes over 10 TPM in 5% TP
R1over10in100pTP = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R1),Zm2_TP]>10)/length(Zm2_TP)  == 1)) #44 genes over 10 TPM in 100% TP

R1over100in5pTP = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R1),Zm2_TP]>100)/length(Zm2_TP)  > 0.05)) #80 genes over 100 TPM in 5% TP
R1over100in100pTP = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R1),Zm2_TP]>100)/length(Zm2_TP)  == 1)) #19 genes over 100 TPM in 100% TP

R1over10in5pTP_shed = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R1),Zm2_TP_shed]>10)/length(Zm2_TP_shed)  > 0.05)) #167 genes over 10 TPM in 5% TP_shed
R1over10in100pTP_shed = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R1),Zm2_TP_shed]>10)/length(Zm2_TP_shed)  == 1)) #31 genes over 10 TPM in 100% TP_shed

R1over100in5pTP_shed = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R1),Zm2_TP_shed]>100)/length(Zm2_TP_shed)  > 0.05)) #37 genes over 100 TPM in 5% TP_shed
R1over100in100pTP_shed = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R1),Zm2_TP_shed]>100)/length(Zm2_TP_shed)  == 1)) #16 genes over 100 TPM in 100% TP_shed



R2over10in5pTP = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R2),Zm2_TP]>10)/length(Zm2_TP)  > 0.05)) #510 genes over 10 TPM in 5% TP
R2over10in100pTP = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R2),Zm2_TP]>10)/length(Zm2_TP)  == 1)) #53 genes over 10 TPM in 100% TP

R2over100in5pTP = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R2),Zm2_TP]>100)/length(Zm2_TP)  > 0.05)) #114 genes over 100 TPM in 5% TP
R2over100in100pTP = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R2),Zm2_TP]>100)/length(Zm2_TP)  == 1)) #15 genes over 100 TPM in 100% TP

R2over10in5pTP_shed = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R2),Zm2_TP_shed]>10)/length(Zm2_TP_shed)  > 0.05)) #242 genes over 10 TPM in 5% TP_shed
R2over10in100pTP_shed = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R2),Zm2_TP_shed]>10)/length(Zm2_TP_shed)  == 1)) #53 genes over 10 TPM in 100% TP_shed

R2over100in5pTP_shed = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R2),Zm2_TP_shed]>100)/length(Zm2_TP_shed)  > 0.05)) #55 genes over 100 TPM in 5% TP_shed
R2over100in100pTP_shed = names(which(rowSums(Zm2_TPM[rownames(ZmChr4R2),Zm2_TP_shed]>100)/length(Zm2_TP_shed)  == 1)) #28 genes over 100 TPM in 100% TP_shed

plot(ZmChr4R1[R1over10in5pTP,2], rowMeans(Zm2_TPM[R1over10in5pTP,Zm2_TP]), pch=19, cex=2)
plot(ZmChr4R2[R2over10in5pTP,2], rowMeans(Zm2_TPM[R2over10in5pTP,Zm2_TP]), pch=19, cex=2)

plot(ZmChr4R1[R1over10in5pTP,2], rowMeans(Zm2_TPM[R1over10in5pTP,Zm2_TP]), pch=19, cex=2, ylim=c(0,500))
abline(h=10)

plot(ZmChr4R2[R2over10in5pTP,2], rowMeans(Zm2_TPM[R2over10in5pTP,Zm2_TP]), pch=19, cex=2, ylim=c(0,500))
abline(h=10)
