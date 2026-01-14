setwd('C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R')
#Arabidopsis thaliana and lyrata CoGe: https://genomevolution.org/r/9fb74

#Table: https://genomevolution.org/coge/data/diags/24739/36190/24739_36190.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.ks

AtNSdata = read.table('Arabidopsis_thaliana_lyrata_CoGe.tab', skip = 2, sep = '\t', as.is = T)
AtNSdata = cbind(AtNSdata, sub("[.].*", "", sub("[|].*", "", sub(".*AT", "ID=gene-AT", AtNSdata[,4]))), sub("[|].*", '', sub(".*CDS:", "CDS:", AtNSdata[,8])), as.numeric(sub('.*[|]','',AtNSdata[,4])))
colnames(AtNSdata)[13:15] = c('thalianaGene','lyrataGene','Identity')


## Remove putative orthologs that are not reciprocal best hits ##
removes = rep(FALSE, nrow(AtNSdata))
for (sgene in as.character(unique(AtNSdata[duplicated(AtNSdata[,13]),13]))) {
	tr = which(AtNSdata[,13] == sgene)
	removes[tr[rank(-AtNSdata[tr,15]) != 1]] = TRUE
}
for (mgene in as.character(unique(AtNSdata[duplicated(AtNSdata[,14]),14]))) {
	tr = which(AtNSdata[,14] == mgene)
	removes[tr[rank(-AtNSdata[tr,15]) != 1]] = TRUE
}
AtNSdata = AtNSdata[!removes,]
##

AtNSdata = AtNSdata[AtNSdata[,15] >= 80,]  # Require at least 80% identity between orthologs when calculating dNdS
AtNSdata = AtNSdata[-which(is.na(as.numeric(AtNSdata[,1]))),] # Remove genes that lack calculated synonymous and non-synonymous rates or are undefined
AtNSdata = AtNSdata[-which(as.numeric(AtNSdata[,1]) == 0),] # Remove genes with 0 synonymous subs
AtNS = as.numeric(AtNSdata[,2]) / as.numeric(AtNSdata[,1])
names(AtNS) = AtNSdata[,13]

load("HapExpScore10k.RData")
library('ComplexHeatmap')
library(circlize)
col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))
#Heatmap(HapExpScore, col = col_fun, cluster_columns=T, cluster_rows=T, show_row_names = FALSE, show_column_names = FALSE)

pseudocount = 1*10^6/quantile(colSums(D), p = .1)
A2 = sweep(D, 2, colSums(D), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(D[,colnames(A2b)] >= 10) >= 10, ]  # Require each gene to have at least 10 UMIs in at least 10 cells

dNdSgenes = names(AtNS[which(names(AtNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d

A3 = A2d[dNdSgenes,over30k]
AtNS2 = AtNS[dNdSgenes]
plot(rowMeans(A3),AtNS2)
plot(rowMeans(A3),AtNS2, ylim=c(0,2))

MonoSamp = names(stages[grep("mono", stages)])
BiSamp = names(stages[-grep("mono", stages)])
plot(rowMeans(A3[,MonoSamp]), rowMeans(A3[,BiSamp]))
plot(rowMeans(D[dNdSgenes,MonoSamp]), rowMeans(D[dNdSgenes,BiSamp]))

#MonoGs = names(which(rowMeans(D[dNdSgenes,MonoSamp]) >= 4*rowMeans(D[dNdSgenes,BiSamp]))) # Genes with 4x mean UMIs from MonoSamp compared to BiSamp
#BiGs = names(which(rowMeans(D[dNdSgenes,BiSamp]) >= 4*rowMeans(D[dNdSgenes,MonoSamp]))) # Genes with 4x mean UMIs from BiSamp compared to MonoSamp
MonoGs = names(which(apply(D[dNdSgenes,MonoSamp], 1, max) >= 4*apply(D[dNdSgenes,BiSamp], 1, max))) # Genes with 4x max UMIs from MonoSamp compared to BiSamp
BiGs = names(which(apply(D[dNdSgenes,BiSamp], 1, max) >= 4*apply(D[dNdSgenes,MonoSamp], 1, max))) # Genes with 4x max UMIs from BiSamp compared to MonoSamp


plot(rowMeans(D[dNdSgenes,MonoSamp]), rowMeans(D[dNdSgenes,BiSamp]), pch=19, cex=3)
points(rowMeans(D[c(MonoGs,BiGs),MonoSamp]), rowMeans(D[c(MonoGs,BiGs),BiSamp]), pch=19, cex=3, col="red2")

plot(rowMeans(A3[,MonoSamp]), rowMeans(A3[,BiSamp]), pch=19, cex=3)
points(rowMeans(A3[c(MonoGs,BiGs),MonoSamp]), rowMeans(A3[c(MonoGs,BiGs),BiSamp]), pch=19, cex=3, col="red2")

boxplot(list(MonoGs = AtNS[MonoGs], BiGs = AtNS[BiGs])) # These gene groups look to have the same dN/dS

plot(apply(D[dNdSgenes,MonoSamp], 1, max), apply(D[dNdSgenes,BiSamp], 1, max), pch=19, cex=3)
points(apply(D[c(MonoGs,BiGs),MonoSamp], 1, max), apply(D[c(MonoGs,BiGs),BiSamp], 1, max), pch=19, cex=3, col="red2")

HapNonHap = list("Not or lowly expressed when haploid" = AtNS[-which(names(AtNS) %in% rownames(A2d))], "Expressed when haploid" = AtNS2)
boxplot(HapNonHap, ylim = c(0,1)) # Some outliers at really high values
#length(AtNS2) = 2567
#length(AtNS[-which(names(AtNS) %in% rownames(A2d))]) = 17074
HapNonHap2 = data.frame("NS" = c(AtNS[-which(names(AtNS) %in% rownames(A2d))], AtNS2), Class = c(rep(1, times = length(AtNS[-which(names(AtNS) %in% rownames(A2d))])), rep(2, times = length(AtNS2))))
pairwise.wilcox.test(HapNonHap2$NS, HapNonHap2$Class, p.adjust = 'holm')  # p-values
#data:  HapNonHap$NS and HapNonHap$Class
#  1
#2 <2e-16 #If this was done right, I think that means a very significant result


fano = apply(A2d, 1, var)/rowMeans(A2d)  # fano factor is a measure of gene variance
hmat = A2d[rank(-fano[rownames(A2d)]) <= 500,over30k]
minmax = function(x) {
	sweep(x - log(pseudocount, 10), 1, apply(x - log(pseudocount, 10), 1, max), '/')
}

HapFrac_col = colorRamp2(c(0.5, 1), c("white", "purple4"))
UMI_col = colorRamp2(c(4.4, 5.3), c("white", "forestgreen"))




#Rice CoGe: https://genomevolution.org/coge/SynMap.pl

#Table: https://genomevolution.org/coge/data/diags/16890/51405/16890_51405.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.ks

RiceNSdata = read.table('Rice_CoGe.tab', skip = 2, sep = '\t', as.is = T)
RiceNSdata = cbind(RiceNSdata, sub("[|].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[|].*", '', sub(".*OL", "OL", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
colnames(RiceNSdata)[13:15] = c('sativaGene','longiGene','Identity')

## Remove putative orthologs that are not reciprocal best hits ##
removes = rep(FALSE, nrow(RiceNSdata))
for (sgene in as.character(unique(RiceNSdata[duplicated(RiceNSdata[,13]),13]))) {
	tr = which(RiceNSdata[,13] == sgene)
	removes[tr[rank(-RiceNSdata[tr,15]) != 1]] = TRUE
}
for (mgene in as.character(unique(RiceNSdata[duplicated(RiceNSdata[,14]),14]))) {
	tr = which(RiceNSdata[,14] == mgene)
	removes[tr[rank(-RiceNSdata[tr,15]) != 1]] = TRUE
}
RiceNSdata = RiceNSdata[!removes,]
##

RiceNSdata = RiceNSdata[RiceNSdata[,15] >= 80,]  # Require at least 80% identity between orthologs when calculating dNdS
RiceNSdata = RiceNSdata[-which(is.na(as.numeric(RiceNSdata[,1]))),] # Remove genes that lack calculated synonymous and non-synonymous rates or are undefined
RiceNSdata = RiceNSdata[-which(as.numeric(RiceNSdata[,1]) == 0),] # Remove genes with 0 synonymous subs
RiceNS = as.numeric(RiceNSdata[,2]) / as.numeric(RiceNSdata[,1])
names(RiceNS) = RiceNSdata[,13]

install.packages('riceidconverter')
library('riceidconverter')
RiceNSdata$japonicaGene = 
RiceConvert = RiceIDConvert(RiceNSdata[,13],'MSU', toType = 'SYMBOL')
#Not a neat conversion...probably better to redo the processing using the gff that matches

load("Rice_7-8_2025_with_meta.RData")
library('ComplexHeatmap')
library(circlize)
col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))
#Heatmap(HapExpScore, col = col_fun, cluster_columns=T, cluster_rows=T, show_row_names = FALSE, show_column_names = FALSE)

pseudocount = 1*10^6/quantile(colSums(D), p = .1)
A2 = sweep(D, 2, colSums(D), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(D[,colnames(A2b)] >= 10) >= 10, ]  # Require each gene to have at least 10 UMIs in at least 10 cells

dNdSgenes = names(AtNS[which(names(AtNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d

A3 = A2d[dNdSgenes,over30k]
AtNS2 = AtNS[dNdSgenes]
plot(rowMeans(A3),AtNS2)
plot(rowMeans(A3),AtNS2, ylim=c(0,2))

MonoSamp = names(stages[grep("mono", stages)])
BiSamp = names(stages[-grep("mono", stages)])
plot(rowMeans(A3[,MonoSamp]), rowMeans(A3[,BiSamp]))
plot(rowMeans(D[dNdSgenes,MonoSamp]), rowMeans(D[dNdSgenes,BiSamp]))

#MonoGs = names(which(rowMeans(D[dNdSgenes,MonoSamp]) >= 4*rowMeans(D[dNdSgenes,BiSamp]))) # Genes with 4x mean UMIs from MonoSamp compared to BiSamp
#BiGs = names(which(rowMeans(D[dNdSgenes,BiSamp]) >= 4*rowMeans(D[dNdSgenes,MonoSamp]))) # Genes with 4x mean UMIs from BiSamp compared to MonoSamp
MonoGs = names(which(apply(D[dNdSgenes,MonoSamp], 1, max) >= 4*apply(D[dNdSgenes,BiSamp], 1, max))) # Genes with 4x max UMIs from MonoSamp compared to BiSamp
BiGs = names(which(apply(D[dNdSgenes,BiSamp], 1, max) >= 4*apply(D[dNdSgenes,MonoSamp], 1, max))) # Genes with 4x max UMIs from BiSamp compared to MonoSamp


plot(rowMeans(D[dNdSgenes,MonoSamp]), rowMeans(D[dNdSgenes,BiSamp]), pch=19, cex=3)
points(rowMeans(D[c(MonoGs,BiGs),MonoSamp]), rowMeans(D[c(MonoGs,BiGs),BiSamp]), pch=19, cex=3, col="red2")

plot(rowMeans(A3[,MonoSamp]), rowMeans(A3[,BiSamp]), pch=19, cex=3)
points(rowMeans(A3[c(MonoGs,BiGs),MonoSamp]), rowMeans(A3[c(MonoGs,BiGs),BiSamp]), pch=19, cex=3, col="red2")

boxplot(list(MonoGs = AtNS[MonoGs], BiGs = AtNS[BiGs])) # These gene groups look to have the same dN/dS

plot(apply(D[dNdSgenes,MonoSamp], 1, max), apply(D[dNdSgenes,BiSamp], 1, max), pch=19, cex=3)
points(apply(D[c(MonoGs,BiGs),MonoSamp], 1, max), apply(D[c(MonoGs,BiGs),BiSamp], 1, max), pch=19, cex=3, col="red2")

HapNonHap = list("Not or lowly expressed when haploid" = AtNS[-which(names(AtNS) %in% rownames(A2d))], "Expressed when haploid" = AtNS2)
boxplot(HapNonHap, ylim = c(0,1)) # Some outliers at really high values
#length(AtNS2) = 2567
#length(AtNS[-which(names(AtNS) %in% rownames(A2d))]) = 17074
HapNonHap2 = data.frame("NS" = c(AtNS[-which(names(AtNS) %in% rownames(A2d))], AtNS2), Class = c(rep(1, times = length(AtNS[-which(names(AtNS) %in% rownames(A2d))])), rep(2, times = length(AtNS2))))
pairwise.wilcox.test(HapNonHap2$NS, HapNonHap2$Class, p.adjust = 'holm')  # p-values
#data:  HapNonHap$NS and HapNonHap$Class
#  1
#2 <2e-16 #If this was done right, I think that means a very significant result

