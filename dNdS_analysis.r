setwd('C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R')
#Arabidopsis thaliana and lyrata CoGe: https://genomevolution.org/r/9fb74
#Arabidopsis thaliana Col-0 (thale cress) (v10, id24739): masked repeats 50x
#Arabidopsis lyrata (EnsemblPlantsv1) (vEnsemblPlantsv1, id36190): Masked by source        ###5 MYA
#Capsella rubella (vCaprub_01, id56621): masked using RepeatMasker      ###10-14 MYA
#Brassica rapa (v1.5, id28890): NCBI WindowMasker (Hard)     ###17-20 MYA

#Move forward with one of the three options:
    #Arabidopsis lyrata comparison
    AtNSdata = read.table('Arabidopsis_thaliana_lyrata_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    AtNSdata = cbind(AtNSdata, sub("[.].*", "", sub("[|].*", "", sub(".*AT", "ID=gene-AT", AtNSdata[,4]))), sub("[|].*", '', sub(".*CDS:", "CDS:", AtNSdata[,8])), as.numeric(sub('.*[|]','',AtNSdata[,4])))
    colnames(AtNSdata)[13:15] = c('thalianaGene','lyrataGene','Identity')

    #Capsella rubella comparison
    AtNSdata = read.table('Arabidopsis_thaliana_Capsella_rubella_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    AtNSdata = cbind(AtNSdata, sub("[.].*", "", sub("[|].*", "", sub(".*AT", "ID=gene-AT", AtNSdata[,4]))), sub("[|].*", '', sub(".*cds-", "cds-", AtNSdata[,8])), as.numeric(sub('.*[|]','',AtNSdata[,4])))
    colnames(AtNSdata)[13:15] = c('thalianaGene','CapsellaGene','Identity')

    #Brassica rapa comparison
    AtNSdata = read.table('Arabidopsis_thaliana_Brassica_rapa_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    AtNSdata = cbind(AtNSdata, sub("[.].*", "", sub("[|].*", "", sub(".*AT", "ID=gene-AT", AtNSdata[,4]))), sub("[|].*", '', sub(".*Bra", "Bra", AtNSdata[,8])), as.numeric(sub('.*[|]','',AtNSdata[,4])))
    colnames(AtNSdata)[13:15] = c('thalianaGene','BrassicaGene','Identity')

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
#Oryza sativa japonica (Rice) (v7, id16890): masked
#Oryza longistaminata (Oryza longistaminata acc. IRGC110404) (vV2.1, id51405): masked using RepeatMasker        ###2 Mya
#Leersia perrieri (Gramene v1.4 Barbazuk Lab AS Project) (v1.4, id60020): unmasked      ###20 Mya
#Brachypodium distachyon (Brachypodium_distachyon_v3.0.dna_rm.toplevel.fa) (vv3, id52736): masked       ###50 Mya
##Remapping rice using files from https://rice.uga.edu/download_osa1r7.shtml

#Move forward with one of the three options:
    RiceNSdata = read.table('Rice_Oryza_longistaminata_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    RiceNSdata = cbind(RiceNSdata, sub("[.].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[|].*", '', sub(".*OL", "OL", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
    colnames(RiceNSdata)[13:15] = c('RiceGene','longiGene','Identity')

    RiceNSdata = read.table('Rice_Leersia_perrieri_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    RiceNSdata = cbind(RiceNSdata, sub("[.].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[|].*", '', sub(".*OL", "OL", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
    colnames(RiceNSdata)[13:15] = c('RiceGene','LeersiaGene','Identity')

    RiceNSdata = read.table('Rice_Brachypodium_distachyon_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    RiceNSdata = cbind(RiceNSdata, sub("[.].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[|].*", '', sub(".*OL", "OL", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
    colnames(RiceNSdata)[13:15] = c('RiceGene','BrachyGene','Identity')    

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
load("Rice_RedoAnnotations.RData")
library('ComplexHeatmap')
library(circlize)
col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))
#Heatmap(HapExpScore, col = col_fun, cluster_columns=T, cluster_rows=T, show_row_names = FALSE, show_column_names = FALSE)

pseudocount = 1*10^6/quantile(colSums(D2), p = .1)
A2 = sweep(D2, 2, colSums(D2), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(D2[,colnames(A2b)] >= 3) >= 10, ]  # Require each gene to have at least 10 UMIs in at least 10 cells

dNdSgenes = names(RiceNS[which(names(RiceNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d

over30k = names(which(colSums(D) >25000))
over30k = over30k[which(Rice_meta[over30k,12] == "N")] #getting rid of the 2 no cell controls


A3 = A2d[dNdSgenes,over30k]
RiceNS2 = RiceNS[dNdSgenes]
plot(rowMeans(A3),RiceNS2)
plot(rowMeans(A3),RiceNS2, ylim=c(0,2))


BINR = function (xx, bin = 10^6) 
{
    bin = as.numeric(Rice_genes[, 1]) * 10^6 + round(Rice_genes[, 2]/bin)
    out = by(xx, bin, colSums)
    out2 = t(matrix(unlist(out), nrow = ncol(g1)))
    colnames(out2) = colnames(g1)
    rownames(out2) = names(out)
    return(out2)
}

library(ggplot2)
library(ggpubr)
g1_bin = BINR(g1)
g2_bin = BINR(g2)
AlleleFrac_bin = g1_bin/(g1_bin + g2_bin)
AlleleFrac_bin[(g1_bin+g2_bin) < 10] = NA #remove bins with <10 genoinformative transcripts
R_binUse = which(abs(rowMeans(AlleleFrac_bin, na.rm=T) - .5) < .4)  # Exclude bins with >90% of all transcripts mapping to the same allele across all samples
AlleleFrac_bin[-R_binUse,] = NA
FracMono_all = 100*colMeans(abs(AlleleFrac_bin - .5) >= .3, na.rm=T)


MonoSamp = names(which(FracMono_all[over30k] > 10))
BiSamp = names(which(FracMono_all[over30k] < 10))
plot(rowMeans(A3[,MonoSamp]), rowMeans(A3[,BiSamp]))
plot(rowMeans(D2[dNdSgenes,MonoSamp]), rowMeans(D2[dNdSgenes,BiSamp]))

NotInBi = names(which(rowSums(D2[dNdSgenes,BiSamp] > 0) == 0))
MonoGs = names(which(rowSums(D2[NotInBi,MonoSamp] > 0) >= 3))

MonoGs = names(which(rowSums(D2[dNdSgenes,MonoSamp] > 10) >= 3))
NotInMono = names(which(rowMeans(D2[dNdSgenes,MonoSamp] > 5) == 0))
BiGs = names(which(rowSums(D2[NotInMono,BiSamp] > 0) >= 3))

RiceStages = Rice_meta[over30k,9]
names(RiceStages) = over30k
RiceStages[names(RiceStages[grep("tetrad", RiceStages)])] = "Tetrad"
RiceStages[names(which(RiceStages[BiSamp] == "UM"))] = "UM_bi"
RiceStages[names(which(RiceStages[MonoSamp] == "UM"))] = "UM_mono"

meanExp = matrix(0, ncol = length(unique(RiceStages)), nrow = nrow(A2))
rownames(meanExp) = rownames(A2)
colnames(meanExp) = c('Tetrad', 'UM_bi', 'UM_mono', 'UM/BM', 'BM', 'BM/Tri', 'Tri')
for (stage in colnames(meanExp)) {
	meanExp[,stage] = rowMeans(A2[,RiceStages == stage])
}
meanExp = meanExp[dNdSgenes,]

MonoGs = rownames(meanExp)[rowSums(meanExp[,3:7] >= 100) > 0]
MonoGs = names(which(apply(meanExp[MonoGs,3:7], 1, max)/apply(meanExp[MonoGs,1:2], 1, max) >= 2))
BiGs = rownames(meanExp)[(rowSums(meanExp[,1:2] >= 100) > 0) & (rowSums(meanExp[,3:7] >= 100) == 0)]


IncGs = rownames(meanExp)[rowSums(meanExp[,1:7] >= 100) > 0]
meanExp2 = meanExp[IncGs,]
MonoGs = rownames(meanExp2)[rowSums(meanExp2[,3:7] >= 100) > 0]
BiGs = rownames(meanExp2)[rowSums(meanExp2[,3:7] >= 100) == 0]


MonoGs = rownames(meanExp)[rowSums(meanExp[,3:7] >= 0) > 0]
MonoGs = names(which(apply(meanExp[MonoGs,3:7], 1, max) > meanExp[,2])) #used this for HapNonHap
BiGs = names(which(apply(meanExp[,3:7], 1, max) < meanExp[,2]))
BiGs = names(which(apply(meanExp[,1:2], 1, mean) > apply(meanExp[,3:7], 1, mean)))

#####BiGs is not selected ideally in any of these. When looking at meanExp[BiGs,] many genes look like they are haploid expressed, just later than the first set

plot(rowMeans(D2[dNdSgenes,MonoSamp]), rowMeans(D2[dNdSgenes,BiSamp]), pch=19, cex=3)
points(rowMeans(D2[BiGs,MonoSamp]), rowMeans(D2[BiGs,BiSamp]), pch=19, cex=3, col="green2")
points(rowMeans(D2[MonoGs,MonoSamp]), rowMeans(D2[MonoGs,BiSamp]), pch=19, cex=3, col="purple3")

plot(rowMeans(A2[IncGs,MonoSamp]), rowMeans(A2[IncGs,BiSamp]), pch=19, cex=3)
points(rowMeans(A2[BiGs,MonoSamp]), rowMeans(A2[BiGs,BiSamp]), pch=19, cex=3, col="green2")
points(rowMeans(A2[MonoGs,MonoSamp]), rowMeans(A2[MonoGs,BiSamp]), pch=19, cex=3, col="purple3")

boxplot(list(BiGs = RiceNS[BiGs], MonoGs = RiceNS[MonoGs]))
boxplot(list(BiGs = RiceNS[BiGs], MonoGs = RiceNS[MonoGs]), ylim=c(0,1))
boxplot(list(BiGs = RiceNS[names(which(RiceNS[BiGs] < 1))], MonoGs = RiceNS[names(which(RiceNS[MonoGs] < 1))])) 

plot(RiceNS[IncGs],apply(meanExp[IncGs,1:7], 1, max), xlim=c(0,1), pch=19, cex=3)


plot(apply(D2[IncGs,MonoSamp], 1, max), apply(D2[IncGs,BiSamp], 1, max), pch=19, cex=3)
points(apply(D2[c(MonoGs,BiGs),MonoSamp], 1, max), apply(D2[c(MonoGs,BiGs),BiSamp], 1, max), pch=19, cex=3, col="red2")

RiceNS3 = RiceNS2[MonoGs]
HapNonHap = list("Not or lowly expressed when haploid" = RiceNS[-which(names(RiceNS) %in% names(RiceNS3))], "Expressed when haploid" = RiceNS3)
boxplot(HapNonHap, ylim = c(0,1)) # Some outliers at really high values
HapNonHap = list("Not or lowly expressed when haploid" = RiceNS[names(which(RiceNS[-which(names(RiceNS) %in% names(RiceNS3))] < 1))], "Expressed when haploid" = RiceNS3[names(which(RiceNS3 < 1))])
boxplot(HapNonHap) # Some outliers at really high values


#length(RiceNS3) = 785
#length(RiceNS[-which(names(RiceNS) %in% names(RiceNS3))]) = 7947
HapNonHap2 = data.frame("NS" = c(AtNS[-which(names(AtNS) %in% rownames(A2d))], AtNS2), Class = c(rep(1, times = length(AtNS[-which(names(AtNS) %in% rownames(A2d))])), rep(2, times = length(AtNS2))))
pairwise.wilcox.test(HapNonHap2$NS, HapNonHap2$Class, p.adjust = 'holm')  # p-values
#data:  HapNonHap$NS and HapNonHap$Class
#  1
#2 <2e-16 #If this was done right, I think that means a very significant result



#Tomato CoGe
#Solanum lycopersicum (tomato) (v4, id57792): NCBI WindowMasker (Hard)
#Solanum tuberosum group Phureja DM1-3 516 R44 (potato) (v4.03, id52025): unmasked
#Capsicum sp. KIBR KC317 (cm334 Chromosome) (vv1.5, id25101): NCBI WindowMasker (Hard)
#Nicotiana tabacum Nitab-v4.5 (v1.0, id57796): NCBI WindowMasker (Hard)

#Move forward with one of the three options:
    TomNSdata = read.table('Tomato_Potato_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("[|].*", "", sub(".*LOC", "LOC", TomNSdata[,4])), sub("[|].*", '', sub(".*OL", "OL", TomNSdata[,8])), as.numeric(sub('.*[|]','',TomNSdata[,4])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','PotatoGene','Identity')

    TomNSdata = read.table('Tomato_Pepper_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("[|].*", "", sub(".*LOC", "LOC", TomNSdata[,4])), sub("[|].*", '', sub(".*OL", "OL", TomNSdata[,8])), as.numeric(sub('.*[|]','',TomNSdata[,4])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','PepperGene','Identity')

    TomNSdata = read.table('Tomato_Tobacco_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("[|].*", "", sub(".*LOC", "LOC", TomNSdata[,4])), sub("[|].*", '', sub(".*OL", "OL", TomNSdata[,8])), as.numeric(sub('.*[|]','',TomNSdata[,4])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','TobaccoGene','Identity')

## Remove putative orthologs that are not reciprocal best hits ##
removes = rep(FALSE, nrow(TomNSdata))
for (sgene in as.character(unique(TomNSdata[duplicated(TomNSdata[,13]),13]))) {
	tr = which(TomNSdata[,13] == sgene)
	removes[tr[rank(-TomNSdata[tr,15]) != 1]] = TRUE
}
for (mgene in as.character(unique(TomNSdata[duplicated(TomNSdata[,14]),14]))) {
	tr = which(TomNSdata[,14] == mgene)
	removes[tr[rank(-TomNSdata[tr,15]) != 1]] = TRUE
}
TomNSdata = TomNSdata[!removes,]
##

TomNSdata = TomNSdata[TomNSdata[,15] >= 80,]  # Require at least 80% identity between orthologs when calculating dNdS
TomNSdata = TomNSdata[-which(is.na(as.numeric(TomNSdata[,1]))),] # Remove genes that lack calculated synonymous and non-synonymous rates or are undefined
TomNSdata = TomNSdata[-which(as.numeric(TomNSdata[,1]) == 0),] # Remove genes with 0 synonymous subs
TomNS = as.numeric(TomNSdata[,2]) / as.numeric(TomNSdata[,1])
names(TomNS) = TomNSdata[,13]
