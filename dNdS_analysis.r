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
#AtNSdata = AtNSdata[-which(as.numeric(AtNSdata[,1]) == 0),] # Remove genes with 0 synonymous subs #skip for capsella and brassica
AtNS = as.numeric(AtNSdata[,2]) / as.numeric(AtNSdata[,1])
names(AtNS) = AtNSdata[,13]

load("HapExpScore10k.RData")
library('ComplexHeatmap')
library(circlize)

pseudocount = 1*10^6/quantile(colSums(D), p = .1)
A2 = sweep(D, 2, colSums(D), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(D[,colnames(A2b)] >= 10) >= 10, ]  # Require each gene to have at least 10 UMIs in at least 10 cells

MonoSamp = names(stages[grep("mono", stages)])
BiSamp = names(stages[-grep("mono", stages)])
plot(rowMeans(A3[,MonoSamp]), rowMeans(A3[,BiSamp]))
plot(rowMeans(D[dNdSgenes,MonoSamp]), rowMeans(D[dNdSgenes,BiSamp]))

meanExp = matrix(0, ncol = length(unique(stages)), nrow = nrow(A2))
rownames(meanExp) = rownames(A2)
colnames(meanExp) = c('tetrads', 'UM/BM_bi', 'UM/BM_mono', 'BM/Tri', 'Tri_bi', 'Tri_mono')
for (stage in colnames(meanExp)) {
	meanExp[,stage] = rowMeans(A2[,stages == stage])
}

library('dplyr')
AtNS0 = AtNS[which(names(AtNS) %in% rownames(meanExp))]
meanExp0 = meanExp[names(AtNS0),]
boxplot(list("ND" = AtNS0[names(which(meanExp0[,2] < 1))], "1-10" = AtNS0[which(between(meanExp0[,2], 1,10))], "11-29" = AtNS0[which(between(meanExp0[,2], 10,30))], "30-99" = AtNS0[which(between(meanExp0[,2], 30,100))], ">=100" = AtNS0[names(which(meanExp0[,2] >= 100))]), ylim=c(0,1), main="UM/BM biallelic")
boxplot(list("ND" = AtNS0[names(which(meanExp0[,3] < 1))], "1-10" = AtNS0[which(between(meanExp0[,3], 1,10))], "11-29" = AtNS0[which(between(meanExp0[,3], 10,30))], "30-99" = AtNS0[which(between(meanExp0[,3], 30,100))], ">=100" = AtNS0[names(which(meanExp0[,3] >= 100))]), ylim=c(0,1), main="UM/BM monoallelic")
boxplot(list("ND" = AtNS0[names(which(meanExp0[,5] < 1))], "1-10" = AtNS0[which(between(meanExp0[,5], 1,10))], "11-29" = AtNS0[which(between(meanExp0[,5], 10,30))], "30-99" = AtNS0[which(between(meanExp0[,5], 30,100))], ">=100" = AtNS0[names(which(meanExp0[,5] >= 100))]), ylim=c(0,1), main="Tri biallelic")
boxplot(list("ND" = AtNS0[names(which(meanExp0[,6] < 1))], "1-10" = AtNS0[which(between(meanExp0[,6], 1,10))], "11-29" = AtNS0[which(between(meanExp0[,6], 10,30))], "30-99" = AtNS0[which(between(meanExp0[,6], 30,100))], ">=100" = AtNS0[names(which(meanExp0[,6] >= 100))]), ylim=c(0,1), main="Tri monoallelic")


summary(apply(meanExp, 1, mean))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#    0.00     0.00     1.01    26.10    12.03 64989.87

dNdSgenes = names(AtNS[which(names(AtNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d
meanExp2 = meanExp[dNdSgenes,]
IncGs = rownames(meanExp2)[rowSums(meanExp2[,1:6] >= 100) > 0]
#meanExp2 = meanExp2[IncGs,]
#MonoGs = rownames(meanExp2)[rowSums(meanExp2[,c(3,4,6)] >= 100) > 0]
#BiGs = rownames(meanExp2)[rowSums(meanExp2[,c(3,4,6)] >= 100) == 0]

#dNdSComp = data.frame("NS" = c(AtNS, AtNS[IncGs], AtNS[BiGs], AtNS[MonoGs]), Class = c(rep("AllGs", times = length(AtNS)), rep("IncGs", times = length(IncGs)), rep("BiGs", times = length(BiGs)), rep("MonoGs", times = length(MonoGs))))
#boxplot(NS ~ Class, data=dNdSComp, ylim=c(0,1)) 
#table(dNdSComp$Class)
# AllGs   BiGs  IncGs MonoGs 
# 19641    154   1647   1493
#pairwise.wilcox.test(dNdSComp$NS, dNdSComp$Class, p.adjust = 'holm')  # p-values


#Trying using Y_samp (which is how haploid expressed a gene is for a sample, from ~0.5-1 typically) for setting MonoGs and BiGs
#Y_samp = ((ratio - .5)*Hap_all)[filtG,]+0.5 #This is how Y_samp was calculated in At_early_haploid_search
IncGs2 = IncGs[which(IncGs %in% rownames(Y_samp))]
#plot(rowMeans(Y_samp[IncGs2,MonoSamp], na.rm=T), rowMeans(Y_samp[IncGs2,BiSamp], na.rm=T))
summary(rowMeans(Y_samp[IncGs2,BiSamp], na.rm=T)) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.2250  0.5483  0.5917  0.5938  0.6297  1.0000       1
summary(rowMeans(Y_samp[IncGs2,MonoSamp], na.rm=T))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.2000  0.8283  0.8982  0.8340  0.9300  1.0000      25
MonoGs2 = names(which(rowMeans(Y_samp[IncGs2,MonoSamp], na.rm=T) > 0.6))
BiGs2 = names(which(rowMeans(Y_samp[IncGs2,], na.rm=T) <= 0.6))
#plot(AtNS[IncGs2], rowMeans(Y_samp[IncGs2,], na.rm=T))

dNdSComp2 = data.frame("dNdS" = c(AtNS, AtNS[IncGs2], AtNS[BiGs2], AtNS[MonoGs2]), Class = c(rep("AllGs", times = length(AtNS)), rep("IncGs", times = length(IncGs2)), rep("BiGs", times = length(BiGs2)), rep("MonoGs", times = length(MonoGs2))))
boxplot(dNdS ~ Class, data=dNdSComp2, ylim=c(0,1), main="0.6 haploid expression score threshold") 
table(dNdSComp2$Class)
# AllGs   BiGs  IncGs MonoGs #0.6
# 19641    114    587    489
# AllGs   BiGs  IncGs MonoGs #0.5
# 19641     12    587    522
pairwise.wilcox.test(dNdSComp2$dNdS, dNdSComp2$Class, p.adjust = 'holm')  # p-values






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
    RiceNSdata = cbind(RiceNSdata, sub("[.].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[|].*", '', sub(".*CDS:", "", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
    colnames(RiceNSdata)[13:15] = c('RiceGene','LeersiaGene','Identity')

    RiceNSdata = read.table('Rice_Brachypodium_distachyon_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    RiceNSdata = cbind(RiceNSdata, sub("[.].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[|].*", '', sub(".*CDS:", "", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
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
#RiceNSdata = RiceNSdata[-which(as.numeric(RiceNSdata[,1]) == 0),] # Remove genes with 0 synonymous subs #don't use for Lp or Bd
RiceNS = as.numeric(RiceNSdata[,2]) / as.numeric(RiceNSdata[,1])
names(RiceNS) = RiceNSdata[,13]


load("Rice_7-8_2025_with_meta.RData")
load("RiceRemap_7-8_2025.RData")
library('ComplexHeatmap')
library(circlize)
col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))
#Heatmap(HapExpScore, col = col_fun, cluster_columns=T, cluster_rows=T, show_row_names = FALSE, show_column_names = FALSE)

pseudocount = 1*10^6/quantile(colSums(D2), p = .1)
A2 = sweep(D2, 2, colSums(D2), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(D2[,colnames(A2b)] >= 5) >= 3, ]  # Require each gene to have at least 5 UMIs in at least 3 cells
over30k = names(which(colSums(D2) >25000))
over30k = over30k[which(Rice_meta[over30k,12] == "N")] #getting rid of the 2 no cell controls


#A3 = A2d[dNdSgenes,over30k]
#RiceNS2 = RiceNS[dNdSgenes]
#plot(rowMeans(A3),RiceNS2)
#plot(rowMeans(A3),RiceNS2, ylim=c(0,2))


BINR = function (xx, bin = 10^6) 
{
    bin = as.numeric(Rice_genes2[, 1]) * 10^6 + round(Rice_genes2[, 2]/bin)
    out = by(xx, bin, colSums)
    out2 = t(matrix(unlist(out), nrow = ncol(g1_2)))
    colnames(out2) = colnames(g1_2)
    rownames(out2) = names(out)
    return(out2)
}

library(ggplot2)
library(ggpubr)
g1_bin = BINR(g1_2)
g2_bin = BINR(g2_2)
AlleleFrac_bin = g1_bin/(g1_bin + g2_bin)
AlleleFrac_bin[(g1_bin+g2_bin) < 10] = NA #remove bins with <10 genoinformative transcripts
R_binUse = which(abs(rowMeans(AlleleFrac_bin, na.rm=T) - .5) < .4)  # Exclude bins with >90% of all transcripts mapping to the same allele across all samples
AlleleFrac_bin[-R_binUse,] = NA
FracMono_all = 100*colMeans(abs(AlleleFrac_bin - .5) >= .3, na.rm=T)

RiceStages = Rice_meta[over30k,9]
names(RiceStages) = over30k
RiceBiSamp = names(which(FracMono_all[over30k] < 20))
RiceMonoSamp = names(which(FracMono_all[over30k] > 20))
RiceStages[names(RiceStages[grep("tetrad", RiceStages)])] = "Tetrad"
RiceStages[names(which(RiceStages[RiceBiSamp] == "UM"))] = "UM_bi"
RiceStages[names(which(RiceStages[RiceMonoSamp] == "UM"))] = "UM_mono"

meanExp = matrix(0, ncol = length(unique(RiceStages)), nrow = nrow(A2))
rownames(meanExp) = rownames(A2)
colnames(meanExp) = c('Tetrad', 'UM_bi', 'UM_mono', 'UM/BM', 'BM', 'BM/Tri', 'Tri')
for (stage in colnames(meanExp)) {
	meanExp[,stage] = rowMeans(A2[,RiceStages == stage])
}

#summary(apply(meanExp, 1, mean))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#    0.000     0.000     0.083    17.921     3.786 30871.167


library('dplyr')
RiceNS0 = RiceNS[which(names(RiceNS) %in% rownames(meanExp))]
meanExp0 = meanExp[names(RiceNS0),]
boxplot(list("ND" = RiceNS0[names(which(meanExp0[,1] < 1))], "1-10" = RiceNS0[which(between(meanExp0[,1], 1,10))], "11-29" = RiceNS0[which(between(meanExp0[,1], 10,30))], "30-99" = RiceNS0[which(between(meanExp0[,1], 30,100))], ">=100" = RiceNS0[names(which(meanExp0[,1] >= 100))]), ylim=c(0,1), main="Tetrad")
boxplot(list("ND" = RiceNS0[names(which(meanExp0[,3] < 1))], "1-10" = RiceNS0[which(between(meanExp0[,3], 1,10))], "11-29" = RiceNS0[which(between(meanExp0[,3], 10,30))], "30-99" = RiceNS0[which(between(meanExp0[,3], 30,100))], ">=100" = RiceNS0[names(which(meanExp0[,3] >= 100))]), ylim=c(0,1), main="UM monoallelic")
boxplot(list("ND" = RiceNS0[names(which(meanExp0[,5] < 1))], "1-10" = RiceNS0[which(between(meanExp0[,5], 1,10))], "11-29" = RiceNS0[which(between(meanExp0[,5], 10,30))], "30-99" = RiceNS0[which(between(meanExp0[,5], 30,100))], ">=100" = RiceNS0[names(which(meanExp0[,5] >= 100))]), ylim=c(0,1), main="BM")
boxplot(list("ND" = RiceNS0[names(which(meanExp0[,7] < 1))], "1-10" = RiceNS0[which(between(meanExp0[,7], 1,10))], "11-29" = RiceNS0[which(between(meanExp0[,7], 10,30))], "30-99" = RiceNS0[which(between(meanExp0[,7], 30,100))], ">=100" = RiceNS0[names(which(meanExp0[,7] >= 100))]), ylim=c(0,1), main="Tri")


dNdSgenes = names(RiceNS[which(names(RiceNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d
meanExp2 = meanExp[dNdSgenes,]
IncGs = rownames(meanExp2)[rowSums(meanExp2[,1:7] >= 100) > 0]
meanExp2 = meanExp2[IncGs,]
MonoGs = rownames(meanExp2)[rowSums(meanExp2[,3:7] >= 100) > 0]
BiGs = rownames(meanExp2)[rowSums(meanExp2[,3:7] >= 100) == 0]
dNdSComp = data.frame("dNdS" = c(RiceNS, RiceNS[IncGs], RiceNS[BiGs], RiceNS[MonoGs]), 
    Class = c(rep("AllGs", times = length(RiceNS)), rep("IncGs", times = length(IncGs)), 
    rep("BiGs", times = length(BiGs)), rep("MonoGs", times = length(MonoGs))))
boxplot(dNdS ~ Class, data=dNdSComp, ylim=c(0,1), main="Rice and O. longistaminata") 
table(dNdSComp$Class)
# AllGs   BiGs  IncGs MonoGs #O. longistaminata
#  8732     42    653    611
# AllGs   BiGs  IncGs MonoGs #L. perrieri
# 16629     84   1235   1151
# AllGs   BiGs  IncGs MonoGs #B. distachyon
# 12019     69   1005    936
#pairwise.wilcox.test(dNdSComp$dNdS, dNdSComp$Class, p.adjust = 'holm')  # p-values


#Trying using how haploid expressed a gene is for a sample for setting MonoGs and BiGs. Didn't go well for rice, too few genes in BiGs category
RiceRatio = g1_2/(g1_2+g2_2)
RiceRatio[(g1_2+g2_2) < 5] = NA  # Remove measurements with under 5 genoinformative transcripts
RiceRatio = RiceRatio[,over30k] # Remove samples with <25,000 UMIs (over30k name is legacy, just hasn't been changed)
RiceRatio[which(abs(rowMeans(RiceRatio, na.rm=T)-0.5) >=0.3),] = NA #Remove genes with > 80% transcripts matching one allele

GeneHapScore = abs(RiceRatio - .5)+0.5 #This should go from 0.5-1
IncGs2 = IncGs[which(IncGs %in% rownames(GeneHapScore))]
summary(rowMeans(GeneHapScore[IncGs2,RiceBiSamp], na.rm=T))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.5000  0.6013  0.6547  0.6740  0.7000  1.0000    1623
summary(rowMeans(GeneHapScore[IncGs2,RiceMonoSamp], na.rm=T))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.5624  0.9203  0.9683  0.9245  0.9861  1.0000    1499 
MonoGs2 = names(which(rowMeans(GeneHapScore[IncGs2,RiceMonoSamp], na.rm=T) > 0.65))
BiGs2 = names(which(rowMeans(GeneHapScore[IncGs2,], na.rm=T) <= 0.65))
dNdSComp2 = data.frame("dNdS" = c(RiceNS, RiceNS[IncGs2], RiceNS[BiGs2], RiceNS[MonoGs2]), Class = c(rep("AllGs", times = length(RiceNS)), rep("IncGs", times = length(IncGs2)), rep("BiGs", times = length(BiGs2)), rep("MonoGs", times = length(MonoGs2))))
boxplot(dNdS ~ Class, data=dNdSComp2, ylim=c(0,1), main="0.65 haploid expression score threshold") 
table(dNdSComp2$Class)
# AllGs  IncGs MonoGs #0.5 no BiGs come through #O. longistaminata
#  8732    653    111
# AllGs   BiGs  IncGs MonoGs #0.65 #O. longistaminata
# 8732      6    653    107
#pairwise.wilcox.test(dNdSComp2$dNdS, dNdSComp2$Class, p.adjust = 'holm')  # p-values




#Tomato CoGe
#Solanum lycopersicum (tomato) (v4, id56999): unmasked #13 chr from ITAG (connected to Sol Genomics)
#Solanum tuberosum group Phureja DM1-3 516 R44 (potato) (v4.03, id52025): unmasked #Clonally propagated
#Capsicum annuum_cvCM334 (v1.0, id54650) #Diploid #Tomato_Pepper_CM334_CoGe.tab #
#Nicotiana tabacum Nitab-v4.5 (v1.0, id57796): NCBI WindowMasker (Hard) #Polyploid
#Capsicum annuum var. glabriusculum (v1.0, id54649) (wild chili pepper chiltepin) #Tomato_Pepper_glab_CoGe.tab #13 chr from Sol Genomics
#Solanum pimpinellifolium (v1.0, id54664): scaffolds #Tomato_pimpinellifolium_CoGe.tab #only 823 genes w/ dNdS data
#Solanum pennellii (v2.0, id28941) #Tomato_pennellii_CoGe.tab
#Solanum melongena (Japanese eggplant) (v1, id24647) #diploid #Tomato_Eggplant_CoGe.tab
#Nicotiana attenuata (vv2.4, id33947): unmasked #Diploid #Tomato_Nicotiana_attenuata_CoGe.tab #chr and scaffolds, 2311 Gs
#Nicotiana tomentosiformis (v1.0, id54668) #Diploid #Tomato_Nicotiana_tomentosiformis_CoGe.tab #all scaffolds, 1721 Gs
#Petunia inflata (v1.0, id54665) #Diploid #Tomato_Petunia_inflata_CoGe.tab
#Petunia axillaris subsp. axillaris (v1.0, id54663) #Diploid
#Coffea canephora (v1.0, id54651) #Diploid #has chr #Tomato_Coffea_canephora_CoGe.tab
#Coffea eugenioides (vCeug_1.0, id62727) #Diploid #not chr lv
#Cucumis melo var. cantalupo (v3.5.1, id52079) #Diploid #Tomato_Melon_CoGe.tab

#Move forward with one of the options:
    TomNSdata = read.table('Tomato_pimpinellifolium_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("\\..*", "", sub(".*CDS:Sol", "Sol", TomNSdata[,8])), sub("[|].*", '', sub(".*SPi", "SPi", TomNSdata[,4])), as.numeric(sub('.*[|]','',TomNSdata[,8])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','pimpiGene','Identity')

    TomNSdata = read.table('Tomato_pennellii_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("\\..*", "", sub(".*CDS:Sol", "Sol", TomNSdata[,8])), sub("[|].*", '', sub(".*SPi", "SPi", TomNSdata[,4])), as.numeric(sub('.*[|]','',TomNSdata[,8])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','pennelliiGene','Identity')

    TomNSdata = read.table('Tomato_Eggplant_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("\\..*", "", sub(".*CDS:Sol", "Sol", TomNSdata[,8])), sub("[|].*", '', sub(".*CDS:Sme", "Sme", TomNSdata[,4])), as.numeric(sub('.*[|]','',TomNSdata[,8])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','EggplantGene','Identity')

    TomNSdata = read.table('Tomato_Pepper_CM334_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("\\..*", "", sub(".*CDS:Sol", "Sol", TomNSdata[,8])), sub("[|].*", '', sub(".*CDS:CA", "CA", TomNSdata[,4])), as.numeric(sub('.*[|]','',TomNSdata[,8])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','PepperCM334Gene','Identity')

    TomNSdata = read.table('Tomato_Pepper_glab_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("\\..*", "", sub(".*CDS:Sol", "Sol", TomNSdata[,8])), sub("[|].*", '', sub(".*Ca", "Ca", TomNSdata[,4])), as.numeric(sub('.*[|]','',TomNSdata[,8])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','PepperChiltepinGene','Identity')

    TomNSdata = read.table('Tomato_Nicotiana_attenuata_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("\\..*", "", sub(".*CDS:Sol", "Sol", TomNSdata[,8])), sub("[|].*", '', sub(".*NIAT", "NIAT", TomNSdata[,4])), as.numeric(sub('.*[|]','',TomNSdata[,8])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','N.attenuataGene','Identity')

    TomNSdata = read.table('Tomato_Nicotiana_tomentosiformis_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("\\..*", "", sub(".*CDS:Sol", "Sol", TomNSdata[,8])), sub("[|].*", '', sub(".*cds", "cds", TomNSdata[,4])), as.numeric(sub('.*[|]','',TomNSdata[,8])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','N.tomentoGene','Identity')

    TomNSdata = read.table('Tomato_Petunia_inflata_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("\\..*", "", sub(".*CDS:Sol", "Sol", TomNSdata[,8])), sub("[|].*", '', sub(".*Peinf", "Peinf", TomNSdata[,4])), as.numeric(sub('.*[|]','',TomNSdata[,8])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','PetuniaGene','Identity')

    TomNSdata = read.table('Tomato_Coffea_canephora_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    TomNSdata = cbind(TomNSdata, sub("\\..*", "", sub(".*CDS:Sol", "Sol", TomNSdata[,8])), sub("[|].*", '', sub(".*GSCOCT", "GSCOCT", TomNSdata[,4])), as.numeric(sub('.*[|]','',TomNSdata[,8])))
    colnames(TomNSdata)[13:15] = c('TomatoGene','CoffeaGene','Identity')


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
#TomNSdata = TomNSdata[-which(as.numeric(TomNSdata[,1]) == 0),] # Remove genes with 0 synonymous subs 
TomNS = as.numeric(TomNSdata[,2]) / as.numeric(TomNSdata[,1])
names(TomNS) = TomNSdata[,13]

load("TomatoRemap_11_2024-7_2025_with_meta.RData")

pseudocount = 1*10^6/quantile(colSums(D2), p = .1)
A2 = sweep(D2, 2, colSums(D2), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(D2[,colnames(A2b)] >= 5) >= 3, ]  # Require each gene to have at least 5 UMIs in at least 3 cells
over30k = names(which(colSums(D2) >25000))
over30k = over30k[which(Sl_meta[over30k,14] == "N")] #getting rid of no cell controls


BINSl = function (xx, bin = 10^6) 
{
    bin = as.numeric(Sl_genes2[, 1]) * 10^6 + round(Sl_genes2[, 2]/bin)
    out = by(xx, bin, colSums)
    out2 = t(matrix(unlist(out), nrow = ncol(g1_2)))
    colnames(out2) = colnames(g1_2)
    rownames(out2) = names(out)
    return(out2)
}

library(ggplot2)
library(ggpubr)
g1_bin = BINSl(g1_2)
g2_bin = BINSl(g2_2)
AlleleFrac_bin = g1_bin/(g1_bin + g2_bin)
AlleleFrac_bin[(g1_bin+g2_bin) < 10] = NA #remove bins with <10 genoinformative transcripts
Sl_binUse = which(abs(rowMeans(AlleleFrac_bin, na.rm=T) - .5) < .4)  # Exclude bins with >90% of all transcripts mapping to the same allele across all samples
AlleleFrac_bin[-Sl_binUse,] = NA
FracMono_all = 100*colMeans(abs(AlleleFrac_bin - .5) >= .3, na.rm=T)
FracMono_all[names(which(colSums(!is.na(AlleleFrac_bin))<10))] = NA

TomStages = Sl_meta[over30k,8]
#TomStages = sub('tetrad_sep', 'tetrad', TomStages)
names(TomStages) = over30k
#TomBiSamp = names(which(FracMono_all[over30k] < 20))
#TomMonoSamp = names(which(FracMono_all[over30k] > 20))
#TomStages[names(which(TomStages[TomBiSamp] == "UM"))] = "UM_bi"
#TomStages[names(which(TomStages[TomMonoSamp] == "UM"))] = "UM_mono"


meanExp = matrix(0, ncol = length(unique(TomStages)), nrow = nrow(A2))
rownames(meanExp) = rownames(A2)
colnames(meanExp) = c('Meiosis', 'tetrad', 'UM', 'UM/BM', 'BM', 'BM_shed')
#colnames(meanExp) = c('Meiosis', 'tetrad', 'UM_bi', 'UM_mono', 'UM/BM', 'BM', 'BM_shed')
for (stage in colnames(meanExp)) {
	meanExp[,stage] = rowMeans(A2[,TomStages == stage])
}

#summary(apply(meanExp, 1, mean))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#     0.00      0.03      1.42     29.80     10.39 232873.56

#library('dplyr')
TomNS0 = TomNS[which(names(TomNS) %in% rownames(meanExp))]
meanExp0 = meanExp[names(TomNS0),]
boxplot(list("ND" = TomNS0[names(which(meanExp0[,1] < 1))], "1-10" = TomNS0[which(between(meanExp0[,1], 1,10))], "11-29" = TomNS0[which(between(meanExp0[,1], 10,30))], "30-99" = TomNS0[which(between(meanExp0[,1], 30,100))], ">=100" = TomNS0[names(which(meanExp0[,1] >= 100))]), ylim=c(0,1), main="Meiosis")
boxplot(list("ND" = TomNS0[names(which(meanExp0[,2] < 1))], "1-10" = TomNS0[which(between(meanExp0[,2], 1,10))], "11-29" = TomNS0[which(between(meanExp0[,2], 10,30))], "30-99" = TomNS0[which(between(meanExp0[,2], 30,100))], ">=100" = TomNS0[names(which(meanExp0[,2] >= 100))]), ylim=c(0,1), main="Tetrad")
boxplot(list("ND" = TomNS0[names(which(meanExp0[,3] < 1))], "1-10" = TomNS0[which(between(meanExp0[,3], 1,10))], "11-29" = TomNS0[which(between(meanExp0[,3], 10,30))], "30-99" = TomNS0[which(between(meanExp0[,3], 30,100))], ">=100" = TomNS0[names(which(meanExp0[,3] >= 100))]), ylim=c(0,1), main="UM")
boxplot(list("ND" = TomNS0[names(which(meanExp0[,4] < 1))], "1-10" = TomNS0[which(between(meanExp0[,4], 1,10))], "11-29" = TomNS0[which(between(meanExp0[,4], 10,30))], "30-99" = TomNS0[which(between(meanExp0[,4], 30,100))], ">=100" = TomNS0[names(which(meanExp0[,4] >= 100))]), ylim=c(0,1), main="UM/BM")
boxplot(list("ND" = TomNS0[names(which(meanExp0[,5] < 1))], "1-10" = TomNS0[which(between(meanExp0[,5], 1,10))], "11-29" = TomNS0[which(between(meanExp0[,5], 10,30))], "30-99" = TomNS0[which(between(meanExp0[,5], 30,100))], ">=100" = TomNS0[names(which(meanExp0[,5] >= 100))]), ylim=c(0,1), main="BM")
boxplot(list("ND" = TomNS0[names(which(meanExp0[,6] < 1))], "1-10" = TomNS0[which(between(meanExp0[,6], 1,10))], "11-29" = TomNS0[which(between(meanExp0[,6], 10,30))], "30-99" = TomNS0[which(between(meanExp0[,6], 30,100))], ">=100" = TomNS0[names(which(meanExp0[,6] >= 100))]), ylim=c(0,1), main="BM_shed")
#boxplot(list("ND" = TomNS0[names(which(meanExp0[,7] < 1))], "1-10" = TomNS0[which(between(meanExp0[,7], 1,10))], "11-29" = TomNS0[which(between(meanExp0[,7], 10,30))], "30-99" = TomNS0[which(between(meanExp0[,7], 30,100))], ">=100" = TomNS0[names(which(meanExp0[,7] >= 100))]), ylim=c(0,1), main="BM_shed")



IncGs = rownames(meanExp)[rowSums(meanExp[,2:6] >= 100) > 0]
meanExp2 = meanExp[IncGs,]
MonoGs = rownames(meanExp2)[rowSums(meanExp2[,3:6] >= 100) > 0]
BiGs = rownames(meanExp2)[rowSums(meanExp2[,3:6] >= 100) == 0] 
dNdSgenes = names(TomNS[which(names(TomNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d
IncGs = IncGs[which(IncGs %in% dNdSgenes)] 
BiGs = BiGs[which(BiGs %in% dNdSgenes)] 
MonoGs = MonoGs[which(MonoGs %in% dNdSgenes)] 


CM334_TomNS = TomNS #6231 genes w/ dNdS data
CM334_BiGs = BiGs #52 BiGs w/ dNdS data         #227 total      #22.9% w/ dNdS data
CM334_IncGs = IncGs #543 IncGs w/ dNdS data     #2094 total     #25.9% w/ dNdS data
CM334_MonoGs = MonoGs #491 MonoGs w/ dNdS data  #1867 total     #26.3% w/ dNdS data

Chiltepin_TomNS = TomNS #6704 genes w/ dNdS data
Chiltepin_BiGs = BiGs #61 BiGs w/ dNdS data         #227 total      #26.9% w/ dNdS data
Chiltepin_IncGs = IncGs #544 IncGs w/ dNdS data     #2094 total     #26.0% w/ dNdS data
Chiltepin_MonoGs = MonoGs #483 MonoGs w/ dNdS data  #1867 total     #25.9% w/ dNdS data


Both_TomNS = CM334_TomNS[which(CM334_TomNS %in% Chiltepin_TomNS)] #529
Both_IncGs = CM334_IncGs[which(CM334_IncGs %in% Chiltepin_IncGs)] #284
Both_MonoGs = CM334_MonoGs[which(CM334_MonoGs %in% Chiltepin_MonoGs)] #256
Both_BiGs = CM334_BiGs[which(CM334_BiGs %in% Chiltepin_BiGs)] #28

uniqueCM334_BiGs = CM334_BiGs[-which(CM334_BiGs %in% Both_BiGs)]
uniqueChiltepin_BiGs = Chiltepin_BiGs[-which(Chiltepin_BiGs %in% Both_BiGs)]

dNdSComp = data.frame("dNdS" = c(CM334_TomNS[Both_BiGs], CM334_TomNS[uniqueCM334_BiGs], 
    Chiltepin_TomNS[Both_BiGs], Chiltepin_TomNS[uniqueChiltepin_BiGs]), 
    Class = c(rep("CM334 Matching BiGs", times = length(Both_BiGs)), rep("CM334 Unique BiGs", times = length(uniqueCM334_BiGs)), 
    rep("Chiltepin Matching BiGs", times = length(Both_BiGs)), rep("Chiltepin Unique BiGs", times = length(uniqueChiltepin_BiGs))))
boxplot(dNdS ~ Class, data=dNdSComp, ylim=c(0,1), main="Tomato vs. two peppers") 
table(dNdSComp$Class)






dNdSgenes = names(TomNS[which(names(TomNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d
meanExp2 = meanExp[dNdSgenes,]
IncGs = rownames(meanExp2)[rowSums(meanExp2[,2:7] >= 100) > 0]
meanExp2 = meanExp2[IncGs,]
MonoGs = rownames(meanExp2)[rowSums(meanExp2[,4:7] >= 100) > 0]
BiGs = rownames(meanExp2)[rowSums(meanExp2[,4:7] >= 100) == 0]

dNdSComp = data.frame("dNdS" = c(TomNS, TomNS[IncGs], TomNS[BiGs], TomNS[MonoGs]), 
    Class = c(rep("AllGs", times = length(TomNS)), rep("IncGs", times = length(IncGs)), 
    rep("BiGs", times = length(BiGs)), rep("MonoGs", times = length(MonoGs))))
#boxplot(dNdS ~ Class, data=dNdSComp, ylim=c(0,1), main="Tomato and Eggplant") 
table(dNdSComp$Class)
#pairwise.wilcox.test(dNdSComp$dNdS, dNdSComp$Class, p.adjust = 'holm')  # p-values


Egg_TomNS = TomNS #4652
Egg_BiGs = BiGs #133
Egg_IncGs = IncGs #424
Egg_MonoGs = MonoGs #291

CM334_TomNS = TomNS #6231
CM334_BiGs = BiGs #179
CM334_IncGs = IncGs #579
CM334_MonoGs = MonoGs #400

Chiltepin_TomNS = TomNS #6704
Chiltepin_BiGs = BiGs #224
Chiltepin_IncGs = IncGs #625
Chiltepin_MonoGs = MonoGs #401


Both_TomNS = CM334_TomNS[which(CM334_TomNS %in% Chiltepin_TomNS)] #529
Both_IncGs = CM334_IncGs[which(CM334_IncGs %in% Chiltepin_IncGs)] #308
Both_MonoGs = CM334_MonoGs[which(CM334_MonoGs %in% Chiltepin_MonoGs)] #204
Both_BiGs = CM334_BiGs[which(CM334_BiGs %in% Chiltepin_BiGs)] #104

summary(Egg_TomNS[Egg_BiGs])
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
 0.00918  0.02582  0.10705  1.62293  0.20213 99.05589


summary(CM334_TomNS[CM334_BiGs])
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  0.0000   0.0181   0.1081   9.8581   0.2652 116.0000
summary(Chiltepin_TomNS[Chiltepin_BiGs])
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
 0.00087  0.02243  0.10629  0.86084  0.24567 98.99495

summary(CM334_TomNS[Both_BiGs])
    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
 0.00000   0.03488   0.12893   5.66703   0.23709 105.00000
summary(Chiltepin_TomNS[Both_BiGs])
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
 0.00993  0.03829  0.12740  1.17071  0.22737 98.99495

summary(CM334_TomNS[CM334_BiGs]>1)
   Mode   FALSE    TRUE 
logical     159      20 
summary(Chiltepin_TomNS[Chiltepin_BiGs]>1)
   Mode   FALSE    TRUE 
logical     215       9

summary(CM334_TomNS[Both_BiGs]>1) 
   Mode   FALSE    TRUE 
logical      97       7
summary(Chiltepin_TomNS[Both_BiGs]>1)       
   Mode   FALSE    TRUE 
logical     101       3

summary(CM334_TomNS[CM334_MonoGs]>1)
   Mode   FALSE    TRUE 
logical     375      25
summary(Chiltepin_TomNS[Chiltepin_MonoGs]>1)
   Mode   FALSE    TRUE 
logical     392       9

summary(CM334_TomNS[Both_MonoGs]>1) 
   Mode   FALSE    TRUE 
logical     200       4
summary(Chiltepin_TomNS[Both_MonoGs]>1)
   Mode   FALSE    TRUE 
logical     203       1



#Trying using how haploid expressed a gene is for a sample for setting MonoGs and BiGs. Didn't go well for rice, too few genes in BiGs category
TomRatio = g1_2/(g1_2+g2_2)
TomRatio[(g1_2+g2_2) < 5] = NA  # Remove measurements with under 5 genoinformative transcripts
TomRatio = TomRatio[,over30k] # Remove samples with <10,000 UMIs (over30k name is legacy, just hasn't been changed)
TomRatio[which(abs(rowMeans(TomRatio, na.rm=T)-0.5) >=0.3),] = NA #Remove genes with > 80% transcripts matching one allele
GeneHapScore = abs(TomRatio - .5)+0.5 #This should go from 0.5-1

dNdSgenes = names(TomNS[which(names(TomNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d
meanExp2 = meanExp[dNdSgenes,]
IncGs = rownames(meanExp2)[rowSums(meanExp2[,1:7] >= 100) > 0]
meanExp2 = meanExp2[IncGs,]
IncGs2 = IncGs[which(IncGs %in% rownames(GeneHapScore))]
summary(rowMeans(GeneHapScore[IncGs2,TomBiSamp], na.rm=T))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.5000  0.6027  0.6530  0.6528  0.6932  0.8333     176
summary(rowMeans(GeneHapScore[IncGs2,TomMonoSamp], na.rm=T))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.6679  0.9708  0.9868  0.9730  0.9974  1.0000     113 
MonoGs2 = names(which(rowMeans(GeneHapScore[IncGs2,TomMonoSamp], na.rm=T) > 0.65))
BiGs2 = IncGs2[-which(IncGs2 %in% MonoGs2)]
#BiGs2 = unique(c(BiGs2, names(which(rowSums(Sl_tissues2 >= 100) > 0))))
#BiGs2 = BiGs2[-which(BiGs2 %in% MonoGs2)]
#IncGs2 = c(BiGs2, MonoGs2)
#BiGs2 = names(which(rowMeans(GeneHapScore[IncGs2,], na.rm=T) <= 0.65))
dNdSComp2 = data.frame("dNdS" = c(TomNS, TomNS[IncGs2], TomNS[BiGs2], TomNS[MonoGs2]), Class = c(rep("AllGs", times = length(TomNS)), rep("IncGs", times = length(IncGs2)), rep("BiGs", times = length(BiGs2)), rep("MonoGs", times = length(MonoGs2))))
boxplot(dNdS ~ Class, data=dNdSComp2, ylim=c(0,1), main="0.65 haploid expression score threshold") 
table(dNdSComp2$Class)
# AllGs   BiGs  IncGs MonoGs #0.65 #Potato
#  8971     13    738    281
# AllGs   BiGs  IncGs MonoGs #0.65 #Pepper CM334
#  6231     10    593    245
# AllGs   BiGs  IncGs MonoGs #0.65 #Tobacco
# 5710     10    566    216
#pairwise.wilcox.test(dNdSComp2$dNdS, dNdSComp2$Class, p.adjust = 'holm')  # p-values


load("Sl_tissues.RData")
dNdSgenes = names(TomNS[which(names(TomNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d
meanExp2 = meanExp[dNdSgenes,]
GametoGs = names(which(rowSums(meanExp2[,2:7] >= 100) > 0)) 
dNdSgenes_tissues = names(TomNS[which(names(TomNS) %in% rownames(Sl_tissues_TPM))])
Sl_tissues2 = Sl_tissues_TPM[dNdSgenes_tissues,c(1:3,8)]
SporoGs = unique(c(names(which(meanExp2[,1] >= 100)), names(which(rowSums(Sl_tissues2 >= 100) > 0))))
#SporoGs = names(which(rowSums(Sl_tissues2 >= 100) > 0))
IncGs = SporoGs[which(SporoGs %in% GametoGs)]
MonoGs = names(which(rowSums(meanExp2[IncGs,4:7] >= 100) > 0))
BiGs = IncGs[-which(IncGs %in% MonoGs)]
dNdSComp = data.frame("dNdS" = c(TomNS, TomNS[IncGs], TomNS[BiGs], TomNS[MonoGs]), 
    Class = c(rep("AllGs", times = length(TomNS)), rep("IncGs", times = length(IncGs)), 
    rep("BiGs", times = length(BiGs)), rep("MonoGs", times = length(MonoGs))))
boxplot(dNdS ~ Class, data=dNdSComp, ylim=c(0,1), main="Tomato and Eggplant") 
table(dNdSComp$Class)
# AllGs   BiGs  IncGs MonoGs #Potato
#  8971    113    488    375
# AllGs   BiGs  IncGs MonoGs #Pepper CM334
#  6231     15    363    348
#pairwise.wilcox.test(dNdSComp$dNdS, dNdSComp$Class, p.adjust = 'holm')  # p-values




dNdSgenes = names(TomNS[which(names(TomNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d
meanExp2 = meanExp[dNdSgenes,]
IncGs = rownames(meanExp2)[rowSums(meanExp2[,1:7] >= 100) > 0]
meanExp2 = meanExp2[IncGs,]
MonoGs = rownames(meanExp2)[rowSums(meanExp2[,4:7] >= 100) > 0]
BiGs = rownames(meanExp2)[rowSums(meanExp2[,4:7] >= 100) == 0]

dNdSgenes_tissues = names(TomNS[which(names(TomNS) %in% rownames(Sl_tissues_TPM))])
Sl_tissues2 = Sl_tissues_TPM[dNdSgenes_tissues,c(1:3,8)]
BiGs = unique(c(names(which(rowSums(meanExp2[,4:7] >= 100) == 0)), names(which(rowSums(Sl_tissues2 >= 100) > 0))))
#BiGs = names(which(rowSums(Sl_tissues2 >= 100) > 0))
BiGs = BiGs[-which(BiGs %in% MonoGs)]
IncGs = c(MonoGs, BiGs)
dNdSComp = data.frame("dNdS" = c(TomNS, TomNS[IncGs], TomNS[BiGs], TomNS[MonoGs]), 
    Class = c(rep("AllGs", times = length(TomNS)), rep("IncGs", times = length(IncGs)), 
    rep("BiGs", times = length(BiGs)), rep("MonoGs", times = length(MonoGs))))
boxplot(dNdS ~ Class, data=dNdSComp, ylim=c(0,1), main="Tomato and Pepper") 
table(dNdSComp$Class)

#pairwise.wilcox.test(dNdSComp$dNdS, dNdSComp$Class, p.adjust = 'holm')  # p-values



dNdSgenes = names(TomNS[which(names(TomNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d
meanExp2 = meanExp[dNdSgenes,]
IncGs = rownames(meanExp2)[rowSums(meanExp2[,1:7] >= 100) > 0]
meanExp2 = meanExp2[IncGs,]
MonoGs = rownames(meanExp2)[rowSums(meanExp2[,4:7] >= 100) > 0]
BiGs = rownames(meanExp2)[rowSums(meanExp2[,4:7] >= 100) == 0]

dNdSgenes = names(TomNS[which(names(TomNS) %in% rownames(A2d))]) # Genes with dn/ds data that are in A2d
meanExp2 = meanExp[dNdSgenes,]
GametoGs = names(which(rowSums(meanExp2[,2:7] >= 100) > 0)) 
dNdSgenes_tissues = names(TomNS[which(names(TomNS) %in% rownames(Sl_tissues_TPM))])
Sl_tissues2 = Sl_tissues_TPM[dNdSgenes_tissues,c(1:3,8)]
SporoGs = unique(c(names(which(meanExp2[,1] >= 100)), names(which(rowSums(Sl_tissues2 >= 100) > 0))))
#SporoGs = names(which(rowSums(Sl_tissues2 >= 100) > 0))
IncGs = SporoGs[which(SporoGs %in% GametoGs)]
MonoGs = names(which(rowSums(meanExp2[IncGs,4:7] >= 100) > 0))
BiGs = IncGs[-which(IncGs %in% MonoGs)]
dNdSComp = data.frame("dNdS" = c(TomNS, TomNS[IncGs], TomNS[BiGs], TomNS[MonoGs]), 
    Class = c(rep("AllGs", times = length(TomNS)), rep("IncGs", times = length(IncGs)), 
    rep("BiGs", times = length(BiGs)), rep("MonoGs", times = length(MonoGs))))
boxplot(dNdS ~ Class, data=dNdSComp, ylim=c(0,1), main="Tomato and Eggplant") 
table(dNdSComp$Class)
