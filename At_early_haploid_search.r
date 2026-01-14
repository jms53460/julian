#Trying to look for early haploid genes
setwd('C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R')
load("At_data_10_2025.RData")

library(readxl)
At_Stages <- read_excel("C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R/At_Stages2.xlsx")

At_meta <- read_excel("C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R/At_meta2.xlsx")
At_meta = as.data.frame(At_meta)
rownames(At_meta) = At_meta$Sample
At_meta[grep("tetrad", At_meta[,9]),9] = "tetrads"
At_meta[which(At_meta[,8] == "35"),8] = "0.35"

calcFracs = function(stage) {
	F_Col0 = ratio[,stages == stage]

	frac = abs(F_Col0 - 0.5) + .5

	keeps = rowSums(!is.na(frac)) >= 5  # Only consider genes with at least 10 genoinformative transcripts in at least 5 cells 
	frac = rowMeans(frac[keeps, ], na.rm=T)
	meanBias = rowMeans(F_Col0, na.rm=T)[names(frac)]
	data.frame(gene=names(frac),bias=meanBias,frac=frac,stage=stage)
}


over30k = names(which(colSums(D) >25000))

stages = as.character(At_meta[over30k,9])
stages = factor(stages, levels = c('tetrads','UM','UM/BM','BM','BM/Tri','Tri'))
names(stages) = over30k

library(dplyr)

ratio = g1/(g1+g2)
ratio[(g1+g2) < 5] = NA  # Remove measurements with under 5 genoinformative transcripts
ratio = ratio[,over30k]
ratio[which(abs(rowMeans(ratio, na.rm=T)-0.5) >=0.4),] = NA

#calcFracs2 = function(stage) {
	
#    F_Col0 = ratio[,stages == stage]
#    frac = HapFrac_subset

#	keeps = rowSums(!is.na(frac)) >= 10  # Only consider genes with at least 10 genoinformative transcripts in at least 10 cells 
#	frac = rowMeans(frac[keeps, ], na.rm=T)
#	meanBias = rowMeans(F_Col0, na.rm=T)[names(frac)]
#	data.frame(gene=names(frac),bias=meanBias,frac=frac,stage=stage)
#}

calcFracs2 = function(stage){
    stages = stages[lowHapFrac]
    stages[which(stages == "UM" | stages == "BM")] = "UM/BM"
    filtG = names(which(rowSums(!is.na(ratio)) >= 10))
    filtS = names(stages[stages == stage])
    X = rowMeans(ratio[filtG, filtS], na.rm=T)
    X = X[which(!is.na(X))]
    Y = rowMeans(((ratio - .5)*Hap_all)[filtG, filtS], na.rm=T)+0.5
    Y = Y[names(X)]
    Y = Y[which(!is.na(Y))]
    X = X[names(Y)]
    data.frame(gene=names(X), bias=X, frac=Y, stage=stage)
}

MonoCols = c('#D4D4D4', '#BCBDDC', '#FF9955', '#D45500')
ggFnFracs = function(stageS) {
	stage = fracs$stage == stageS
	ggscatter(fracs[stage,], x = 'bias', y = 'frac', color = 'class', size = 1) + scale_color_manual(values = c('#d4d4d4', MonoCols[-1])) + theme(legend.position = 'none', plot.title = element_text(hjust = 0.5), axis.title = element_blank()) + ylim(c(0,100)) + scale_x_continuous(breaks = seq(0,100,10), labels=c('0', rep('',4), '50', rep('',4), '100'), limits = c(0,100)) + ggtitle(stageS)
}

#lowHapFrac = which(HapFrac_all > 0) 
#lowHapFrac = which(HapFrac_all < 0.6)
lowHapFrac = which(between(HapFrac_all, 0.6, 0.8))

fracs = Reduce(function(...) merge(...,all=T), lapply(c("UM/BM", "Tri"), calcFracs2))
#fracs = Reduce(function(...) merge(...,all=T), lapply(c("UM","UM/BM", "BM", "Tri"), calcFracs2))
#fracs = Reduce(function(...) merge(...,all=T), lapply(c("BM", "Tri"), calcFracs2))

#fracs = fracs[which(!is.na(fracs$bias)),]

fracs = fracs[order(fracs$stage != 'Tri', fracs$stage),]
fracs$class = factor(c('<.8','<.8','exc','>.8','exc','>.95'), levels = c('exc','<.8','>.8','>.95'))[2*(fracs$frac > .95) + 2*(fracs$frac > .8) + 1 + (abs(fracs$bias - .5) < (fracs$frac-.5)*.9)]
fracs$frac = fracs$frac*100
fracs$bias = fracs$bias*100

library(ggplot2)
library(ggpubr)


svg('At_volcanos2.svg', width=4, height=2)
ggarrange(ggFnFracs('UM/BM'), ggFnFracs('Tri'), ncol = 2, align = 'h')
dev.off()


svg('At_volcanos.svg', width=6, height=2)
ggarrange(ggFnFracs('UM'), ggFnFracs('UM/BM'), ggFnFracs('BM'), ggFnFracs('Tri'), ncol = 4, align = 'h')
dev.off()


svg('At_volcano_UM_BM.svg', width=2, height=2)
ggFnFracs('UM/BM')
dev.off()


fracByStage = matrix(data=NA, nrow=length(rownames(D)), ncol=4)
rownames(fracByStage) = rownames(D)
colnames(fracByStage) = c("UM", "UM/BM", "BM", "Tri")
for (stage in c("UM", "UM/BM", "BM", "Tri")){
    stageGenes = fracs[fracs$stage == stage,1]
    stageFracs = fracs[fracs$stage == stage,3]
    stageFracs[which(fracs[fracs$stage == stage,5] == "exc")] = NA
    fracByStage[stageGenes,stage] = stageFracs
}

pairs(fracByStage)
#pairs function, can simplify with unique and round


####Trying bootstrap
#lowHapFrac = which(HapFrac_all > 0) 
#lowHapFrac = which(HapFrac_all < 0.6)
lowHapFrac = which(between(HapFrac_all, 0.6, 0.8))

over30k = names(which(colSums(D) >25000))
stages = as.character(At_meta[over30k,9])
stages = factor(stages, levels = c('tetrads','UM','UM/BM','BM','BM/Tri','Tri'))
names(stages) = over30k

ratio = g1/(g1+g2)
ratio[(g1+g2) < 5] = NA  # Remove measurements with under 5 genoinformative transcripts
ratio = ratio[,over30k]
ratio[which(abs(rowMeans(ratio, na.rm=T)-0.5) >=0.4),] = NA #Remove genes with > 90% transcripts matching one allele

BootSamps = stages[lowHapFrac]
BootSamps[which(BootSamps == "UM" | BootSamps == "BM")] = "UM/BM"
BootSamps = names(BootSamps[BootSamps == "UM/BM"])
filtG = names(which(rowSums(!is.na(ratio[,BootSamps])) >= 5))
filtG = names(which(rowSums(!is.na(Hap_all[filtG,BootSamps])) >= 5))


Boot = function(x){
    samps = sample(BootSamps, replace=TRUE)
    X = rowMeans(ratio[filtG, samps], na.rm=T)
    Y = rowMeans(((ratio - .5)*Hap_all)[filtG, samps], na.rm=T)+0.5
   cbind(X,Y)
}
set.seed(1)
Boots = sapply(1:100, Boot)
Boots = lapply(list(X = Boots[1:length(filtG),], Y = Boots[(1:length(filtG)) + length(filtG),]), function(xx) { 
        rownames(xx) = filtG
        return(xx)
    })

HapExpScore = (Boots$Y - abs(Boots$X - 0.5) - 0.5) / (1 - abs(Boots$X - 0.5) - 0.5)
summary(rowMeans(HapExpScore, na.rm=T))
library('ComplexHeatmap')
library(circlize)
col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))
Heatmap(HapExpScore, col = col_fun, cluster_columns=T, cluster_rows=T, show_row_names = FALSE, show_column_names = FALSE)
MeanHapExpScore = rowMeans(HapExpScore, na.rm=T)



#HapExpScore per gene per sample
over30k = names(which(colSums(D) >25000))
stages = as.character(At_meta[over30k,9])
stages = factor(stages, levels = c('tetrads','UM','UM/BM','BM','BM/Tri','Tri'))
names(stages) = over30k
NonPollen = names(stages[which(stages == "tetrads" | stages == "UM" | stages == "UM/BM" | stages == "BM")])

ratio = g1/(g1+g2)
ratio[(g1+g2) < 5] = NA  # Remove measurements with under 5 genoinformative transcripts
ratio = ratio[,over30k]
#ratio[which(abs(rowMeans(ratio, na.rm=T)-0.5) >=0.3),] = NA #Remove genes with > 80% transcripts matching one allele

filtG = names(which(rowSums(!is.na(ratio[,NonPollen])) >= 5))
filtG = names(which(rowSums(!is.na(Hap_all[filtG,NonPollen])) >= 5))

X_samp = ratio[filtG,]
Y_samp = ((ratio - .5)*Hap_all)[filtG,]+0.5
Ranked500G = rownames(X_samp[order(rowSums(is.na(X_samp))),])[1:500]
#Y[,1:length(over30k)] = rowMeans(((ratio - .5)*Hap_all)[filtG,], na.rm=T)+0.5
#HapExpScore = (Y - abs(X - 0.5)) / (1 - abs(X - 0.5))
#summary(rowMeans(HapExpScore, na.rm=T))
#library('ComplexHeatmap')
#library(circlize)
col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))
Heatmap(X_samp[Ranked500G,], name = '% transcipts from Col-0 allele', col = col_fun, cluster_columns=F, cluster_rows=F, show_row_names = FALSE, show_column_names = FALSE, use_raster=F)
Heatmap(Y_samp[Ranked500G,], name = '% transcripts matching haplotype', col = col_fun, cluster_columns=F, cluster_rows=F, show_row_names = FALSE, show_column_names = FALSE, use_raster=F)


BIN2 = function (xx, bin = 10^6) 
{
    bin = as.numeric(genes[, 1]) * 10^6 + round(genes[, 2]/bin)
    out = by(xx, bin, colSums)
    out2 = t(matrix(unlist(out), nrow = ncol(g1)))
    colnames(out2) = colnames(g1)
    rownames(out2) = names(out)
    return(out2)
}

library(ggplot2)
library(ggpubr)
g1_bin = BIN2(g1)
g2_bin = BIN2(g2)
AlleleFrac_bin = g1_bin/(g1_bin + g2_bin)
AlleleFrac_bin[(g1_bin+g2_bin) < 10] = NA #remove bins with <10 genoinformative transcripts
binUse = which(abs(rowMeans(AlleleFrac_bin, na.rm=T) - .5) < .4)  # Exclude bins with >90% of all transcripts mapping to the same allele across all samples
AlleleFrac_bin[-binUse,] = NA
FracMono_all = 100*colMeans(abs(AlleleFrac_bin - .5) >= .3, na.rm=T)

pseudocount = 1*10^6/quantile(colSums(D), p = .1)
A2 = sweep(D, 2, colSums(D), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform


A2d = A2b[rowSums(D[,colnames(A2b)] >= 10) >= 10, ]  # Require each gene to have at least 10 UMIs in at least 10 cells
fano = apply(A2d, 1, var)/rowMeans(A2d)  # fano factor is a measure of gene variance
hmat = A2d[rank(-fano[rownames(A2d)]) <= 500,over30k]
minmax = function(x) {
	sweep(x - log(pseudocount, 10), 1, apply(x - log(pseudocount, 10), 1, max), '/')
}

library(circlize)
HapFrac_col = colorRamp2(c(0.5, 1), c("white", "purple4"))
UMI_col = colorRamp2(c(4.4, 5.3), c("white", "forestgreen"))

source('Pseudotime Velocity Functions.R')

stages = as.character(At_meta[over30k,9])
names(stages) = over30k
stages[which(stages == "UM" | stages == "BM")] = "UM/BM"
stages[names(which(HapFrac_all[names(which(stages == "UM/BM"))] < .7))] = "UM/BM_bi"
stages[names(which(HapFrac_all[names(which(stages == "UM/BM"))] > .7))] = "UM/BM_mono"
stages[names(which(HapFrac_all[names(which(stages == "Tri"))] < .7))] = "Tri_bi"
stages[names(which(HapFrac_all[names(which(stages == "Tri"))] > .7))] = "Tri_mono"
stages = factor(stages, levels = c('tetrads','UM/BM_bi','UM/BM_mono','BM/Tri','Tri_bi', 'Tri_mono'))

pseudotime = function (x, stages) 
{
    PCs = pcaFn(x)
    pT = principal_curve(PCs)$lambda
    pT = pT * sign(cor(pT, stages[rownames(x)]))
    names(pT) = rownames(x)
    return(pT)
}

pVel = function (xx, nboot = 200) 
{
    xx = t(xx)
    stages = as.numeric(stages)
    names(stages) = rownames(xx)
    set.seed(1)
    out = NULL
    for (i in 1:nboot) {
        samps = sample(rownames(xx), replace = T)
        pT = pseudotime(xx[samps, ], stages)
        #pT = pseudotime(xx[samps, ])


        #pT[, 2] = pT[, 2] - mean(apply(pT, 1, diff), na.rm = T)
        #pT = rowMeans(pT, na.rm = T)
        pT = 100 * (pT - quantile(pT, 0.1))/diff(quantile(pT, c(0.1, 0.9)))
        names(pT) = samps
        out = c(out, pT)
        if (i == 1) {
            cat("0% |---------|---------|---------|---------| 100%\n   .")
        }
        if ((i%%round(nboot/40)) == 0) {
            cat(".")
        }
    }
    pT = by(out, names(out), mean)[rownames(xx)]
    pV = c(rep(NA, 4), velFn(pT[order(pT)]), rep(NA, 4))
    pV = pV/median(pV, na.rm = T)
    return(list(pT = pT, pV = pV))
}

pseudos = pVel(A2d[rank(-fano) <= 2000,over30k], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords = order(stages, pseudos$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime

hmat2 = A2d[,over30k][rank(-fano[rownames(A2d)]) <= 500,ords]

library(seriation)
o1 = seriate(dist(t(scale(t(hmat2)))), method = "OLO")

Heatmap(minmax(hmat2), name = 'expression
level', 
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    HapFrac = HapFrac_all[over30k][ords], UMIcounts = log(colSums(D[,over30k][,ords]),10), Yellow = At_meta[over30k,13][ords],
    col = list(Stage = c("tetrads" = "#EB1E2C", "UM/BM_bi" = "#F9A729", "UM/BM_mono" = "#F9D23C", "BM/Tri" = "#5FBB68", 
    "Tri_bi" = "#64CDCC", "Tri_mono" = "#A4A4D5"), HapFrac = HapFrac_col, UMIcounts = UMI_col, Yellow = c("Y" = "yellow", "N" = "white"))),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_rows=as.dendrogram(o1[[1]]), cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)

col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))
Heatmap(X_samp[Ranked500G,], name = '% transcipts from Col-0 allele', col = col_fun, cluster_columns=F, cluster_rows=F, show_row_names = FALSE, show_column_names = FALSE, use_raster=F)
Heatmap(Y_samp[Ranked500G,], name = '% transcripts matching haplotype', col = col_fun, cluster_columns=F, cluster_rows=F, show_row_names = FALSE, show_column_names = FALSE, use_raster=F)

HapFrac_row = 100*rowMeans(abs(HapFrac_genes_all[filtG,] - .5) >= .3, na.rm=T)
HapFrac_row[which(rowSums(!is.na(Y_samp[filtG,])) < 5)] = NA
#HapFrac_row_col = colorRamp2(c(0, 25, 50, 75, 100), c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))
HapFrac_row_col = colorRamp2(c(0, 100), c("white", "purple4"))
HapFrac_row_bi = 100*rowMeans(abs(HapFrac_genes_all[filtG,names(which(HapFrac_all < .65))] - .5) >= .3, na.rm=T)
HapFrac_row_bi[which(rowSums(!is.na(Y_samp[filtG,names(which(HapFrac_all < .65))])) < 10)] = NA


Heatmap(X_samp[Ranked500G,ords], name = '% transcipts 
from Col-0 allele', 
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    HapFrac = HapFrac_all[over30k][ords], UMIcounts = log(colSums(D[,over30k][,ords]),10), Yellow = At_meta[over30k,13][ords],
    col = list(Stage = c("tetrads" = "#EB1E2C", "UM/BM_bi" = "#F9A729", "UM/BM_mono" = "#F9D23C", "BM/Tri" = "#5FBB68", 
    "Tri_bi" = "#64CDCC", "Tri_mono" = "#A4A4D5"), HapFrac = HapFrac_col, UMIcounts = UMI_col, Yellow = c("Y" = "yellow", "N" = "white"))),
    col = col_fun, cluster_rows=F, cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)

Heatmap(X_samp[Cands,ords], name = '% transcipts 
from Col-0
allele', 
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    HapFrac = HapFrac_all[over30k][ords], UMIcounts = log(colSums(D[,over30k][,ords]),10), Yellow = At_meta[over30k,13][ords],
    col = list(Stage = c("tetrads" = "#EB1E2C", "UM/BM_bi" = "#F9A729", "UM/BM_mono" = "#F9D23C", "BM/Tri" = "#5FBB68", 
    "Tri_bi" = "#64CDCC", "Tri_mono" = "#A4A4D5"), HapFrac = HapFrac_col, UMIcounts = UMI_col, Yellow = c("Y" = "yellow", "N" = "white"))),
    col = col_fun, cluster_rows=F, cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)

Heatmap(Y_samp[Cands,ords], name = '% transcipts 
matching
haplotype', 
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    HapFrac = HapFrac_all[over30k][ords], UMIcounts = log(colSums(D[,over30k][,ords]),10), Yellow = At_meta[over30k,13][ords],
    col = list(Stage = c("tetrads" = "#EB1E2C", "UM/BM_bi" = "#F9A729", "UM/BM_mono" = "#F9D23C", "BM/Tri" = "#5FBB68", 
    "Tri_bi" = "#64CDCC", "Tri_mono" = "#A4A4D5"), HapFrac = HapFrac_col, UMIcounts = UMI_col, Yellow = c("Y" = "yellow", "N" = "white"))),
    col = col_fun, cluster_rows=F, cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)

Heatmap(Y_samp[,ords], name = '% transcipts 
matching 
haplotype', right_annotation = rowAnnotation("HapFrac row bi" = HapFrac_row_bi, 
    col = list("HapFrac row bi" = HapFrac_row_col)),
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    HapFrac = HapFrac_all[over30k][ords], UMIcounts = log(colSums(D[,over30k][,ords]),10), Yellow = At_meta[over30k,13][ords],
    col = list(Stage = c("tetrads" = "#EB1E2C", "UM/BM_bi" = "#F9A729", "UM/BM_mono" = "#F9D23C", "BM/Tri" = "#5FBB68", 
    "Tri_bi" = "#64CDCC", "Tri_mono" = "#A4A4D5"), HapFrac = HapFrac_col, UMIcounts = UMI_col, Yellow = c("Y" = "yellow", "N" = "white"))),
    col = col_fun, cluster_rows=F, cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)

plot(colMeans(X_samp[,ords], na.rm=T), colMeans(Y_samp[,ords], na.rm=T), cex=2, pch=19)



#Calc HapExpScore from UM/BM_bi and UM/BM_mono
ratio = g1/(g1+g2)
ratio[(g1+g2) < 5] = NA  # Remove measurements with under 5 genoinformative transcripts
ratio = ratio[,over30k] # Remove samples with <25,000 UMIs (over30k name is legacy, just hasn't been changed)
#ratio[which(abs(rowMeans(ratio, na.rm=T)-0.5) >=0.3),] = NA #Remove genes with > 80% transcripts matching one allele

BiasedG = names(which(abs(rowMeans(ratio, na.rm=T)-0.5) >=0.3))
BiasedFiltG = BiasedG[which(BiasedG %in% filtG)]

BootSamps = names(stages[stages == "UM/BM_bi"])
filtG = names(which(rowSums(!is.na(ratio[,BootSamps])) >= 5))
filtG = names(which(rowSums(!is.na(Hap_all[filtG,BootSamps])) >= 5))

Boot = function(x){
    samps = sample(BootSamps, replace=TRUE)
    X = rowMeans(ratio[filtG, samps], na.rm=T)
    Y = rowMeans(((ratio - .5)*Hap_all)[filtG, samps], na.rm=T)+0.5
   cbind(X,Y)
}
set.seed(1)
Boots = sapply(1:10000, Boot)
Boots = lapply(list(X = Boots[1:length(filtG),], Y = Boots[(1:length(filtG)) + length(filtG),]), function(xx) { 
        rownames(xx) = filtG
        return(xx)
    })

HapExpScore = (Boots$Y - abs(Boots$X - 0.5) - 0.5) / (1 - abs(Boots$X - 0.5) - 0.5)
Frac95CI = apply(HapExpScore, 1, quantile, probs=c(0.025,0.975), na.rm=T)

save(HapExpScore, Frac95CI, Boots, Boot, ratio, D, g1, g2, genes, BootSamps, filtG, X_samp, Y_samp, At_meta, stages, ords, over30k, file = "HapExpScore10k.RData")
load("HapExpScore10k.RData")

summary(rowMeans(HapExpScore, na.rm=T))
#library('ComplexHeatmap')
#library(circlize)
#col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))
#Heatmap(HapExpScore, col = col_fun, cluster_columns=T, cluster_rows=T, show_row_names = FALSE, show_column_names = FALSE)
plot(rowMeans(HapExpScore, na.rm=T), cex=2, pch=19)
hist(rowMeans(HapExpScore, na.rm=T), breaks=100)

Frac95CI = apply(HapExpScore, 1, quantile, probs=c(0.025,0.975), na.rm=T)
plot(Frac95CI['2.5%',LateBiCand], Frac95CI['97.5%',LateBiCand], cex=2, pch=19, xlim=c(0,1), ylim=c(0,1))
abline(a=0, b=1, col="red2")

plot(Frac95CI['2.5%',], Frac95CI['97.5%',], cex=2, pch=19)
abline(a=0, b=1, col="red2")

summary(Frac95CI['97.5%',]-Frac95CI['2.5%',])
hist(Frac95CI['97.5%',]-Frac95CI['2.5%',], breaks=100)    
hist(Frac95CI['97.5%',]-Frac95CI['2.5%',], breaks=100, xlim=c(0,1))



BootSamps = names(stages[stages == "UM/BM_bi"])
filtG = names(which(rowSums(!is.na(ratio[,BootSamps])) >= 5))
filtG = names(which(rowSums(!is.na(Hap_all[filtG,BootSamps])) >= 5))

Boot = function(x){
    samps = sample(BootSamps, replace=TRUE)
    X = rowMeans(ratio[filtG, samps], na.rm=T)
    Y = rowMeans(((ratio - .5)*Hap_all)[filtG, samps], na.rm=T)+0.5
   cbind(X,Y)
}

set.seed(1)
a1 = proc.time()
Boots = sapply(1:(10), Boot)
Boots = lapply(list(X = Boots[1:length(filtG),], Y = Boots[(1:length(filtG)) + length(filtG),]), function(xx) { 
        rownames(xx) = filtG
        return(xx)
    })
a2 = proc.time()
a2-a1 
#   user  system elapsed
#   0.98    0.40    1.69

HapExpScore = (Boots$Y - abs(Boots$X - 0.5) - 0.5) / (1 - abs(Boots$X - 0.5) - 0.5)
Frac95CI = apply(HapExpScore, 1, quantile, probs=c(0.025,0.975), na.rm=T)
#EarlyMonoCand = names(which(Frac95CI['2.5%',] >= 0.1))
EarlyMonoCand = names(which(rowMeans(HapExpScore <= 0, na.rm=T) == 0))
hist(rowMeans(HapExpScore[EarlyMonoCand,] <= 0, na.rm=T))
hist(rowMeans(HapExpScore <= 0, na.rm=T)) #p value. For 971 samples significance is .05/971 = ~5x10^-5


#Alt Boot only including samples with info for each gene
BootSamps = names(stages[stages == "UM/BM_bi"])
filtG = names(which(rowSums(!is.na(ratio[,BootSamps])) >= 5))
filtG = names(which(rowSums(!is.na(Hap_all[filtG,BootSamps])) >= 5))
SampsWithInfo = !is.na(ratio[filtG,BootSamps])
BootAlt = function(x){
    X = matrix(data=NA, nrow=length(filtG), ncol=1)
    Y = X
    for(g in 1:length(filtG)){
        samps = sample(names(which(SampsWithInfo[g,])), replace=TRUE)
        X[g] = mean(ratio[filtG, samps], na.rm=T)
        Y[g] = mean(((ratio - .5)*Hap_all)[filtG, samps], na.rm=T)+0.5
    }
    cbind(X,Y)
}
set.seed(1)
a1 = proc.time()
BootAlts = sapply(1:(10), BootAlt)
BootAlts = lapply(list(X = BootAlts[1:length(filtG),], Y = BootAlts[(1:length(filtG)) + length(filtG),]), function(xx) { 
        rownames(xx) = filtG
        return(xx)
    })
a2 = proc.time()
a2-a1
#   user  system elapsed
# 275.22  158.05  437.50

    


BootSamps = names(stages[stages == "UM/BM_mono"])
filtG = names(which(rowSums(!is.na(ratio[,BootSamps])) >= 5))
filtG = names(which(rowSums(!is.na(Hap_all[filtG,BootSamps])) >= 5))
set.seed(1)
Boots = sapply(1:2000, Boot)
Boots = lapply(list(X = Boots[1:length(filtG),], Y = Boots[(1:length(filtG)) + length(filtG),]), function(xx) { 
        rownames(xx) = filtG
        return(xx)
    })
HapExpScore = (Boots$Y - abs(Boots$X - 0.5) - 0.5) / (1 - abs(Boots$X - 0.5) - 0.5)
Frac95CI = apply(HapExpScore, 1, quantile, probs=c(0.025,0.975), na.rm=T)
LateBiCand = names(which(Frac95CI['97.5%',] <= 0.3))
LateBiCand = names(which(Frac95CI['2.5%',LateBiCand] >= -0.3))
hist(rowMeans(HapExpScore[LateBiCand,] > 0.3, na.rm=T))

Cands = c(EarlyMonoCand, LateBiCand)
Cands = c(EarlyMonoCand, BiasedFiltG)
Heatmap(X_samp[Cands,ords], name = '% transcipts 
from Col-0
allele', 
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    HapFrac = HapFrac_all[over30k][ords], UMIcounts = log(colSums(D[,over30k][,ords]),10), Yellow = At_meta[over30k,13][ords],
    col = list(Stage = c("tetrads" = "#EB1E2C", "UM/BM_bi" = "#F9A729", "UM/BM_mono" = "#F9D23C", "BM/Tri" = "#5FBB68", 
    "Tri_bi" = "#64CDCC", "Tri_mono" = "#A4A4D5"), HapFrac = HapFrac_col, UMIcounts = UMI_col, Yellow = c("Y" = "yellow", "N" = "white"))),
    col = col_fun, cluster_rows=F, cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)

Heatmap(Y_samp[Cands,ords], name = '% transcipts 
matching
haplotype', 
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    HapFrac = HapFrac_all[over30k][ords], UMIcounts = log(colSums(D[,over30k][,ords]),10), Yellow = At_meta[over30k,13][ords],
    col = list(Stage = c("tetrads" = "#EB1E2C", "UM/BM_bi" = "#F9A729", "UM/BM_mono" = "#F9D23C", "BM/Tri" = "#5FBB68", 
    "Tri_bi" = "#64CDCC", "Tri_mono" = "#A4A4D5"), HapFrac = HapFrac_col, UMIcounts = UMI_col, Yellow = c("Y" = "yellow", "N" = "white"))),
    col = col_fun, cluster_rows=F, cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)




#calc p for each gene. Use simulation to figure out how many samples need to have gene info for significance
nP = 3
simSamps = t(sapply(1:1000, function(i) { sample(0:1, replace = T, nP) }))
rowMeans(abs(simSamps - 0.5)+ 0.5)
t.test(abs(simSamps[1,]-.5), mu=0.5)

simBoot = function(x, nP=3){
    simSamps = t(sapply(1:1000, function(i) { sample(0:1, replace = T, nP) }))
    samps = sample(1:nP, replace = T, nP)
    X = rowMeans(simSamps[,samps])
    Y = rowMeans((simSamps - .5)[,samps])+.5
    cbind(X,Y)
}

simBoots = sapply(1:100, simBoot)
simBoots = lapply(list(X = simBoots[1:1000,], Y = simBoots[1:1000 + 1000,]), function(xx) {
    rownames(xx) = 1:1000
    return(xx)
})
simBoots$Y - abs(simBoots$X - 0.5)

summary(apply((simBoots$Y - abs(simBoots$X - 0.5)), 1, mean)) 
Frac95CI = apply((simBoots$Y - abs(simBoots$X - 0.5)), 1, quantile, probs=c(0.025,0.975), na.rm=T)
plot(Frac95CI['2.5%',], Frac95CI['97.5%',], cex=2, pch=19, xlim=c(0,1), ylim=c(0,1))
abline(a=0, b=1, col="red2")
plot(Frac95CI['97.5%',]-Frac95CI['2.5%',], apply(fracsBootMod, 1, mean), cex=2, pch=19, xlim=c(0,0.5), ylim=c(0,1))


BootFracs = function(x){
    samps = sample(BootSamps, replace=TRUE)
    X = rowMeans(ratio[filtG, samps], na.rm=T)
    Y = rowMeans(((ratio - .5)*Hap_all)[filtG, samps], na.rm=T)+0.5
    SampFracs = data.frame(gene=names(X), bias=X, frac=Y)
    SampFrac = SampFracs$frac
    names(SampFrac) = SampFracs$gene
    return(SampFrac)
}

BootFracsMod = function(x){
    samps = sample(BootSamps, replace=TRUE)
    X = rowMeans(ratio[filtG, samps], na.rm=T)
    Y = rowMeans(((ratio - .5)*Hap_all)[filtG, samps], na.rm=T)+0.5
    SampFracs = data.frame(gene=names(X), bias=X, frac=Y)
    SampFrac = SampFracs$frac - abs(SampFracs$bias - 0.5)
    names(SampFrac) = SampFracs$gene
    return(SampFrac)
}

fracsBoot = sapply(1:100, BootFracs)

fracsBootMod = sapply(1:100, BootFracsMod)
summary(apply(fracsBootMod, 1, mean)) 
Frac95CI = apply(fracsBootMod, 1, quantile, probs=c(0.025,0.975), na.rm=T)
plot(Frac95CI['2.5%',], Frac95CI['97.5%',], cex=2, pch=19, xlim=c(0,1), ylim=c(0,1))
abline(a=0, b=1, col="red2")
plot(Frac95CI['97.5%',]-Frac95CI['2.5%',], apply(fracsBootMod, 1, mean), cex=2, pch=19, xlim=c(0,0.5), ylim=c(0,1))


summary((apply(fracsBoot, 1, sd)))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
 0.0115  0.0527  0.0669  0.0724  0.0850  0.2453    1442

 For comparison, here is the frac value summary from the whole group of samples:
summary(calcFracs2("UM/BM")$frac)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0003568 0.6000000 0.7000000 0.6828712 0.7712821 1.0000000

And here is the frac value summary produced with BootFracs:
fracsBootMean = apply(fracsBoot, 1, mean)
summary(fracsBootMean)                                  
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.3988  0.6326  0.7040  0.6798  0.7470  0.9008

Frac95CI = apply(fracsBoot, 1, quantile, probs=c(0.025,0.975), na.rm=T)
plot(Frac95CI['2.5%',], Frac95CI['97.5%',], cex=2, pch=19, xlim=c(0,1), ylim=c(0,1))
abline(a=0, b=1, col="red2")

plot(Frac95CI['97.5%',]-Frac95CI['2.5%',], fracsBootMean, cex=2, pch=19, xlim=c(0,0.5), ylim=c(0,1))

#which(rowSums(is.na(fracsBoot)) < 80)

> sd(zxBoot) -> sem
> quantile(zxBoot, p = c(.025,.975)) -> ci95
> mean(zxBoot <= 0) -> p value for true population mean less than 0

duplicated(fracs[which(fracs$class == ">.8"),1]) #Tried this with <0.6 HapFrac, no duplicates
#with 0.55-0.75 HapFrac, ID=gene-AT1G07670 is a duplicate in UM and BM. 
fracs[which(fracs$class == ">.95"),1] %in% fracs[which(fracs$class == ">.8"),1]
#with 0.55-0.75 HapFrac, no matches between >.8 and >.95
#with 0.6-0.8 HapFrac, 5 matches between >.8 and >.95, ID=gene-AT1G44191 and ID=gene-AT5G65687 are >.95 in UM/BM stage and >.8 in UM and BM stages


fracs[which(fracs$class == ">.95"),][which(fracs[which(fracs$class == ">.95"),1] %in% fracs[which(fracs$class == ">.8"),1]),]