#Prep_for_ortholog_comp
###In local R terminal
setwd('C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R')
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(circlize)
library(seriation)
library('ComplexHeatmap')
source('Pseudotime Velocity Functions.R')
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
minmax = function(x) {
	sweep(x - log(pseudocount, 10), 1, apply(x - log(pseudocount, 10), 1, max), '/')
}


###Rice prep
load("RiceRemap_7-8_2025.RData")

Rice_Stages <- read_excel("C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R/Rice_Stages.xlsx")
Rice_meta <- Rice_Stages[rep(row.names(Rice_Stages), each = 8), ]
Rice_meta = as.data.frame(Rice_meta)
rownames(Rice_meta) = paste(rep(Rice_Stages$Plate, each = 8), '_', rep((as.numeric(sub('[A-Z].+','',Rice_Stages$Celseq_primer)) - 1)*8, each = 8) + rep(1:8, nrow(Rice_Stages)), 's', sep ='')

Rice_meta[which(Rice_meta[,10] == "shed"),12] = "TP_shed"
Rice_meta[grep("sep only", Rice_meta[,10]),12] = "tetrad_sep"
Rice_meta[grep("sep,", Rice_meta[,10])[c(6:8,13:16,23:24)],12] = "tetrad_sep"

Rice_meta$No_cell2 = unlist(lapply(Rice_Stages$No_cell, function(xx) {       
        out = rep('N', 8)
        out[as.numeric(strsplit(xx,',')[[1]])] = 'Y'
        return(out) }))

D2 = D2[,rownames(Rice_meta)]
g1_2 = g1_2[,rownames(Rice_meta)]
g2_2 = g2_2[,rownames(Rice_meta)]

BINR = function (xx, bin = 10^6) 
{
    bin = as.numeric(Rice_genes2[, 1]) * 10^6 + round(Rice_genes2[, 2]/bin)
    out = by(xx, bin, colSums)
    out2 = t(matrix(unlist(out), nrow = ncol(g1_2)))
    colnames(out2) = colnames(g1_2)
    rownames(out2) = names(out)
    return(out2)
}

g1_bin = BINR(g1_2)
g2_bin = BINR(g2_2)
AlleleFrac_bin = g1_bin/(g1_bin + g2_bin)
AlleleFrac_bin[(g1_bin+g2_bin) < 10] = NA #remove bins with <10 genoinformative transcripts
R_binUse = which(abs(rowMeans(AlleleFrac_bin, na.rm=T) - .5) < .4)  # Exclude bins with >90% of all transcripts mapping to the same allele across all samples
AlleleFrac_bin[-R_binUse,] = NA
FracMono_all = 100*colMeans(abs(AlleleFrac_bin - .5) >= .3, na.rm=T)
FracMono_all[names(which(colSums(!is.na(AlleleFrac_bin))<10))] = NA # Exclude samples with <10 informative bins

over10k = names(which(colSums(D2) >15000))
over10k = over10k[which(Rice_meta[over10k,14] == "N")] #getting rid of the no cell controls (4 at 20k, 26 at 10k) #14 for newer table, 12 for old table

stages = Rice_meta[over10k,12]
names(stages) = over10k
RiceBiSamp = names(which(FracMono_all[over10k] <= 40))
RiceMonoSamp = names(which(FracMono_all[over10k] > 40))
stages[names(stages[grep("early UM", stages)])] = "UM"
stages[names(stages[grep("mid UM", stages)])] = "UM"
stages[names(stages[grep("late UM", stages)])] = "UM"
stages[grep("tetrad", stages)] = "Tetrad"
stages[names(which(stages[RiceBiSamp] == "UM"))] = "UM Bi"
stages[names(which(stages[RiceMonoSamp] == "UM"))] = "UM Mono"
#stages[stages == "TP_shed"] = "TP"
stages = factor(stages, levels = c("Tetrad", "UM Bi", "UM Mono", "UM/BM", "BM", "BM/TP", "TP", "TP_shed"))


Anther_Frac = data.frame(FracMono = FracMono_all[over10k], Anther = as.factor(Rice_meta[over10k,5]), Stage = as.character(stages))
Anther_Frac$Range = NA
for (x in 1:88){
    Anther_Frac[which(Anther_Frac$Anther == unique(Anther_Frac$Anther)[x]),4] = diff(range(Anther_Frac[Anther_Frac$Anther == unique(Anther_Frac$Anther)[x],1]))
}
Anther_Frac$Stage[which(Anther_Frac$Stage == "UM Bi" | Anther_Frac$Stage == "UM Mono")] = "UM"       
table(table(Anther_Frac$Anther[which(Anther_Frac$Stage == "UM")]))
# 0  1  2  3  4  5  6 
#47 18  6  5  5  4  3 #18 anthers with only one sample


ggplot(Anther_Frac, aes(x = Anther, y = FracMono, fill = Stage, color = Stage)) + geom_violin(linewidth=0.1, adjust = 0.3) + labs(title = "Rice", 
x = "Anther", y = "% Monoallelic
    Transcripts") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none", 
    panel.grid=element_blank(), axis.ticks.x = element_blank()) + scale_fill_manual(values=c("Tetrad" = "#EB1E2C",
    "UM" = "#F9A729","UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5")) + scale_color_manual(values=c("Tetrad" = "#EB1E2C",
    "UM" = "#F9A729", "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5")) + ylim(0,100)


stages2 = as.character(stages)
names(stages2) = names(stages)
stages2[rownames(Anther_Frac[Anther_Frac$Stage == "UM" & Anther_Frac$Range > 30,])] = "UM Mix"
stages2 = factor(stages2, levels = c("Tetrad", "UM Bi", "UM Mix", "UM Mono", "UM/BM", "BM", "BM/TP", "TP", "TP_shed"))


pseudocount = 1*10^6/quantile(colSums(D2), p = .1)
A2 = sweep(D2, 2, colSums(D2), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(D2[,colnames(A2b)] >= 5) >= 3, ]  # Require each gene to have at least 5 UMIs in at least 3 cells

fano = apply(A2d, 1, var)/rowMeans(A2d)  # fano factor is a measure of gene variance
pseudos = pVel(A2d[rank(-fano) <= 2000,over10k], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords = order(stages, pseudos$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime
hmap = A2d[,over10k][rank(-fano[rownames(A2d)]) <= 500,ords] #for plotting heatmap
o1 = seriate(dist(t(scale(t(hmap)))), method = "OLO") #clustering genes in heatmap
#R37-90_7s, R47-102_9s

stage_means <- sapply(levels(stages), function(s) {
    rowMeans(A2d[, stages == s, drop = FALSE], na.rm = TRUE)
})

colnames(stage_means) <- levels(stages)
cor(stage_means[, "UM Bi"], stage_means[, "TP"], method = "spearman")
stage_z <- t(scale(t(stage_means)))

endpoint_score <- (
    stage_z[, "UM Bi"] +
    stage_z[, "TP"]
) / 2 -
(
    stage_z[, "UM Mono"] +
    stage_z[, "UM/BM"] +
    stage_z[, "BM"] +
    stage_z[, "BM/TP"]
) / 4


endpoint_score <- (
    stage_z[, "UM Bi"] +
    stage_z[, "BM"]
) / 2 -
(
    stage_z[, "UM Mono"] +
    stage_z[, "UM/BM"]
) / 2



problem_genes <- names(sort(endpoint_score[genes_original], decreasing = TRUE))
genes_filtered <- setdiff(genes_original,problem_genes[1:500])



stage_cor <- apply(stage_means, 1, function(x) {
    cor(x, 1:8, method = "spearman")
})

genes_original <- names(sort(fano, decreasing = TRUE))[1:2000]
#rank(-fano) <= 2000
genes_filtered0 <- names(stage_cor[stage_cor > 0.2 | stage_cor < -0.2])
fano2 = fano[genes_filtered0]
genes_filtered = names(sort(fano2, decreasing = TRUE))[1:2000]

pseudos_filtered <- pVel(
    A2d[genes_filtered, over10k],
    nboot = 80
)
boxplot(pseudos_filtered$pT ~ stages,las = 2,xlab = "Stage",ylab = "Pseudotime")

boxplot(pseudos$pT ~ stages,las = 2,xlab = "Stage",ylab = "Pseudotime")

boxplot(Os_pT_stage$pT ~ Os_pT_stage$stages,las = 2,xlab = "Stage",ylab = "Pseudotime")


flagged_genes <- c()

for(gene in colnames(t(A2d[rank(-fano) <= 2000,over10k][,ords]))) {
  # Fit a smooth spline to smooth out cell-level dropout noise
  fit <- predict(smooth.spline(pseudos$pT, t(A2d[rank(-fano) <= 2000,over10k][,ords])[,gene]), x = pseudos$pT)
  
  # Calculate differences to find peaks and troughs
  dy <- diff(fit$y)
  
  # Count how many times the direction changes sign
  sign_changes <- sum(diff(sign(dy)) != 0)
  
  # If it flips direction more than twice, it's an alternating pattern
  if(sign_changes > 2) {
    flagged_genes <- c(flagged_genes, gene)
  }
}

# Remove flagged genes from your downstream trajectory analysis
clean_matrix <- expression_matrix[, !(colnames(expression_matrix) %in% flagged_genes)]




MeanByStage = simplify2array(by(t(A2d[,over10k]), stages, colMeans))
MeanByStage = MeanByStage[,1:7]
MeanDiffByStage = sapply(2:length(colnames(MeanByStage)), function(x){ abs(MeanByStage[,x] - MeanByStage[,x-1])})
plot(MeanDiffByStage[names(which(rank(-fano[rownames(A2d)]) <= 2000)),][1,])

MeanByStage2 = MeanByStage[rank(-fano) <= 2000,]

plot(MeanByStage2[1,], ylim = c(2.6,4.45))
sapply(201:300, function(x){
    lines(MeanByStage2[x,])
})

hist(sapply(1:2000, function(x){ abs(MeanByStage2[x,5] - MeanByStage2[x,2])/(max(MeanByStage2[x,]) - min(MeanByStage2[x,])) }))
keepers = which(sapply(1:2000, function(x){ abs(MeanByStage2[x,7] - MeanByStage2[x,2])/(max(MeanByStage2[x,]) - min(MeanByStage2[x,])) }) < 0.9)
fano2 = fano[rank(-fano) <= 2000][keepers]

pseudos = pVel(A2d[names(fano2),over10k], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords = order(stages, pseudos$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime
hmap = A2d[names(fano2),over10k][rank(-fano2) <= 500,ords]
o1 = seriate(dist(t(scale(t(hmap)))), method = "OLO")


Os_Z <- t(scale(t(A2d[rank(-fano[rownames(A2d)]) <= 500,over10k])))

# 1. Calculate the distance matrix (Euclidean distance is standard)
gene_dist <- dist(Os_Z, method = "euclidean")

# 2. Perform hierarchical clustering (complete linkage is commonly used)
gene_hclust <- hclust(gene_dist, method = "complete")

# 3. Cut the tree into a specific number of clusters (e.g., k = 6)
gene_clusters <- cutree(gene_hclust, k = 15)

# View cluster assignments
head(gene_clusters)

simplify2array(by(t(A2d[,over10k]), stages, colMeans))
MeanByStage3 = simplify2array(by(t(Os_Z), stages, colMeans))[,1:7]

MeanByStage[rank(-fano) <= 2000,]

plot(MeanByStage3[1,], ylim = c(2.6,4.45))
sapply(201:300, function(x){
    lines(MeanByStage3[x,])
})

row_split = gene_clusters

minmax2 = function(x) {
sweep(x - 1, 1, apply(x - 1, 1, max), '/')
}

Heatmap(minmax(hmap), row_split = gene_clusters, name = 'Z-score', #heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")), #column_split = split,
    top_annotation = HeatmapAnnotation(Stage = stages[ords], "UMIcounts" = log(colSums(D2[,over10k][,ords]),10),
    "% Monoallelic" = FracMono_all[over10k][ords], 
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM Bi" = "#F9A729", "UM Mono" = "#F9A729",
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", "TP_shed" = "gray"), "% Monoallelic" = FracMono_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_rows=T, cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)
#try removing cluster1
keepers = names(which(gene_clusters > 1))
fano2 = fano[rank(-fano) <= 500][keepers]
pseudos = pVel2(A2d[names(fano2),over10k], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords = order(stages, pseudos$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime
hmap = A2d[names(fano2),over10k][rank(-fano2) <= 500,ords]
o1 = seriate(dist(t(scale(t(hmap)))), method = "OLO")

ords2 = order(pseudos$pT)
hmap = A2d[,over10k][rank(-fano[rownames(A2d)]) <= 500,ords2] #for plotting heatmap

#hmap2 = A2d[names(fano),over10k][rank(-fano) <= 500,ords2]

Heatmap(minmax(hmap), name = 'expression
level', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")), #column_split = split,
    top_annotation = HeatmapAnnotation(Stage = stages2[ords2], "UMIcounts" = log(colSums(D2[,over10k][,ords2]),10),
    "% Monoallelic" = FracMono_all[over10k][ords2], "Anther length" = Rice_meta[over10k,6][ords2],
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM Bi" = "#F9A729", "UM Mix" = "#333", "UM Mono" = "#F9A729",
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D9", "TP_shed" = "gray"), "% Monoallelic" = FracMono_col,
    "UMIcounts" = UMI_col, "Anther length" = Anther_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_rows=as.dendrogram(o1[[1]]), cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)




FracMono_col = colorRamp2(c(0, 100), c("white", "purple4"))
UMI_col = colorRamp2(c(4.0, 5.55), c("white", "forestgreen")) #summary(log(colSums(D2[,over10k]),10)) #to choose min and max

#manually splitting using table(stages) and heatmaps to guide
split = factor(rep(c("Os_Tetrad1","Os_Tetrad2","Os_UMbi1","Os_UMbi2","Os_UMbi3","Os_UMmono","Os_UM/BM","Os_BM","Os_BM/TP","Os_TP","Os_TP_shed"), c(23,10,8,32,30,33,10,24,10,124,3)), levels = c("Os_Tetrad1","Os_Tetrad2","Os_UMbi1","Os_UMbi2","Os_UMbi3","Os_UMmono","Os_UM/BM","Os_BM","Os_BM/TP","Os_TP","Os_TP_shed"))

Heatmap(minmax(hmap), name = 'expression
level', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")), #column_split = split,
    top_annotation = HeatmapAnnotation(Stage = stages[ords], "UMIcounts" = log(colSums(D2[,over10k][,ords]),10),
    "% Monoallelic" = FracMono_all[over10k][ords], 
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM Bi" = "#F9A729", "UM Mono" = "#F9A729",
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", "TP_shed" = "gray"), "% Monoallelic" = FracMono_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_rows=as.dendrogram(o1[[1]]), cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)

Heatmap(cor(hmap, use = 'pairwise.complete.obs', method = 'pearson'), 
name = 'Expression_cor', right_annotation = rowAnnotation(Stage = stages[ords], "UMIcounts" = log(colSums(D2[,ords]),10),
    "% Monoallelic" = FracMono_all[over10k][ords], 
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM Bi" = "#F9A729", "UM Mono" = "#F9A729",
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", "TP_shed" = "gray"), "% Monoallelic" = FracMono_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), show_row_names = FALSE, show_column_names = FALSE)


#ords2 = order(pseudos$pT)

pseudos2 = pVel2(A2d[rank(-fano) <= 2000,names(stages2)], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords2 = order(stages2, pseudos2$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime

#ords2 = order(stages2, Rice_meta[over10k,6])
hmap2 = A2d[,over10k][rank(-fano[rownames(A2d)]) <= 500,ords2]
o1_2 = seriate(dist(t(scale(t(hmap2)))), method = "OLO")

FracMono_col = colorRamp2(c(0, 100), c("white", "purple4"))
UMI_col = colorRamp2(c(4.0, 5.55), c("white", "forestgreen")) #summary(log(colSums(D2[,over10k]),10)) #to choose min and max
Anther_col = colorRamp2(c(0.8, 2), c("white", "black"))

svg("Rice_ExpHmap_read.svg", height = 60, width=8)
Heatmap(minmax(hmap2), name = 'expression
level', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")), #column_split = split,
    top_annotation = HeatmapAnnotation(Stage = stages2[ords2], "UMIcounts" = log(colSums(D2[,over10k][,ords2]),10),
    "% Monoallelic" = FracMono_all[over10k][ords2], "Anther length" = Rice_meta[over10k,6][ords2],
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM Bi" = "#F9A729", "UM Mix" = "#333", "UM Mono" = "#F9A729",
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D9", "TP_shed" = "gray"), "% Monoallelic" = FracMono_col,
    "UMIcounts" = UMI_col, "Anther length" = Anther_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_rows=as.dendrogram(o1_2[[1]]), cluster_columns=F, show_row_names = TRUE, show_column_names = FALSE)
dev.off()

Heatmap(minmax(hmap[,rownames(Os_pT_stage)]), name = 'expression
level', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")), #column_split = split,
    top_annotation = HeatmapAnnotation(Stage = stages2[rownames(Os_pT_stage)], "UMIcounts" = log(colSums(D2[,over10k][,rownames(Os_pT_stage)]),10),
    "% Monoallelic" = FracMono_all[over10k][rownames(Os_pT_stage)], "Anther length" = Rice_meta[over10k,][rownames(Os_pT_stage),6],
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM Bi" = "#F9A729", "UM Mix" = "#333", "UM Mono" = "#F9A729",
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D9", "TP_shed" = "gray"), "% Monoallelic" = FracMono_col,
    "UMIcounts" = UMI_col, "Anther length" = Anther_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_rows=as.dendrogram(o1[[1]]), cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)



splitSamps = data.frame(split = split, samp = names(stages[ords]))
calcMeanExp = function(x){
    rowMeans(A2[,splitSamps[splitSamps[,1] == x,2]]) #mean TPM of each gene for each group
}
calcMeanExpL = function(x){
    rowMeans(log(A2[,splitSamps[splitSamps[,1] == x,2]]+100,10)) #mean log(TPM+100,10) of each gene for each group
}

###Objects to use downstream in comparisons:
Os_meanExp = sapply(levels(split)[which(table(split)>=3)], calcMeanExp) #apply meanExp only to the groups with >=3 samples
Os_meanExpL = sapply(levels(split)[which(table(split)>=3)], calcMeanExpL) #apply meanExpL only to the groups with >=3 samples
Os_UMIs = D2[,names(stages)] #min 10k UMIs, excluding no cell controls, 307 samples
Os_g1 = g1_2[,names(stages)]
Os_g2 = g2_2[,names(stages)]
Os_Meta = Rice_meta[names(stages),]
Os_stages = stages
Os_AlleleFracBin = AlleleFrac_bin[,names(stages)]
Os_FracMono = FracMono_all[names(stages)]
Os_TPM = A2[,names(stages)]
Os_genes = Rice_genes2 #need to make sure Os_genes is correct as it was run set to =genes

Os_pT_stage = data.frame(pT = pseudos$pT[ords], stages = as.character(stages[ords]))
Os_pT_stage[grep("TP_shed", Os_pT_stage$stage),2] = "TP"
Os_pT_mean2 = sapply(c("Tetrad", "UM Bi", "UM Mono", "UM/BM", "BM", "BM/TP", "TP"), function(x){ mean(Os_pT_stage[Os_pT_stage$stages == x,1]) })

Os_pT_stage[grep("UM ", Os_pT_stage$stage),2] = "UM"
Os_pT_mean = sapply(c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "TP"), function(x){ mean(Os_pT_stage[Os_pT_stage$stages == x,1]) })
Os_pT_stage = data.frame(pT = pseudos$pT[ords], stages = as.character(stages[ords]))


save(Os_meanExp,Os_meanExpL, Os_UMIs, Os_g1, Os_g2, Os_Meta, Os_stages, Os_AlleleFracBin, Os_FracMono, Os_TPM, 
    Os_genes,Os_pT_mean,Os_pT_stage,Os_pT_mean2, file="Rice.RData")

pVel2 = function (xx, nboot = 200) 
{
    xx = t(xx)
    stages = as.numeric(stages)
    names(stages) = rownames(xx)
    M1_samples = names(which(stages <= 3))
    G_samples = names(which(stages >= 3))
    set.seed(1)
    out = NULL
    for (i in 1:nboot) {
        samps = sample(rownames(xx), replace = T)
        pT_M1 = pseudotime(xx[samps[samps %in% M1_samples], ], 
            stages)
        pT_G = pseudotime(xx[samps[samps %in% G_samples], ], 
            stages)
        pT = cbind(pT_M1[samps], pT_G[samps])
        pT[, 2] = pT[, 2] - mean(apply(pT, 1, diff), na.rm = T)
        pT = rowMeans(pT, na.rm = T)
        pT = 100 * (pT - quantile(pT, 0.1))/diff(quantile(pT, 
            c(0.1, 0.9)))
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

stages[stages == "TP_shed"] = "TP"
stages = as.character(stages)
names(stages) = over10k
stages[stages == "UM Bi"][Rice_meta[names(stages[stages =="UM Bi"]) ,6] < 1.2] = "UM Bi Early"
stages = factor(stages, levels = c("Tetrad", "UM Bi Early", "UM Bi", "UM Mono", "UM/BM", "BM", "BM/TP", "TP"))
stages[grep("UM ", stages)] = "UM"
stages = factor(stages, levels = c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "TP"))


pseudos2 = pVel2(A2d[rank(-fano) <= 2000,names(stages)], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords2 = order(stages, pseudos2$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime
hmap2 = A2d[,over10k][rank(-fano[rownames(A2d)]) <= 500,ords2]
o1_2 = seriate(dist(t(scale(t(hmap2)))), method = "OLO")


Heatmap(minmax(hmap2), name = 'expression
level', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")), #column_split = split,
    top_annotation = HeatmapAnnotation(Stage = stages2[ords2], "UMIcounts" = log(colSums(D2[,over10k][,ords2]),10),
    "% Monoallelic" = FracMono_all[over10k][ords2], "Anther length" = Rice_meta[over10k,6][ords2],
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM Bi" = "#F9A729", "UM Mix" = "#333", "UM Mono" = "#F9A729",
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D9", "TP_shed" = "gray"), "% Monoallelic" = FracMono_col,
    "UMIcounts" = UMI_col, "Anther length" = Anther_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_rows=as.dendrogram(o1_2[[1]]), cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)



Os_pT_stage2 = data.frame(pT = pseudos2$pT[ords2], stages = as.character(stages[ords2]))
#Os_pT_stage2[grep("TP_shed", Os_pT_stage2$stage),2] = "TP"
Os_pT_mean2 = sapply(c("Tetrad", "UM Bi Early", "UM Bi", "UM Mono", "UM/BM", "BM", "BM/TP", "TP"), function(x){ mean(Os_pT_stage2[Os_pT_stage2$stages == x,1]) })
Os_pT_stage2 = data.frame(pT = pseudos2$pT[ords2], stages = as.character(stages2[ords2]))

> Os_pT_mean2 #2
   Tetrad     UM Bi   UM Mono     UM/BM        BM     BM/TP        TP 
-17.21129  38.98524  16.77699  26.04991  34.83275  60.56024  87.39270 
> Os_pT_mean2 #3
   Tetrad     UM Bi   UM Mono     UM/BM        BM     BM/TP        TP 
-14.84839  36.35049  26.00854  36.74260  45.33976  64.10619  88.82583
> Os_pT_mean2 #4
   Tetrad     UM Bi   UM Mono     UM/BM        BM     BM/TP        TP 
-14.11345  38.95979  27.48823  32.27331  41.47032  61.89538  88.30722
> Os_pT_mean2 #5
   Tetrad     UM Bi   UM Mono     UM/BM        BM     BM/TP        TP 
-13.16329  31.42625  28.87352  39.31211  45.61984  64.74208  89.09802
> Os_pT_mean2 #6
   Tetrad     UM Bi   UM Mono     UM/BM        BM     BM/TP        TP 
-13.52773  34.04528  26.18639  33.62023  42.33483  59.43265  87.53719





###Maize prep
load('MaizeF1_4_2026.RData')
load('MaizeRemap_8-10_2025.RData')

D = cbind(D2,D)
g1 = cbind(g1_2,g1)
g2 = cbind(g2_2,g2)
reads_UMIs = t(cbind(t(reads_UMIs2),t(reads_UMIs)))

Stages <- read_excel("C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R/Maize_Stages.xlsx")
Meta <- Stages[rep(row.names(Stages), each = 8), ]
Meta = as.data.frame(Meta)
Meta$Sample = paste(rep(Stages$Plate, each = 8), '_', rep((as.numeric(sub('[A-Z].+','',Stages$Celseq_primer)) - 1)*8, each = 8) + rep(1:8, nrow(Stages)), 's', sep ='')
rownames(Meta) = Meta$Sample

Meta$No_cell2 = unlist(lapply(Stages$No_cell, function(xx) {       
        out = rep('N', 8)
        out[as.numeric(strsplit(xx,',')[[1]])] = 'Y'
        return(out) }))

D = D[,rownames(Meta)]
g1 = g1[,rownames(Meta)]
g2 = g2[,rownames(Meta)]
reads_UMIs = reads_UMIs[rownames(Meta),]

BINZ = function (xx, bin = 10^7) 
{
    bin = as.numeric(Zm_genes[, 1]) * 10^7 + round(Zm_genes[, 2]/bin)
    out = by(xx, bin, colSums)
    out2 = t(matrix(unlist(out), nrow = ncol(g1)))
    colnames(out2) = colnames(g1)
    rownames(out2) = names(out)
    return(out2)
}

g1_bin = BINZ(g1)
g2_bin = BINZ(g2)
AlleleFrac_bin = g1_bin/(g1_bin + g2_bin)
AlleleFrac_bin[(g1_bin+g2_bin) < 10] = NA #remove bins with <10 genoinformative transcripts
binUse = which(abs(rowMeans(AlleleFrac_bin, na.rm=T) - .5) < .4)  # Exclude bins with >90% of all transcripts mapping to the same allele across all samples
AlleleFrac_bin[-binUse,] = NA
FracMono_all = 100*colMeans(abs(AlleleFrac_bin - .5) >= .3, na.rm=T)
FracMono_all[names(which(colSums(!is.na(AlleleFrac_bin))<10))] = NA

No_cell = rownames(Meta[(which(Meta[,15] == "Y")),])
over10k = names(which(colSums(D) >10000))
ETOH_Remove = names(which(colSums(D[,over10k[grep("ETOH", over10k)]]) < 50000))
over10k = over10k[-which(over10k %in% ETOH_Remove)]
over10k = over10k[-which(over10k %in% No_cell)] #2 no cell made it through
over10k = over10k[-grep("BBB", over10k)] #154

#excluding contamination suspect samples: c("M25-36_82s", "M25-36_86s", "M13-24_11s", "M13-24_30s", "M37-48_ETOH_74s") #maybe "M13-24_57s"
#over10k = over10k[-which(over10k %in% c("M25-36_82s", "M25-36_86s", "M13-24_11s", "M13-24_30s", "M37-48_ETOH_74s"))] #149

stages = as.character(Meta[over10k,13])
names(stages) = over10k
stages = factor(stages, levels = c('UM','UM/BM','BM','BM/TP','TP','TP_shed'))
table(stages)

pseudocount = 1*10^6/quantile(colSums(D), p = .1)
A2 = sweep(D, 2, colSums(D), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(D[,colnames(A2b)] >= 10) >= 10, ]  # Require each gene to have at least 10 UMIs in at least 10 cells

fano = apply(A2d, 1, var)/rowMeans(A2d)  # fano factor is a measure of gene variance
pseudos = pVel(A2d[rank(-fano) <= 2000,over10k], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords = order(stages, pseudos$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime
hmap = A2d[,over10k][rank(-fano[rownames(A2d)]) <= 500,ords]
o1 = seriate(dist(t(scale(t(hmap)))), method = "OLO")

FracMono_col = colorRamp2(c(0, 100), c("white", "purple4"))
UMI_col = colorRamp2(c(4.0, 6.0), c("white", "forestgreen")) #summary(log(colSums(D[,over10k]),10)) #to choose min and max

#manually splitting using table(stages) and heatmaps to guide
split = factor(rep(c("Zm2_UM", "Zm2_UM/BM", "Zm2_BM1", "Zm2_BM2", "Zm2_BM/TP", "Zm2_TP", "Zm2_TP_shed"), c(2,1,6,22,19,76,23)), levels = c("Zm2_UM", "Zm2_UM/BM", "Zm2_BM1", "Zm2_BM2", "Zm2_BM/TP", "Zm2_TP", "Zm2_TP_shed"))

Heatmap(minmax(hmap), column_split = split, 
name = 'expression
level (TPM)', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")),
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    "% Monoallelic" = FracMono_all[over10k][ords], 
    col = list(Stage = c("UM" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", 
    "TP_shed" = "#BBC9E5"), "% Monoallelic" = FracMono_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_rows=as.dendrogram(o1[[1]]), cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)

Heatmap(cor(hmap[,over10k][,ords], use = 'pairwise.complete.obs', method = 'pearson'), 
name = 'Expression_cor', right_annotation = rowAnnotation(Stage = stages[ords],
    "% Monoallelic" = FracMono_all[over10k][ords], 
    col = list(Stage = c("UM" = "#F9A729", "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", 
    "TP_shed" = "#BBC9E5"), "% Monoallelic" = FracMono_col)),
col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100))

splitSamps = data.frame(split = split, samp = names(stages[ords]))
calcMeanExp = function(x){
    rowMeans(A2[,splitSamps[splitSamps[,1] == x,2]]) #mean TPM of each gene for each group
}
calcMeanExpL = function(x){
    rowMeans(log(A2[,splitSamps[splitSamps[,1] == x,2]]+100,10)) #mean log(TPM+100,10) of each gene for each group
}

Zm2_meanExp = sapply(levels(split)[which(table(split)>=3)], calcMeanExp) #apply meanExp only to the groups with >=3 samples
Zm2_meanExpL = sapply(levels(split)[which(table(split)>=3)], calcMeanExpL) #apply meanExpL only to the groups with >=3 samples
Zm2_UMIs = D[,over10k]
Zm2_g1 = g1[,over10k]
Zm2_g2 = g2[,over10k]
Zm2_Meta = Meta[over10k,]
Zm2_stages = stages
Zm2_AlleleFracBin = AlleleFrac_bin[,over10k]
Zm2_FracMono = FracMono_all[over10k]
Zm2_TPM = A2[,over10k]

save(Zm2_meanExp, Zm2_meanExpL, Zm2_UMIs, Zm2_g1, Zm2_g2, Zm2_Meta, Zm2_stages, Zm2_AlleleFracBin, Zm2_FracMono, Zm2_TPM, file="Zm2.RData")

###Loading 2021 data
# Load transcript count data and metadata from GEO accession GSE175505
meta = read.csv('GSE175505_Sample_Metadata.csv.gz', row.names=1)  # This file is also identical to Table S1 of Nelms and Walbot (2021)
A = read.csv('GSE175505_Pollen_scRNAseq_transcript_counts.csv.gz', row.names=1)

# Remove samples that did not pass quality control
A = A[,sub('s_.+','s',colnames(A)) %in% (rownames(meta)[meta[,11] == 'Y'])]  
meta = meta[unique(sub('s_.+','s',colnames(A))),]

# Separate data into 3 matrices: the total number of transcripts per precursor (UMIs), the number of genoinformative transcripts matching the B73 allele (g1), and the number of genoinformative transcripts matching the A188 allele (g2)
UMIs = A[,seq(1,ncol(A),4)] + A[,seq(2,ncol(A),4)] + A[,seq(3,ncol(A),4)] + A[,seq(4,ncol(A),4)]
colnames(UMIs) = sub('_cf', '', colnames(UMIs))
g1 = A[,seq(2,ncol(A),4)]
colnames(g1) = sub('_g1', '', colnames(g1))
g2 = A[,seq(3,ncol(A),4)]
colnames(g2) = sub('_g2', '', colnames(g2))

over10k = names(which(colSums(UMIs) > 10000))
A2 = sweep(UMIs, 2, colSums(UMIs), '/')*10^6  # Transcripts per million normalization

#Downloaded table with all gene model versions: https://download.maizegdb.org/Pan-genes/Pan-Zea/
PanGeneTab = read.delim("pan-zea.v4.pan-genes_table.tsv")
V5toV4tab = PanGeneTab[,6:7]
V5toV4tab = V5toV4tab[which(rowSums(V5toV4tab == "NONE") == 0),]
V5toV4tab$Zm00001d = sub("_.*", "", V5toV4tab$Zm00001d)
V5toV4tab$Zm00001eb = sub("_.*", "", V5toV4tab$Zm00001eb)

vConvert = V5toV4tab[which(V5toV4tab$Zm00001d %in% rownames(UMIs)),]
UMIs = UMIs[vConvert[,1],]
g1 = g1[vConvert[,1],]
g2 = g2[vConvert[,1],]
A2 = A2[vConvert[,1],]
rownames(UMIs) = vConvert[,2]
rownames(g1) = vConvert[,2]
rownames(g2) = vConvert[,2]
rownames(A2) = vConvert[,2]

Zm2_UMIs = Zm2_UMIs[rownames(UMIs),] #goes from 39756 genes to 23804 genes, and ~25% of transcripts are lost in both datasets
Zm2_g1 = Zm2_g1[rownames(UMIs),]
Zm2_g2 = Zm2_g2[rownames(UMIs),]
Zm2_meanExp = Zm2_meanExp[rownames(UMIs),]
Zm2_meanExpL = Zm2_meanExpL[rownames(UMIs),]
Zm2_TPM = Zm2_TPM[rownames(UMIs),]

stages = as.character(meta[colnames(UMIs),4])
names(stages) = colnames(UMIs)
stages = stages[over10k] #applying 10k cutoff (Brad's original cutoff was 2k)
stages[stages == "Pollen"] = "TP"
excStage = c(rownames(meta[meta[colnames(UMIs),3] == "Interphase",]),rownames(meta[meta[colnames(UMIs),3] == "Meiosis",]))
stages = stages[-which(names(stages) %in% excStage)]
stages = factor(stages, levels = c("Tetrad", "Early", "Middle", "Late", "BM", "TP")) 

#Tetrad  Early Middle   Late     BM     TP #2k cutoff gives 223 samples
#    22     89     67     19     11     15
#Tetrad  Early Middle   Late     BM     TP #10k cutoff gives 166 samples
#    20     68     44     10      9     15

UMIs = UMIs[,names(stages)]
g1 = g1[,names(stages)]
g2 = g2[,names(stages)]
A2 = A2[,names(stages)]

BINZ = function (xx, bin = 10^7) 
{
    bin = as.numeric(Zm_genes[rownames(UMIs), 1]) * 10^7 + round(Zm_genes[rownames(UMIs), 2]/bin)
    out = by(xx, bin, colSums)
    out2 = t(matrix(unlist(out), nrow = ncol(g1)))
    colnames(out2) = colnames(g1)
    rownames(out2) = names(out)
    return(out2)
}

g1_bin = BINZ(g1)
g2_bin = BINZ(g2)
AlleleFrac_bin = g1_bin/(g1_bin + g2_bin)
AlleleFrac_bin[(g1_bin+g2_bin) < 10] = NA #remove bins with <10 genoinformative transcripts
binUse = which(abs(rowMeans(AlleleFrac_bin, na.rm=T) - .5) < .4)  # Exclude bins with >90% of all transcripts mapping to the same allele across all samples
AlleleFrac_bin[-binUse,] = NA
FracMono_all = 100*colMeans(abs(AlleleFrac_bin - .5) >= .3, na.rm=T)
FracMono_all[names(which(colSums(!is.na(AlleleFrac_bin))<10))] = NA

pseudocount = 1*10^6/quantile(colSums(UMIs), p = .1)
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(UMIs[,colnames(A2b)] >= 10) >= 10, ]  # Require each gene to have at least 10 UMIs in at least 10 cells
fano = apply(A2d, 1, var)/rowMeans(A2d)  # fano factor is a measure of gene variance
pseudos = pVel(A2d[rank(-fano) <= 2000,], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords = order(stages, pseudos$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime
hmap = A2d[rank(-fano[rownames(A2d)]) <= 500,ords]
o1 = seriate(dist(t(scale(t(hmap)))), method = "OLO")

FracMono_col = colorRamp2(c(0, 100), c("white", "purple4"))

#2k cutoff version
#UMI_col = colorRamp2(c(3.15, 5.75), c("white", "forestgreen")) #summary(log(colSums(UMIs),10)) #to choose min and max
#manually splitting using table(stages) and heatmaps to guide
#split = factor(rep(c("Zm_Tetrad", "Zm_UMe1", "Zm_UMe2", "Zm_UMm1", "Zm_UMm2","Zm_UMm3","Zm_UMl", "Zm_BM", "Zm_TP"), c(22,55,34,29,26,12,19,11,15)), levels = c("Zm_Tetrad", "Zm_UMe1", "Zm_UMe2", "Zm_UMm1", "Zm_UMm2","Zm_UMm3","Zm_UMl", "Zm_BM", "Zm_TP"))

#10k cutoff version
UMI_col = colorRamp2(c(3.8, 5.75), c("white", "forestgreen")) #summary(log(colSums(UMIs),10)) #to choose min and max
#split = factor(rep(c("Zm_Tetrad","Zm_UMe","Zm_UMm","Zm_UMl","Zm_BM","Zm_TP"), c(20,68,44,10,9,15)), levels = c("Zm_Tetrad","Zm_UMe","Zm_UMm","Zm_UMl","Zm_BM","Zm_TP"))
split = factor(rep(c("Zm_Tetrad", "Zm_UMe1", "Zm_UMe2", "Zm_UMm1", "Zm_UMm2","Zm_UMm3","Zm_UMl", "Zm_BM", "Zm_TP"), c(20,53,15,17,21,6,10,9,15)), levels = c("Zm_Tetrad", "Zm_UMe1", "Zm_UMe2", "Zm_UMm1", "Zm_UMm2","Zm_UMm3","Zm_UMl", "Zm_BM", "Zm_TP"))

Heatmap(minmax(hmap), name = 'expression
level (TPM)', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")), column_split = split,
    top_annotation = HeatmapAnnotation(Stage = stages[ords], "UMIcounts" = log(colSums(Maize_UMIs[,ords]),10),
    "% Monoallelic" = FracMono_all[ords], 
    col = list(Stage = c("Tetrad" = "#EB1E2C", "Early" = "#F9A729", "Middle" = "#F9A729", "Late" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", 
    "TP_shed" = "#BBC9E5"), "% Monoallelic" = FracMono_col, "UMIcounts" = UMI_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_rows=as.dendrogram(o1[[1]]), cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)

Heatmap(cor(hmap, use = 'pairwise.complete.obs', method = 'pearson'), 
name = 'Expression_cor', right_annotation = rowAnnotation(Stage = stages[ords], "UMIcounts" = log(colSums(Maize_UMIs[,ords]),10), Source = Source[ords],
    "% Monoallelic" = FracMono_all[ords], 
    col = list(Stage = c("Tetrad" = "#EB1E2C", "Early" = "#F9A729", "Middle" = "#F9D23C", "Late" = "limegreen",
     "BM" = "#5FBB68", "TP" = "#A4A4D5", 
    "TP_shed" = "#BBC9E5"), "% Monoallelic" = FracMono_col, "UMIcounts" = UMI_col, "Source" = c("Brad" = "green2", "Julian" = "purple2"))),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), show_row_names = FALSE, show_column_names = FALSE)

splitSamps = data.frame(split = split, samp = names(stages[ords]))
calcMeanExp = function(x){
    rowMeans(A2[,splitSamps[splitSamps[,1] == x,2]]) #mean TPM of each gene for each group
}
calcMeanExpL = function(x){
    rowMeans(log(A2[,splitSamps[splitSamps[,1] == x,2]]+100,10)) #mean log(TPM+100,10) of each gene for each group
}

###Objects to use downstream in comparisons:
Zm_meanExp = sapply(levels(split)[which(table(split)>=3)], calcMeanExp) #apply meanExp only to the groups with >=3 samples
Zm_meanExpL = sapply(levels(split)[which(table(split)>=3)], calcMeanExpL) #apply meanExpL only to the groups with >=3 samples
Zm_UMIs = UMIs
Zm_g1 = g1
Zm_g2 = g2
Zm_Meta = meta[names(stages),]
Zm_stages = stages
Zm_AlleleFracBin = AlleleFrac_bin
Zm_FracMono = FracMono_all
Zm_TPM = A2

Zm_pT_stage = data.frame(pT = pseudos$pT[ords], stages = as.character(stages[ords]))
Zm_pT_mean = sapply(levels(stages), function(x){ mean(Zm_pT_stage[Zm_pT_stage$stages == x,1]) })

save(Zm_meanExp, Zm_meanExpL, Zm_UMIs, Zm_g1, Zm_g2, Zm_Meta, Zm_stages, Zm_AlleleFracBin, Zm_FracMono, Zm_TPM,
    Zm2_meanExp, Zm2_meanExpL, Zm2_UMIs, Zm2_g1, Zm2_g2, Zm2_Meta, Zm2_stages, Zm2_AlleleFracBin, Zm2_FracMono, Zm2_TPM, file="Maize.RData")


load("Maize.RData")
UMIs = cbind(Zm_UMIs, Zm2_UMIs)

stages = c(as.character(Zm_stages), as.character(Zm2_stages))
names(stages) = c(names(Zm_stages), names(Zm2_stages))
stages[rownames(Zm_Meta[grep("mix", Zm_Meta[,5]),])] = "UM/BM"
stages[which(stages == "Early" | stages == "Middle" | stages == "Late")] = "UM"
stages[which(stages == "TP_shed" | stages == "TP")] = "Pollen"
stages = factor(stages, levels = c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "Pollen"))

A2 = sweep(UMIs, 2, colSums(UMIs), '/')*10^6  # Transcripts per million normalization
pseudocount = 1*10^6/quantile(colSums(UMIs), p = .1)
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(UMIs[,colnames(A2b)] >= 10) >= 10, ]  # Require each gene to have at least 10 UMIs in at least 10 cells
fano = apply(A2d, 1, var)/rowMeans(A2d)  # fano factor is a measure of gene variance
pseudos = pVel(A2d[rank(-fano) <= 2000,], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords = order(stages, pseudos$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime

Zm_merge_pT_stage = data.frame(pT = pseudos$pT[ords], stages = as.character(stages[ords]))
Zm_merge_pT_mean = sapply(levels(stages), function(x){ mean(Zm_merge_pT_stage[Zm_merge_pT_stage$stages == x,1]) })
Zm_merge_UMIs = UMIs
Zm_merge_stages = stages
Zm_merge_g1 = cbind(Zm_g1, Zm2_g1)
Zm_merge_g2 = cbind(Zm_g2, Zm2_g2)
Zm_merge_TPM = cbind(Zm_TPM, Zm2_TPM)
Zm_merge_FracMono = c(Zm_FracMono, Zm2_FracMono)

save(Zm_merge_pT_stage, Zm_merge_pT_mean, Zm_merge_UMIs, Zm_merge_stages, Zm_merge_g1, Zm_merge_g2, Zm_merge_TPM, Zm_merge_FracMono, 
    Zm_AlleleFracBin, Zm2_AlleleFracBin, Zm_Meta, Zm2_Meta, file="Maize_merge.RData")



pVel2 = function (xx, nboot = 200) 
{
    xx = t(xx)
    stages = as.numeric(stages)
    names(stages) = rownames(xx)
    M1_samples = names(which(stages <= 3))
    G_samples = names(which(stages >= 3))
    set.seed(1)
    out = NULL
    for (i in 1:nboot) {
        samps = sample(rownames(xx), replace = T)
        pT_M1 = pseudotime(xx[samps[samps %in% M1_samples], ], 
            stages)
        pT_G = pseudotime(xx[samps[samps %in% G_samples], ], 
            stages)
        pT = cbind(pT_M1[samps], pT_G[samps])
        pT[, 2] = pT[, 2] - mean(apply(pT, 1, diff), na.rm = T)
        pT = rowMeans(pT, na.rm = T)
        pT = 100 * (pT - quantile(pT, 0.1))/diff(quantile(pT, 
            c(0.1, 0.9)))
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

pseudos2 = pVel2(A2d[rank(-fano) <= 2000,names(stages)], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords2 = order(stages, pseudos2$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime

Zm_pT_stage2 = data.frame(pT = pseudos2$pT[ords2], stages = as.character(stages[ords2]))
Zm_pT_mean2 = sapply(levels(stages), function(x){ mean(Zm_pT_stage2[Zm_pT_stage2$stages == x,1]) })



###Tomato prep

load("TomatoRemap_11_2024-7_2025_with_meta.RData")

BINSl = function (xx, bin = 10^6) 
{
    bin = as.numeric(Sl_genes2[, 1]) * 10^6 + round(Sl_genes2[, 2]/bin)
    out = by(xx, bin, colSums)
    out2 = t(matrix(unlist(out), nrow = ncol(g1_2)))
    colnames(out2) = colnames(g1_2)
    rownames(out2) = names(out)
    return(out2)
}

g1_bin = BINSl(g1_2)
g2_bin = BINSl(g2_2)
AlleleFrac_bin = g1_bin/(g1_bin + g2_bin)
AlleleFrac_bin[(g1_bin+g2_bin) < 10] = NA #remove bins with <10 genoinformative transcripts
Sl_binUse = which(abs(rowMeans(AlleleFrac_bin, na.rm=T) - .5) < .4)  # was .4 instead of .3, .3 is >80%. Exclude bins with >90% of all transcripts mapping to the same allele across all samples
AlleleFrac_bin[-Sl_binUse,] = NA
FracMono_all = 100*colMeans(abs(AlleleFrac_bin - .5) >= .3, na.rm=T)
FracMono_all[names(which(colSums(!is.na(AlleleFrac_bin))<10))] = NA

over10k = names(which(colSums(D2) >10000))
over10k = over10k[which(Sl_meta[over10k,14] == "N")] #getting rid of no cell controls

stages = Sl_meta[over10k,8]
names(stages) = over10k
stages = stages[-which(is.na(FracMono_all[over10k]))]
stages = stages[-which(stages == "Meiosis")]
stages[grep("tetrad", stages)] = "Tetrad"
stages[names(which(FracMono_all[names(stages[stages == "UM"])] <= 40))] = "UM Bi"
stages[names(which(FracMono_all[names(stages[stages == "UM"])] > 40))] = "UM Mono"
stages = factor(stages, level = c("Tetrad", "UM Bi", "UM Mono", "UM/BM", "BM", "BM_shed"))

pseudocount = 1*10^6/quantile(colSums(D2), p = .1)
A2 = sweep(D2, 2, colSums(D2), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(D2[,colnames(A2b)] >= 5) >= 3, ]  # Require each gene to have at least 5 UMIs in at least 3 cells

fano = apply(A2d, 1, var)/rowMeans(A2d)  # fano factor is a measure of gene variance
pseudos = pVel(A2d[rank(-fano) <= 2000,names(stages)], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords = order(stages, pseudos$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime
hmap = A2d[,names(stages)][rank(-fano[rownames(A2d)]) <= 500,ords]
o1 = seriate(dist(t(scale(t(hmap)))), method = "OLO")

FracMono_col = colorRamp2(c(0, 100), c("white", "purple4"))
UMI_col = colorRamp2(c(4.0, 5.6), c("white", "forestgreen")) #summary(log(colSums(D2[,over10k]),10)) #to choose min and max

#manually splitting using table(stages) and heatmaps to guide
split = factor(rep(c("Sl_Tetrad", "Sl_UMbi", "Sl_UMmono", "Sl_UM/BM", "Sl_BM1","Sl_BM2", "Sl_BM_shed"), c(19,5,4,10,24,23,31)), levels = c("Sl_Tetrad", "Sl_UMbi", "Sl_UMmono", "Sl_UM/BM", "Sl_BM1","Sl_BM2", "Sl_BM_shed"))

Heatmap(minmax(hmap), column_split = split, 
name = 'expression
level (TPM)', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")),
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    "% Monoallelic" = FracMono_all[names(stages)][ords], 
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM Bi" = "#F9A729", "UM Mono" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM_shed" = "#64CDCC"), "% Monoallelic" = FracMono_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_rows=as.dendrogram(o1[[1]]), cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)

Heatmap(cor(hmap, use = 'pairwise.complete.obs', method = 'pearson'), 
name = 'Expression_cor', right_annotation = rowAnnotation(Stage = stages[ords],
    "% Monoallelic" = FracMono_all[names(stages)][ords], 
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM Bi" = "#F9A729", "UM Mono" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM_shed" = "#64CDCC"), "% Monoallelic" = FracMono_col)),
col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), show_row_names = FALSE, show_column_names = FALSE)

splitSamps = data.frame(split = split, samp = names(stages[ords]))
calcMeanExp = function(x){
    rowMeans(A2[,splitSamps[splitSamps[,1] == x,2]]) #mean TPM of each gene for each group
}
calcMeanExpL = function(x){
    rowMeans(log(A2[,splitSamps[splitSamps[,1] == x,2]]+100,10)) #mean log(TPM+100,10) of each gene for each group
}

###Objects to use downstream in comparisons:
Sl_meanExp = sapply(levels(split)[which(table(split)>=3)], calcMeanExp) #apply meanExp only to the groups with >=3 samples
Sl_meanExpL = sapply(levels(split)[which(table(split)>=3)], calcMeanExpL) #apply meanExp only to the groups with >=3 samples
Sl_UMIs = D2[,names(stages)]
Sl_g1 = g1_2[,names(stages)]
Sl_g2 = g2_2[,names(stages)]
Sl_Meta = Sl_meta[names(stages),]
Sl_stages = stages
Sl_AlleleFracBin = AlleleFrac_bin[,names(stages)]
Sl_FracMono = FracMono_all[names(stages)]
Sl_TPM = A2[,names(stages)]
Sl_genes = Sl_genes2

Sl_pT_stage = data.frame(pT = pseudos$pT[ords], stages = as.character(stages[ords]))
Sl_pT_stage[grep("UM ", Sl_pT_stage$stage),2] = "UM"
Sl_pT_mean = sapply(c("Tetrad", "UM", "UM/BM", "BM", "BM_shed"), function(x){ mean(Sl_pT_stage[Sl_pT_stage$stages == x,1]) })
Sl_pT_stage = data.frame(pT = pseudos$pT[ords], stages = as.character(stages[ords]))
Sl_pT_mean2 = sapply(c("Tetrad", "UM Bi", "UM Mono", "UM/BM", "BM", "BM_shed"), function(x){ mean(Sl_pT_stage[Sl_pT_stage$stages == x,1]) })


save(Sl_meanExp, Sl_meanExpL, Sl_UMIs, Sl_g1, Sl_g2, Sl_Meta, Sl_stages, Sl_AlleleFracBin, Sl_FracMono, Sl_TPM, 
    Sl_genes = Sl_genes2, Sl_pT_stage, Sl_pT_mean, Sl_pT_mean2, file="Tomato.RData")


###Arabidopsis prep
load("At_data_2_2026.RData")

over10k = names(which(colSums(D) >10000))
over10k = over10k[which(At_meta[over10k,12] == "N")] #getting rid of the no cell controls (5 at 10k, 2 at 20k, 0 at 25k)

HapFrac = HapFrac_all[over10k]
HapFrac = (HapFrac - 0.5)*200

stages = At_meta[over10k,9]
names(stages) = over10k
stages = sub("tetrads", "Tetrad", stages)
stages = sub("Tri", "TP", stages)
#stages[grep("shed", At_meta[over10k,11])] = "TP_shed"
stages = factor(stages, levels = c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "TP"))#, "TP_shed"))

pseudocount = 1*10^6/quantile(colSums(D), p = .1) #Set pseudocount as 10^6 / 10% quantile of UMI counts per sample
A2 = sweep(D, 2, colSums(D), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(D[,colnames(A2b)] >= 5) >= 3, ]  # Require each gene to have at least 5 UMIs in at least 3 cells

fano = apply(A2d, 1, var)/rowMeans(A2d)  # fano factor is a measure of gene variance
pseudos = pVel(A2d[rank(-fano) <= 2000,names(stages)], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords = order(stages, pseudos$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime
hmap = A2d[,names(stages)][rank(-fano[rownames(A2d)]) <= 500,ords]
o1 = seriate(dist(t(scale(t(hmap)))), method = "OLO")

FracMono_col = colorRamp2(c(0, 100), c("white", "purple4"))
UMI_col = colorRamp2(c(4.0, 5.6), c("white", "forestgreen")) #summary(log(colSums(D2[,over10k]),10)) #to choose min and max

#manually splitting using table(stages) and heatmaps to guide
#split = factor(rep(c("At_Tetrad", "At_UM1","At_UM2","At_UM3", "At_UM/BM", "At_BM", "At_BM/TP","At_TP1","At_TP2","At_TP3","At_TP4","At_TP_shed"), c(7,50,55,33,28,169,8,8,28,70,29,3)), levels = c("At_Tetrad", "At_UM1","At_UM2","At_UM3", "At_UM/BM", "At_BM", "At_BM/TP","At_TP1","At_TP2","At_TP3","At_TP4","At_TP_shed"))

Heatmap(minmax(hmap), #column_split = split,
name = 'expression
level (TPM)', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")),
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    "% Matching Haplotype" = HapFrac[ords],
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM" = "#F9A729", 
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", "TP_shed" = "#BBC9E5"), 
    "% Matching Haplotype" = FracMono_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), 
    cluster_rows=as.dendrogram(o1[[1]]), cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)

Heatmap(cor(hmap, use = 'pairwise.complete.obs', method = 'pearson'), 
name = 'Expression_cor', right_annotation = rowAnnotation(Stage = stages[ords],
    "% Matching Haplotype" = HapFrac[ords],
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM" = "#F9A729", 
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", "TP_shed" = "#BBC9E5"), 
    "% Matching Haplotype" = FracMono_col)),
col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), show_row_names = FALSE, show_column_names = FALSE)

#splitSamps = data.frame(split = split, samp = names(stages[ords]))
#calcMeanExp = function(x){
#    rowMeans(A2[,splitSamps[splitSamps[,1] == x,2]]) #mean TPM of each gene for each group
#}
#calcMeanExpL = function(x){
#    rowMeans(log(A2[,splitSamps[splitSamps[,1] == x,2]]+100,10)) #mean log(TPM+100,10) of each gene for each group
#}

###Objects to use downstream in comparisons:
#At_meanExp = sapply(levels(split)[which(table(split)>=3)], calcMeanExp) #apply meanExp only to the groups with >=3 samples
#At_meanExpL = sapply(levels(split)[which(table(split)>=3)], calcMeanExpL) #apply meanExp only to the groups with >=3 samples
At_UMIs = D[,names(stages)]
At_g1 = g1[,names(stages)]
At_g2 = g2[,names(stages)]
At_Meta = At_meta[names(stages),]
At_stages = stages
At_HapFrac_genes = HapFrac_genes_all[,names(stages)]
At_HapFrac = HapFrac
At_TPM = A2[,names(stages)]


At_pT_stage = data.frame(pT = pseudos$pT[ords], stages = as.character(stages[ords]))
At_pT_stage[grep("TP_shed", At_pT_stage$stage),2] = "TP"
At_pT_mean = sapply(c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "TP"), function(x){ mean(At_pT_stage[At_pT_stage$stages == x,1]) })
At_pT_stage = data.frame(pT = pseudos$pT[ords], stages = as.character(stages[ords]))

save(At_UMIs, At_g1, At_g2, At_Meta, At_stages, At_HapFrac_genes, At_HapFrac, At_TPM, file="Arabidopsis.RData")
#Need to add in At_4_2026.RData data, which has 12 samples over10k (need to calculate HapFrac for these)
save(At_UMIs, At_g1, At_g2, At_Meta, At_stages, At_HapFrac_genes, At_HapFrac, At_TPM, file="Arabidopsis10k.RData")


source('Pseudotime Velocity Functions.R')
pVel2 = function(xx, nboot = 200)
{
    xx = t(xx)

    stages = as.numeric(stages)
    names(stages) = rownames(xx)

    M1_samples = names(which(stages <= 2))
    G_samples  = names(which(stages >= 2))

    set.seed(1)
    out = NULL

    for (i in 1:nboot) {

        # Stratified bootstrap: preserve representation of each stage
        samps = unlist(
            lapply(unique(stages), function(s) {
                ss = names(stages[stages == s])
                sample(ss, length(ss), replace = TRUE)
            })
        )

        M1_boot = samps[samps %in% M1_samples]
        G_boot  = samps[samps %in% G_samples]

        pT_M1 = pseudotime(
            xx[M1_boot, ],
            stages
        )

        pT_G = pseudotime(
            xx[G_boot, ],
            stages
        )

        pT = cbind(
            pT_M1[samps],
            pT_G[samps]
        )

        pT[, 2] = pT[, 2] -
            mean(apply(pT, 1, diff), na.rm = TRUE)

        pT = rowMeans(pT, na.rm = TRUE)

        q = quantile(
            pT,
            c(0.1, 0.9),
            na.rm = TRUE
        )

        if (diff(q) == 0 || !all(is.finite(q))) {
            next
        }

        pT = 100 * (pT - q[1]) / diff(q)

        names(pT) = samps
        out = c(out, pT)

        if (i == 1) {
            cat("0% |---------|---------|---------|---------| 100%\n   .")
        }

        if ((i %% round(nboot/40)) == 0) {
            cat(".")
        }
    }

    pT = by(out, names(out), mean)[rownames(xx)]

    pV = c(
        rep(NA, 4),
        velFn(pT[order(pT)]),
        rep(NA, 4)
    )

    pV = pV / median(pV, na.rm = TRUE)

    return(list(pT = pT, pV = pV))
}

pseudos2 = pVel2(A2d[rank(-fano) <= 2000,names(stages)], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords2 = order(stages, pseudos2$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime

At_pT_stage2 = data.frame(pT = pseudos2$pT[ords2], stages = as.character(stages[ords2]))
At_pT_stage2[grep("TP_shed", At_pT_stage2$stage),2] = "TP"
At_pT_mean2 = sapply(c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "TP"), function(x){ mean(At_pT_stage2[At_pT_stage2$stages == x,1]) })
At_pT_stage2 = data.frame(pT = pseudos2$pT[ords2], stages = as.character(stages[ords2]))

save(Os_pT_stage, Os_pT_mean, Os_pT_mean2, Zm_pT_stage2, Zm_pT_mean2, Sl_pT_stage, Sl_pT_mean, Sl_pT_mean2, At_pT_stage, At_pT_mean, At_pT_stage2, At_pT_mean2, file="Pseudotime_comparison.RData")
#use Os_pT_mean, Zm_pT_mean2, Sl_pT_mean, At_pT_mean2

save(Os_pT_stage, Zm_merge_pT_stage, Sl_pT_stage, At_pT_stage2, file="Pseudotime_comparison2.RData")



pVel2 = function (xx, nboot = 200) 
{
    xx = t(xx)
    stages = as.numeric(stages)
    names(stages) = rownames(xx)
    M1_samples = names(which(stages <= 2))
    G_samples = names(which(stages >= 2))
    set.seed(1)
    out = NULL
    for (i in 1:nboot) {
        samps = sample(rownames(xx), replace = T)
        pT_M1 = pseudotime(xx[samps[samps %in% M1_samples], ], 
            stages)
        pT_G = pseudotime(xx[samps[samps %in% G_samples], ], 
            stages)
        pT = cbind(pT_M1[samps], pT_G[samps])
        pT[, 2] = pT[, 2] - mean(apply(pT, 1, diff), na.rm = T)
        pT = rowMeans(pT, na.rm = T)
        pT = 100 * (pT - quantile(pT, 0.1))/diff(quantile(pT, 
            c(0.1, 0.9)))
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

pseudos2 = pVel2(A2d[rank(-fano) <= 2000,names(stages)], nboot = 80)  # Calculate pseudotime and pseudotime velocity using accessory functions in R script
ords2 = order(stages, pseudos2$pT)  # order samples based on morphologically-defined stage, then order samples within a stage by pseudotime

At_pT_stage2 = data.frame(pT = pseudos2$pT[ords2], stages = as.character(stages[ords2]))
At_pT_stage2[grep("TP_shed", At_pT_stage2$stage),2] = "TP"
At_pT_mean2 = sapply(c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "TP"), function(x){ mean(At_pT_stage2[At_pT_stage2$stages == x,1]) })
At_pT_stage2 = data.frame(pT = pseudos2$pT[ords2], stages = as.character(stages[ords2]))

save(Os_pT_stage, Zm_merge_pT_stage, Sl_pT_stage, At_pT_stage2, file="Pseudotime_comparison2_At10k.RData")



load('At_4_2026.RData')

library(readxl)
At_Stages <- read_excel("C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R/At_Stages_4-21-26.xlsx")

At_meta <- At_Stages[rep(row.names(At_Stages), each = 8), ]
At_meta = as.data.frame(At_meta)

At_meta$Sample = c(paste(rep('A339-350_', times = 96), 1:96, rep('s', times=96), sep=''), paste(rep('A351-362_', times = 96), 1:96, rep('s', times=96), sep=''))
rownames(At_meta) = At_meta$Sample

At_meta$No_cell2 = unlist(lapply(At_Stages$No_Cell, function(xx) {       
        out = rep('N', 8)
        out[as.numeric(strsplit(xx,',')[[1]])] = 'Y'
        return(out) }))

NCCs = At_meta[which(At_meta$No_cell2 == 'Y'),10]
colSums(D[,NCCs])

over10k = names(which(colSums(D) >10000)) #12 samples
over10k = over10k[which(At_meta[over10k,11] == "N")] #getting rid of no cell controls (no NCCs had passed 10k)


