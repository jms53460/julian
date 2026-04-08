Arabidopsis.r
###In local R terminal
setwd('C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R')
load("At_data_2_2026.RData")

library('ComplexHeatmap')
library(circlize)

library(ggplot2)
library(ggpubr)

col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))
#Heatmap(HapExpScore, col = col_fun, cluster_columns=T, cluster_rows=T, show_row_names = FALSE, show_column_names = FALSE)

pseudocount = 1*10^6/quantile(colSums(D), p = .1) #Set pseudocount as 10^6 / 10% quantile of UMI counts per sample
A2 = sweep(D, 2, colSums(D), '/')*10^6  # Transcripts per million normalization
A2b = log(A2+pseudocount,10)  # Log transform
A2d = A2b[rowSums(D[,colnames(A2b)] >= 5) >= 3, ]  # Require each gene to have at least 5 UMIs in at least 3 cells
over30k = names(which(colSums(D) >20000))
over30k = over30k[which(At_meta[over30k,12] == "N")] #getting rid of the 5 no cell controls

fano = apply(A2d, 1, var)/rowMeans(A2d)  # fano factor is a measure of gene variance
hmat = A2d[rank(-fano[rownames(A2d)]) <= 500,over30k]
minmax = function(x) {
	sweep(x - log(pseudocount, 10), 1, apply(x - log(pseudocount, 10), 1, max), '/')
}


stages = At_meta[over30k,9]
names(stages) = over30k
stages = sub("tetrads", "Tetrad", stages)
stages = sub("Tri", "TP", stages)
stages = factor(stages, levels = c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "TP"))


stages = At_meta[over30k,9]
names(stages) = over30k
stages = sub("tetrads", "Tetrad", stages)
stages = sub("Tri", "TP", stages)
Bi = names(which(HapFrac_all[names(stages)] <0.6))
Mono = names(which(HapFrac_all[names(stages)] >=0.6))
stages[names(which(stages[Bi] == "UM"))] = "UM Bi"
stages[names(which(stages[Mono] == "UM"))] = "UM Mono"
stages[names(which(stages[Bi] == "UM/BM"))] = "UM/BM Bi"
stages[names(which(stages[Mono] == "UM/BM"))] = "UM/BM Mono"
stages[names(which(stages[Bi] == "BM"))] = "BM Bi"
stages[names(which(stages[Mono] == "BM"))] = "BM Mono"
stages[names(which(stages[Bi] == "TP"))] = "TP Bi"
stages[names(which(stages[Mono] == "TP"))] = "TP Mono"
stages = factor(stages, levels = c("Tetrad", "UM Bi", "UM Mono", "UM/BM Bi", "UM/BM Mono", "BM Bi", "BM Mono", "BM/TP", "TP Bi", "TP Mono"))


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

HapFrac = HapFrac_all[over30k]
HapFrac = (HapFrac - 0.5)*200
FracMono_col = colorRamp2(c(0, 100), c("white", "purple4"))


CellCycleGs = c('ID=gene-AT3G48750', 'ID=gene-AT3G54180', 'ID=gene-AT2G38620', 'ID=gene-AT1G76540', 'ID=gene-AT1G20930', 'ID=gene-AT5G10270', 'ID=gene-AT5G64960', 'ID=gene-AT1G73690', 'ID=gene-AT1G66750', 'ID=gene-AT1G18040', 'ID=gene-AT5G63610', 
'ID=gene-AT4G28980', 'ID=gene-AT1G44110','ID=gene-AT1G77390','ID=gene-AT5G25380','ID=gene-AT5G11300','ID=gene-AT1G15570','ID=gene-AT1G80370','ID=gene-AT5G43080','ID=gene-AT1G47210','ID=gene-AT1G47220','ID=gene-AT1G47230','ID=gene-AT4G37490',
'ID=gene-AT5G06150','ID=gene-AT3G11520','ID=gene-AT2G26760','ID=gene-AT2G17620','ID=gene-AT4G35620','ID=gene-AT1G20610','ID=gene-AT1G76310','ID=gene-AT1G16330','ID=gene-AT1G70210','ID=gene-AT2G22490','ID=gene-AT4G34160','ID=gene-AT5G67260',
'ID=gene-AT3G50070','ID=gene-AT5G65420','ID=gene-AT5G10440','ID=gene-AT4G37630','ID=gene-AT4G03270','ID=gene-AT5G02110','ID=gene-AT5G27620','ID=gene-AT2G27960','ID=gene-AT2G27970','ID=gene-AT3G48160','ID=gene-AT5G14960','ID=gene-AT3G01330',
'ID=gene-AT5G02470','ID=gene-AT5G03410','ID=gene-AT2G36010','ID=gene-AT5G22220','ID=gene-AT1G47870','ID=gene-AT2G23430','ID=gene-AT3G50630','ID=gene-AT5G48820','ID=gene-AT2G32710','ID=gene-AT3G24810','ID=gene-AT3G19150','ID=gene-AT1G49620',
'ID=gene-AT3G12280','ID=gene-AT1G02970')
CellCycleGs = CellCycleGs[CellCycleGs %in% rownames(A2b)]

hmat3 = A2b[,over30k][CellCycleGs,ords]

DNGs = c('ID=gene-AT5G04560', 'ID=gene-AT2G36490', 'ID=gene-AT3G10010', 'ID=gene-AT5G60250', 'ID=gene-AT5G11800')
hmat4 = A2b[,over30k][DNGs,ords]

E2FGs = c('ID=gene-AT2G36010', 'ID=gene-AT5G22220', 'ID=gene-AT1G47870')
hmat5 = A2b[,over30k][E2FGs,ords]

Heatmap(minmax(hmat5), name = 'expression
level', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")),
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    "% Matching Haplotype" = HapFrac[ords],
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM" = "#F9A729", 
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5"), 
    "% Matching Haplotype" = FracMono_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), cluster_columns=F, show_row_names = T, show_column_names = FALSE)


Heatmap(minmax(hmat2), name = 'expression
level', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")),
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    "% Matching Haplotype" = HapFrac[ords],
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM Bi" = "#F9A729", "UM Mono" = "#F9A729",
    "UM/BM Bi" = "#F9D23C", "UM/BM Mono" = "#F9D23C", "BM Bi" = "#5FBB68", "BM Mono" = "#5FBB68", 
    "BM/TP" = "#64CDCC", "TP Bi" = "#A4A4D5", "TP Mono" = "#A4A4D5"), 
    "% Matching Haplotype" = FracMono_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), 
    cluster_rows=as.dendrogram(o1[[1]]), cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)


#svg('At_expression_2_6_26.svg', width=7, height=4)
Heatmap(minmax(hmat2), name = 'expression
level', heatmap_legend_param = list(at = c(0, 1), labels = c(0, "max")),
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    "% Matching Haplotype" = HapFrac[ords],
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM" = "#F9A729", 
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5"), 
    "% Matching Haplotype" = FracMono_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), 
    cluster_rows=as.dendrogram(o1[[1]]), cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)
#dev.off()


No_cell = At_meta$Sample[which(At_meta[,12] == "Y")]
No_cell_over5k = names(which(colSums(D[,No_cell]) > 5000))


over30k = names(which(colSums(D) >20000))
over30k = over30k[which(At_meta[over30k,12] == "N")] #getting rid of the 5 no cell controls

stages = At_meta[over30k,9]
names(stages) = over30k
stages = stages[-which(stages == "tetrads")]
stages = stages[-which(stages == "BM/Tri")]
stages = sub("Tri", "TP", stages)
stages = factor(stages, levels = c("UM", "UM/BM", "BM", "TP"))
HapFracByStage = data.frame(HapFrac = HapFrac[names(stages)], Stage = stages)

svg('At_HapFracByStage.svg', width=3, height=2)
ggplot(HapFracByStage, aes(x = Stage, y = HapFrac, fill = Stage)) + geom_violin(size=0.5, adjust = 0.3) + labs(title = "Arabidopsis", 
x = "Stage", y = "% Transcripts
    Matching Haplotype") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none", 
    panel.grid=element_blank(), axis.ticks.x = element_blank()) + scale_fill_manual(values=c("UM" = "#F9A729", 
    "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "TP" = "#A4A4D5")) + ylim(0,100)
dev.off()


CellHist = function(stg) {
	gghistogram(HapFracByStage[HapFracByStage$Stage == stg,], 'HapFrac', fill = 'dark gray', bins = 21) + scale_y_continuous(expand = c(0,0)) + xlim(c(-5,105)) + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = 'none', plot.margin = margin(20,0,0,0))
}

svg('At_HapFracByStage_Hist.svg', width=2.5, height=4.5)
ggarrange(CellHist('UM'), CellHist('UM/BM'), CellHist('BM'), CellHist('TP'), nrow = 4, align = 'v')
dev.off()




stages = At_meta[over30k,9]
names(stages) = over30k
stages = sub("tetrads", "Tetrad", stages)
stages = sub("Tri", "TP", stages)
stages = factor(stages, levels = c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "TP"))
BiUM = names(which(HapFrac_all[names(stages[which(stages == 'UM')])] <0.6))
MonoTP = names(which(HapFrac_all[names(stages[which(stages == 'TP')])] >0.9))
HapFrac = HapFrac_genes_all
HapFrac[(g1[,colnames(HapFrac)]+g2[,colnames(HapFrac)]) < 5] = NA
HapFrac = HapFrac*100

#colSums(!is.na(HapFrac[,BiUM]))

ratioHist = function(cell){
gghistogram(data.frame(geneFrac = HapFrac[,cell]), x = 'geneFrac', fill = 'black', bins = 15) + theme(strip.background = element_blank(), strip.text.x = element_blank(), legend.position = 'none') + scale_y_continuous(expand = c(0,0)) + xlab('% Transcripts matching haplotype') + ylab('Genes') + scale_x_continuous(breaks = seq(0,100,10), labels=c('0', rep('',4), '50', rep('',4), '100')) 
}

svg('Arabidopsis_ratioHist_BiUM.svg', width=3.5, height=2)
ratioHist('A182-193_18s') #Biallelic UM
dev.off()

svg('Arabidopsis_ratioHist_MonoTP.svg', width=3.5, height=2)
ratioHist('A242-253_37s') #Monoallelic TP
dev.off()




ratio = g1/(g1+g2)
ratio[(g1+g2) < 10] = NA  # Remove measurements with under 10 genoinformative transcripts
geneUse = which(abs(rowMeans(g1/(g1+g2), na.rm=T) - .5) < .4)  # Exclude genes with >90% of all transcripts mapping to the same allele across all samples
ratio[-geneUse,] = NA
geneFrac = ratio*100

ratioHist = function(cell){
gghistogram(data.frame(geneFrac = geneFrac[,cell]), x = 'geneFrac', fill = 'black', bins = 15) + theme(strip.background = element_blank(), strip.text.x = element_blank(), legend.position = 'none') + scale_y_continuous(expand = c(0,0)) + xlab('% Transcripts matching Col-0') + ylab('Genes') + scale_x_continuous(breaks = seq(0,100,10), labels=c('0', rep('',4), '50', rep('',4), '100')) 
}
svg('Arabidopsis_ratioHist0_BiUM.svg', width=3.2, height=2)
ratioHist('A182-193_18s') #Biallelic UM
dev.off()

svg('Arabidopsis_ratioHist0_MonoTP.svg', width=3.2, height=2)
ratioHist('A242-253_37s') #Monoallelic TP
dev.off()


plotScaleBar = ggplot() + scale_x_continuous(expand=c(0,0), limits = c(-3, max(genes[,2])/10^6 + 3), breaks = seq(0,30,5)) + theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.x=element_line(), plot.margin = margin(0,0,0,0,'cm')) + xlab('Chromosome position (Mb)')

plotChr = function (cell, chr = 1, pad = 3) 
{
    Cdat = data.frame(f_col0 = ratio[, cell], Chr = genes[,1], 
        Position = genes[,2]/10^6)
    Cdat = Cdat[Cdat$Chr == chr, ]
    ggplot(Cdat) + geom_rect(data = data.frame(xmin = -pad, xmax = max(genes[,2])/10^6 + pad, 
        ymin = 0, ymax = 1), aes(xmin = xmin, xmax = xmax, 
        ymin = ymin, ymax = ymax), fill = "#EEEEEE") + geom_point(aes(y = f_col0, 
        x = Position), cex = 1) + geom_hline(yintercept = 0.5, 
        linetype = "dashed") + theme(panel.background = element_blank(), 
        axis.title = element_blank(), panel.border = element_blank(), 
        panel.grid = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) + scale_y_continuous(breaks = seq(0, 
        1, 0.25), labels = c("0%", "", "50%", "", "100%"), limits = c(-0.4, 
        1.05)) + scale_x_continuous(expand = c(0, 0)) + annotate("segment", 
        x = -pad, xend = -pad, y = 0, yend = 1) + theme(plot.margin = margin(0, 
        0, 0.15, 0, "cm"))
}

plotCell = function (cell) 
{
    annotate_figure(ggarrange(plotChr(cell, chr = 1), plotChr(cell, 
        chr = 2), plotChr(cell, chr = 3), plotScaleBar, ncol = 1, nrow = 4, align = "v", 
        heights = c(rep(1, 12), 0.4)), left = text_grob("          % Transcripts from Col-0 allele", 
        rot = 90, size = 10), top = cell)
}

svg('Arabidopsis_plotCell0_BiUM.svg', width=3, height=3)
plotCell('A182-193_18s')
dev.off()

svg('Arabidopsis_plotCell0_MonoTP.svg', width=3, height=3)
plotCell('A242-253_37s')
dev.off()





stages = At_meta[over30k,9]
names(stages) = over30k
stages = sub("tetrads", "Tetrad", stages)
stages = sub("Tri", "TP", stages)
Bi = names(which(HapFrac_all[names(stages)] <0.6))
Mono = names(which(HapFrac_all[names(stages)] >=0.6))
stages[names(which(stages[Bi] == "UM"))] = "UM Bi"
stages[names(which(stages[Mono] == "UM"))] = "UM Mono"
stages[names(which(stages[Bi] == "UM/BM"))] = "UM/BM Bi"
stages[names(which(stages[Mono] == "UM/BM"))] = "UM/BM Mono"
stages[names(which(stages[Bi] == "BM"))] = "BM Bi"
stages[names(which(stages[Mono] == "BM"))] = "BM Mono"
stages[names(which(stages[Bi] == "TP"))] = "TP Bi"
stages[names(which(stages[Mono] == "TP"))] = "TP Mono"
stages = factor(stages, levels = c("Tetrad", "UM Bi", "UM Mono", "UM/BM Bi", "UM/BM Mono", "BM Bi", "BM Mono", "BM/TP", "TP Bi", "TP Mono"))


UM_Mono = names(stages[which(stages == "UM Mono")])
UM_Bi = names(stages[which(stages == "UM Bi")])
UMBM_Mono = names(stages[which(stages == "UM/BM Mono")])
BM_Mono = names(stages[which(stages == "BM Mono")])
BM_Bi = names(stages[which(stages == "BM Bi")])
TP_Mono = names(stages[which(stages == "TP Mono")])
TP_Bi = names(stages[which(stages == "TP Bi")])


summary(colMeans(Hap_all[,over30k], na.rm=T)) #Mean haplotype per sample going from -1 (Ler) to 1 (Col)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.703061 -0.173490 -0.016565 -0.001059  0.153133  0.717311
summary(colMeans(Hap_all[,UM_Mono], na.rm=T))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.68435 -0.20140 -0.03755 -0.01826  0.21395  0.68291

hist(((colMeans(Hap_all[,over30k], na.rm=T))+1)/2)
(summary(colMeans(Hap_all[,over30k], na.rm=T))+1)/2 #Mean haplotype per sample going from 0 (Ler) to 1 (Col)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1485  0.4133  0.4917  0.4995  0.5766  0.8587
hist(((colMeans(Hap_all[,UM_Mono], na.rm=T))+1)/2)
(summary(colMeans(Hap_all[,UM_Mono], na.rm=T))+1)/2
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1578  0.3993  0.4812  0.4909  0.6070  0.8415




Hap_all2 = (Hap_all+1)/2

chr_lengths <- table(genes[,1])[1:5]
chr_ends <- cumsum(chr_lengths)
chr_starts <- c(1, head(chr_ends, -1) + 1)

chr_ranges <- Map(seq, chr_starts, chr_ends)
names(chr_ranges) <- names(chr_lengths)

FillHapGaps_all <- function(Hap_all2, chr_ranges) {
  
  for (chr in names(chr_ranges)) {
    
    rows <- chr_ranges[[chr]]
    chr_mat <- Hap_all2[rows, , drop = FALSE]
    
    for (s in seq_len(ncol(chr_mat))) {
      
      hap_col <- chr_mat[, s]
      known <- which(!is.na(hap_col))
      
      if (length(known) < 2) next #if there are 0 or 1 non NA values, no filling will occur
      
      for (i in seq_len(length(known) - 1)) {
        
        start <- known[i]
        end   <- known[i + 1]
        
        if (end - start > 1) {
          hap_col[(start + 1):(end - 1)] <-
            (hap_col[start] + hap_col[end]) / 2
        }
      }
      
      chr_mat[, s] <- hap_col
    }
    
    Hap_all2[rows, ] <- chr_mat
  }
  
  Hap_all2
}

a1 = proc.time()
Hap_all2 <- FillHapGaps_all(Hap_all2, chr_ranges)
a2 = proc.time()
a2-a1
#   user  system elapsed
#   6.28    0.56    6.91

#FillHapGaps = function(samp, chr) {
#    HapInfoGs = which(!is.na(Hap_all2[which(genes[,1] == chr),samp]))
#    Gs = 1:length(which(genes[,1] == chr))
#    for (x in 1:(length(HapInfoGs)-1)) {
#    Hap_all2[which(genes[,1] == chr),][Gs[Gs > HapInfoGs[x] & Gs < HapInfoGs[x+1]],samp] = mean(c(Hap_all2[which(genes[,1] == chr),][HapInfoGs[x],samp],Hap_all2[which(genes[,1] == chr),][HapInfoGs[x+1],samp]))
#    }
#    return(Hap_all2)
#}

#Hap_all2 = (Hap_all+1)/2
#a1 = proc.time()
#for (samp in UM_Mono){
#    Hap_all2 = FillHapGaps(samp, chr=1)
#    Hap_all2 = FillHapGaps(samp, chr=2)
#    Hap_all2 = FillHapGaps(samp, chr=3)
#    Hap_all2 = FillHapGaps(samp, chr=4)
#    Hap_all2 = FillHapGaps(samp, chr=5)
#}
#a2 = proc.time()
#a2-a1
#   user  system elapsed
#4086.95 1382.53 5560.45 #93 mins for 14 UM_Mono samples, expect 1549 mins or 25.8 hrs for all 234 samples w/ >20k UMIs

Hap_all3 = Hap_all2
Hap_all3[which(Hap_all3 <= 0.2)] = 0
Hap_all3[which(Hap_all3 >= 0.8)] = 1
Hap_all3[Hap_all3 > 0.2 & Hap_all3 < 0.8] = NA

save(Hap_all2, Hap_all3, UM_Mono, UM_Bi, BM_Mono, BM_Bi, TP_Mono, TP_Bi, file='At_hap_update.Rdata')

load('At_hap_update.Rdata')

Hap_all2[rowMeans(!is.na(Hap_all2)) >= 0.9,] #Require 90% of samples to have data for the gene to use it

samp = sample(over30k, 14)
plot(rowMeans(Hap_all2[,samp], na.rm=T), ylim=c(0,1), cex=2, pch=19)
lines(rowMeans(Hap_all2[,samp], na.rm=T), ylim=c(0,1), col='purple2')
abline(v=c(0,chr_ends)+0.5)

samp = sample(TP_Mono, 14)
plot(rowMeans(Hap_all2[,samp], na.rm=T), ylim=c(0,1), cex=2, pch=19)
lines(rowMeans(Hap_all2[,samp], na.rm=T), ylim=c(0,1), col='purple2')
abline(v=c(0,chr_ends)+0.5)


plot(rowMeans(Hap_all2[,over30k], na.rm=T), ylim=c(0,1), cex=2, pch=19)
lines(rowMeans(Hap_all2[,over30k], na.rm=T), ylim=c(0,1), col='purple2')
abline(v=c(0,chr_ends)+0.5)


plot(rowMeans(Hap_all2[,UM_Mono], na.rm=T), ylim=c(0,1), cex=2, pch=19)
lines(rowMeans(Hap_all2[,UM_Mono], na.rm=T), ylim=c(0,1), col='purple2')
abline(v=c(0,chr_ends)+0.5)

plot(rowMeans(Hap_all3[rowMeans(!is.na(Hap_all2)) >= 0.9,over30k], na.rm=T), ylim=c(0,1), cex=2, pch=19)

plot(rowMeans(Hap_all3[rowMeans(!is.na(Hap_all2)) >= 0.9,UM_Mono], na.rm=T), ylim=c(0,1), cex=2, pch=19)
lines(rowMeans(Hap_all3[rowMeans(!is.na(Hap_all2)) >= 0.9,UM_Mono], na.rm=T), ylim=c(0,1), col='purple2')
abline(v=c(0,chr_ends)+0.5)

BootUM_Mono = function(x){
    samps = sample(UM_Mono, replace=TRUE)
    MeanHap = rowMeans(Hap_all3[rowMeans(!is.na(Hap_all2)) >= 0.9,samps], na.rm=T)
    return(MeanHap)
}

set.seed(1)
a1 = proc.time()
BootsUM_Mono = sapply(1:100, BootUM_Mono)
a2 = proc.time()
a2-a1
#   user  system elapsed    100 bootstraps. For 10,000 would need 1.371 hours, for 100,000 would need 13.71 hours
#  29.34   12.69   49.36

set.seed(1)
a1 = proc.time()
BootsUM_Mono = sapply(1:10000, BootUM_Mono)
a2 = proc.time()
a2-a1
#   user  system elapsed    10,000 bootstraps took 1.1056 hours
#2605.29  942.41 3980.17

Frac95CI = apply(BootsUM_Mono, 1, quantile, probs=c(0.025,0.975), na.rm=T)
plot(Frac95CI['2.5%',], Frac95CI['97.5%',], cex=2, pch=19)

plot(rowMeans(BootsUM_Mono), ylim=c(0,1), cex=2, pch=19)
points(Frac95CI['2.5%',], col='blue2', cex=1, pch=19)
points(Frac95CI['97.5%',], col='red2', cex=1, pch=19)

#abline(v=c(0,chr_ends)+0.5) #currently inaccurate for this plot

save(BootsUM_Mono, Frac95CI, file = "10kBootsUM_Mono_At.Rdata")


over5k = names(which(colSums(D) > 5000))
stages = At_meta[over5k,9]
names(stages) = over5k
stages = sub("tetrads", "Tetrad", stages)
stages = sub("Tri", "TP", stages)
Bi = names(which(HapFrac_all[names(stages)] <0.7))
Mono = names(which(HapFrac_all[names(stages)] >=0.7))
stages[names(which(stages[Bi] == "UM"))] = "UM Bi"
stages[names(which(stages[Mono] == "UM"))] = "UM Mono"
stages[names(which(stages[Bi] == "UM/BM"))] = "UM/BM Bi"
stages[names(which(stages[Mono] == "UM/BM"))] = "UM/BM Mono"
stages[names(which(stages[Bi] == "BM"))] = "BM Bi"
stages[names(which(stages[Mono] == "BM"))] = "BM Mono"
stages[names(which(stages[Bi] == "TP"))] = "TP Bi"
stages[names(which(stages[Mono] == "TP"))] = "TP Mono"
stages = factor(stages, levels = c("Tetrad", "UM Bi", "UM Mono", "UM/BM Bi", "UM/BM Mono", "BM Bi", "BM Mono", "BM/TP", "TP Bi", "TP Mono"))


UM_Mono = names(stages[which(stages == "UM Mono")])


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
g1_frac = g1_bin/(g1_bin + g2_bin)
AlleleFrac = g1_frac
AlleleFrac[(g1_bin+g2_bin) < 10] = NA #remove bins with <10 genoinformative transcripts

#Check if the haplotype calls make sense for the monoallelic UMs
plotCellChrHapVsBin = function(cell, chr, maxs = 20) {
	keeps = ((g1[,cell] + g2[,cell]) > 0) & (genes$Chr == chr)
    keeps2 = which(round(as.numeric(rownames(AlleleFrac))/(10^6)) == chr)
    keeps3 = which(!is.na(AlleleFrac[keeps2,cell]))
    keeps4 = rownames(genes[genes$Chr == chr,])

    g1a = g1[keeps,cell]
	g1a[g1a > maxs] = maxs
    g1b = g1[,cell]
    g1b[g1b > maxs] = maxs
	g2a = g2[keeps,cell]
	g2a[g2a > maxs] = maxs
    g2b = g2[,cell]
    g2b[g2b > maxs] = maxs

    par(mar = c(5.1, 5.1, 4.1, 4.1))
	plot(genes$Position[keeps], g1a, pch = 19, cex=0.5, col = 'red', ylab=NA, xlab = 'chromosome position (Mb)', ylim = c(-maxs,maxs), xaxt = 'n', yaxt = 'n', main = cell)
	abline(h=0)
	for (ii in which(keeps)) { lines(rep(genes$Position[ii], times=2), c(g1b[ii], -g2b[ii])) }
	points(genes$Position[keeps], -g2a, pch = 19, cex=0.5, col = 'blue')
	points(genes$Position[keeps], g1a, pch = 19, cex=0.5, col = 'red')
    lines(((1:length(keeps2) - 0.5)*10^6)[keeps3], (AlleleFrac[names(keeps3),cell]*40)-20, col='purple2', lwd = 2, lty='dashed')
    lines(genes[keeps4,2], (Hap_all3[keeps4,cell]*40)-20, col='black', lwd=2)

    axis(1, at=seq(0,max(genes$Position[keeps]),(5*10^6)), labels=seq(0,max(genes$Position[keeps])/(10^6),5))
    axis(2, at=seq(-20,0,10), labels=c(20,10,0), col = 'blue', col.axis = 'blue', las=2)
    axis(2, at=seq(20,0,-10), labels=c(20,10,0), col = 'red', col.axis = 'red', las=2)
    axis(2, at=0, labels=0, las=2)
    axis(4, at=seq(-20,20,10), labels=seq(0,100,25), col = 'purple2', col.axis = 'purple2', las=2)
    mtext("% of transcripts from Col-0 allele", side=4, line=3, col="purple2")
    mtext("transcript counts", side=2, line=4)
    mtext("Col-0", side=2, at=10, line=2.5, col="red")
    mtext("Ler-0", side=2, at=-10, line=2.5, col="blue")
}
plotCellChrHapVsBin(cell = 'A194-205_55s', chr=4)

svg("CheckUM_Mono.svg", width=30, height=20)
par(mfrow = c(4, 6))
for (cell in UM_Mono){
    plotCellChrHapVsBin(cell, chr=1)
}
dev.off()

par(mfrow = c(1,1))


geneUse = which(abs(rowMeans(Hap_all2[,over30k], na.rm=T) - .5) < .4)  # Exclude genes with >90% of all transcripts mapping to the same allele across all samples
Hap_all2[-geneUse,] = NA
geneUse2 = which(rowSums(!is.na(Hap_all2[,UM_Mono])) >= 10) # Exclude genes with <5 Monoallelic UMs with data
Hap_all2[-geneUse2,] = NA
plot(rowMeans(Hap_all2[,UM_Mono], na.rm=T))

plot(rowMeans(Hap_all2[geneUse2[grep("AT1", names(geneUse2))],UM_Mono],na.rm=T))

hist((rowMeans(Hap_all2[,over30k], na.rm=T)))

hist((rowMeans(Hap_all2[,UM_Mono], na.rm=T)))
hist((rowMeans(Hap_all2[,UM_Bi], na.rm=T)))

hist((rowMeans(Hap_all2[,BM_Mono], na.rm=T)))
hist((rowMeans(Hap_all2[,BM_Bi], na.rm=T)))

hist((rowMeans(Hap_all2[,TP_Mono], na.rm=T)))
hist((rowMeans(Hap_all2[,TP_Bi], na.rm=T)))

hist((rowMeans(Hap_all2[,Mono], na.rm=T)), xlim=c(0,1))
summary(rowMeans(Hap_all2[,Mono], na.rm=T))
hist((rowMeans(Hap_all2[,Bi], na.rm=T)), xlim=c(0,1))
summary(rowMeans(Hap_all2[,Bi], na.rm=T))

HapGroup = c(rep("over20k", times=length(rowMeans(Hap_all2[,over30k], na.rm=T))), 
    rep("Mono", times=length(rowMeans(Hap_all2[,Mono], na.rm=T))), rep("Bi", times=length(rowMeans(Hap_all2[,Bi], na.rm=T))), 
    rep("UM Mono", times=length(rowMeans(Hap_all2[,UM_Mono], na.rm=T))), rep("UM Bi", times=length(rowMeans(Hap_all2[,UM_Bi], na.rm=T))), 
    rep("BM Mono", times=length(rowMeans(Hap_all2[,BM_Mono], na.rm=T))), rep("BM Bi", times=length(rowMeans(Hap_all2[,BM_Bi], na.rm=T))),
    rep("TP Mono", times=length(rowMeans(Hap_all2[,TP_Mono], na.rm=T))), rep("TP Bi", times=length(rowMeans(Hap_all2[,TP_Bi], na.rm=T))))
HapGroup = factor(HapGroup, levels=c('over20k','Mono', 'Bi','UM Mono', 'UM Bi', 'BM Mono', 'BM Bi', 'TP Mono', 'TP Bi'))


HapByGene = data.frame(Hap = c(rowMeans(Hap_all2[,over30k], na.rm=T), rowMeans(Hap_all2[,Mono], na.rm=T), 
    rowMeans(Hap_all2[,Bi], na.rm=T), rowMeans(Hap_all2[,UM_Mono], na.rm=T), rowMeans(Hap_all2[,UM_Bi], na.rm=T), 
    rowMeans(Hap_all2[,BM_Mono], na.rm=T), rowMeans(Hap_all2[,BM_Bi], na.rm=T), 
    rowMeans(Hap_all2[,TP_Mono], na.rm=T), rowMeans(Hap_all2[,TP_Bi], na.rm=T)), 
    Group = HapGroup)


HapByGene = data.frame(Hap = c(rowMeans(Hap_all2[,over30k], na.rm=T), rowMeans(Hap_all2[,Mono], na.rm=T), rowMeans(Hap_all2[,Bi], na.rm=T), rowMeans(Hap_all2[,UM_Mono], na.rm=T), rowMeans(Hap_all2[,UM_Bi], na.rm=T), rowMeans(Hap_all2[,BM_Mono], na.rm=T), rowMeans(Hap_all2[,BM_Bi], na.rm=T)), 
    Group = c(rep("over20k", times=length(rowMeans(Hap_all2[,over30k], na.rm=T))), rep("Mono", times=length(rowMeans(Hap_all2[,Mono], na.rm=T))), rep("Bi", times=length(rowMeans(Hap_all2[,Bi], na.rm=T))), rep("UM Mono", times=length(rowMeans(Hap_all2[,UM_Mono], na.rm=T))), rep("UM Bi", times=length(rowMeans(Hap_all2[,UM_Bi], na.rm=T))), rep("BM Mono", times=length(rowMeans(Hap_all2[,BM_Mono], na.rm=T))), rep("BM Bi", times=length(rowMeans(Hap_all2[,BM_Bi], na.rm=T)))))

boxplot(Hap ~ Group, data=HapByGene)
pairwise.wilcox.test(HapByGene$Hap, HapByGene$Group, p.adjust = 'holm')  # p-values







Heatmap(Hap_all[,over30k][ords])

Heatmap(((Hap_all[,over30k][,ords]+1)/2), name = 'Haplotype', heatmap_legend_param = list(at = c(0, 1), labels = c(0, 100)),
    top_annotation = HeatmapAnnotation(Stage = stages[ords],
    "% Matching Haplotype" = HapFrac[ords],
    col = list(Stage = c("Tetrad" = "#EB1E2C", "UM Bi" = "#F9A729", "UM Mono" = "#F9A729",
    "UM/BM Bi" = "#F9D23C", "UM/BM Mono" = "#F9D23C", "BM Bi" = "#5FBB68", "BM Mono" = "#5FBB68", 
    "BM/TP" = "#64CDCC", "TP Bi" = "#A4A4D5", "TP Mono" = "#A4A4D5"), 
    "% Matching Haplotype" = FracMono_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), 
    cluster_rows=F, cluster_columns=F, show_row_names = FALSE, show_column_names = FALSE)



table(stages)
#     BM  BM/Tri tetrads     Tri      UM   UM/BM
#     77       3       4      87      49      14 
table(At_meta[,9])
#     BM  BM/Tri tetrads     Tri      UM   UM/BM
#    304      48       8     304     360      96

#Assuming (4/8)/4 may be the success ratio for separated tetrads, I would need to collect 20 strips for 20 separated tetrads
#If we want these, collect 2 plates w/ 20 strips of sep tetrads, 4 strips of pollen (possibly shedding?)
#1 plate should yield ~10 separated tetrads by this estimate