setwd('C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R')

load('4_2024_At_Spike_ins5.RData')
colnames(D) = sub('At1', 'A111-113', sub('At2', 'A114-116', sub('At3', 'A117-119', sub('At4', 'A120-122', colnames(D)))))
colnames(g1) = sub('At1', 'A111-113', sub('At2', 'A114-116', sub('At3', 'A117-119', sub('At4', 'A120-122', colnames(g1)))))
colnames(g2) = sub('At1', 'A111-113', sub('At2', 'A114-116', sub('At3', 'A117-119', sub('At4', 'A120-122', colnames(g2)))))
D_4 = D
g1_4 = g1
g2_4 = g2

load('7-8_2024_At_4.RData')
D_7_8 = D
g1_7_8 = g1
g2_7_8 = g2


load('11_2024_At.RData')
colnames(D) = sub('S2_L002_', '', sub('S1-8_', '', sub('S3_L002_', '', colnames(D))))
colnames(g1) = sub('S2_L002_', '', sub('S1-8_', '', sub('S3_L002_', '', colnames(g1))))
colnames(g2) = sub('S2_L002_', '', sub('S1-8_', '', sub('S3_L002_', '', colnames(g2))))
D_11 = D
g1_11 = g1
g2_11 = g2

load('12_2024_At.RData')
colnames(D) = sub('S187_L007_', '', sub('S188_L007_', '', colnames(D)))
colnames(g1) = sub('S187_L007_', '', sub('S188_L007_', '', colnames(g1)))
colnames(g2) = sub('S187_L007_', '', sub('S188_L007_', '', colnames(g2)))
D_12 = D
g1_12 = g1
g2_12 = g2


D = cbind(D_4, D_7_8, D_11, D_12)
g1 = cbind(g1_4, g1_7_8, g1_11, g1_12)
g2 = cbind(g2_4, g2_7_8, g2_11, g2_12)


library(readxl)
At_Stages <- read_excel("C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R/At_Stages2.xlsx")

At_meta <- read_excel("C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R/At_meta2.xlsx")
rownames(At_meta) = At_meta$Sample

D = D[rownames(genes),]
g1 = g1[rownames(genes),]
g2 = g2[rownames(genes),]

D = D[,At_meta$Sample]
g1 = g1[,At_meta$Sample]
g2 = g2[,At_meta$Sample]

No_cell = At_meta$Sample[which(At_meta$No_cell_well == 'Y')]



plotChr = function (cell, chr = 1, pad = 3) 
{
    Cdat = data.frame(f_col0 = AlleleFrac[, cell], Chr = floor(as.numeric(rownames(AlleleFrac))/10^6), 
        Position = (as.numeric(rownames(AlleleFrac))%%10^6) + 
            0.5)
    Cdat = Cdat[Cdat$Chr == chr, ]
    ggplot(Cdat) + geom_rect(data = data.frame(xmin = -pad, xmax = max(genes[,2])/10^6 + 
        pad, ymin = 0, ymax = 1), aes(xmin = xmin, xmax = xmax, 
        ymin = ymin, ymax = ymax), fill = "#EEEEEE") + geom_point(aes(y = f_col0, 
        x = Position), cex = 3) + geom_hline(yintercept = 0.5, 
        linetype = "dashed") + theme(panel.background = element_blank(), 
        axis.title = element_blank(), panel.border = element_blank(), 
        panel.grid = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) + scale_y_continuous(breaks = seq(0, 
        1, 0.25), labels = c("0%", "", "50%", "", "100%"), limits = c(-0.4, 
        1.05)) + scale_x_continuous(expand = c(0, 0)) + annotate("segment", 
        x = -pad, xend = -pad, y = 0, yend = 1) + theme(plot.margin = margin(0, 
        0, 0.15, 0, "cm"))
}

plotScaleBar = ggplot() + scale_x_continuous(expand=c(0,0), limits = c(-3, max(genes[,2])/10^6 + 3), breaks = seq(0,300,5)) + theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.x=element_line(), plot.margin = margin(0,0,0,0,'cm')) + xlab('Chromosome position (Mb)')

plotCell2 = function (cell) 
{
    annotate_figure(ggarrange(plotChr(cell, chr = 1), plotChr(cell, 
        chr = 2), plotChr(cell, chr = 3), plotChr(cell, chr = 4), 
        plotChr(cell, chr = 5), plotScaleBar, ncol = 1, nrow = 6, 
        align = "v", heights = c(rep(1, 5), 0.4)), left = text_grob("          % Transcripts from Col-0 allele", 
        rot = 90, size = 10), top = cell)
}

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
AlleleFrac2 = AlleleFrac[,which(colSums(D) >= 10000)] ###491/1120 passed >= 10000, 5 of these are no cell controls



HapCallV4 = function(cell, chr, recs=1.5, mx = (2.5), ss=1, niter = 10000, burnin = 2000) {
	keeps0 = names(which(abs(rowMeans(g1/(g1+g2), na.rm=T) - .5) < .3)) #exclude genes with >80% transcripts matching the same allele
    keeps = names(which(((g1[,cell] + g2[,cell]) > 0) & (genes$Chr == chr)))
    
    g1b = g1[keeps[keeps %in% keeps0],colnames(AlleleFrac2)]
    g2b = g2[keeps[keeps %in% keeps0],colnames(AlleleFrac2)]
    g1bFracMean = rowMeans(g1b/(g1b+g2b), na.rm=T)
    g2bFracMean = rowMeans(g2b/(g1b+g2b), na.rm=T)

	gA = g1b[,cell]
	gB = g2b[,cell]
	ScoreAB = pbinom(gB,gA+gB,g2bFracMean,log.p=T) - pbinom(gA,gA+gB,g1bFracMean,log.p=T)
	k = -log(recs/length(ScoreAB))
	ScoreAB[abs(ScoreAB) > k/mx] = sign(ScoreAB[abs(ScoreAB) > k/mx])*k/mx

	Hap = rep(-sign(sum(ScoreAB)), length(ScoreAB))
	Scores = c(-sum(ScoreAB*Hap) - sum(abs(diff(Hap)/2))*k, rep(NA,niter))

	alpha = 2*ScoreAB*Hap
	beta0 = -c(k, rep(2*k,length(Hap)-2), k)

	l = length(Hap)
	XOs = abs(diff(Hap))
	deltaK = (XOs-1)*k*2
	beta = beta0 + (c(0, XOs) + c(XOs,0))*k
    	HapOut = Hap*0
    
	for (i in 1:niter) {
        ar = alpha + beta - log(runif(l))*ss
		swaps = ar >= 0
		if (any(swaps)) {
			for (j in 2:l) {
				swaps[j] = ar[j] >= swaps[j-1]*deltaK[j-1]
			}
			alpha[swaps] = -alpha[swaps]
			Hap[swaps] = -Hap[swaps]
			XOs = abs(diff(Hap))
			deltaK = (XOs-1)*k*2
			beta = beta0 + (c(0, XOs) + c(XOs,0))*k
			Scores[i] = -sum(ScoreAB*Hap) - sum(XOs)*k/2
			if (i > burnin) { HapOut = HapOut + Hap }
		}
		
	}
	HapOut = HapOut/sum(!is.na(Scores[-c(1:burnin)]))
	Position = genes[names(ScoreAB),2]
    Genes = names(ScoreAB)
    return(list(HapOut=HapOut, Scores=Scores, Pos=Position, Genes=Genes))
}


set.seed(1)

cell = "A230-241_79s"


AlleleFrac_genes = (g1/(g1+g2))

HapCallCell = function(cell){
    HapChr1 = HapCallV4(cell,chr=1,recs=1.5, mx=(2.5), ss=1, niter=10000, burnin=4000)
    HapChr2 = HapCallV4(cell,chr=2,recs=1.5, mx=(2.5), ss=1, niter=10000, burnin=4000)
    HapChr3 = HapCallV4(cell,chr=3,recs=1.5, mx=(2.5), ss=1, niter=10000, burnin=4000)
    HapChr4 = HapCallV4(cell,chr=4,recs=1.5, mx=(2.5), ss=1, niter=10000, burnin=4000)
    HapChr5 = HapCallV4(cell,chr=5,recs=1.5, mx=(2.5), ss=1, niter=10000, burnin=4000)
    HapCell = data.frame(Hap = c(HapChr1$HapOut, HapChr2$HapOut, HapChr3$HapOut, HapChr4$HapOut, HapChr5$HapOut), 
                         Pos = c(HapChr1$Pos, HapChr2$Pos, HapChr3$Pos, HapChr4$Pos, HapChr5$Pos),
                         Genes = c(HapChr1$Genes, HapChr2$Genes, HapChr3$Genes, HapChr4$Genes, HapChr5$Genes))
    rownames(HapCell) = HapCell$Genes
    HapCallGenes = rownames(D[which(rownames(D) %in% rownames(HapCell)),])
    Hap = HapCell$Hap
    names(Hap) = HapCell$Genes
    HapCell$Hap = round(HapCell$Hap)
    HapCell$Hap[which(round(HapCell$Hap) == 0)] = NA
    HapCell$Hap[which(round(HapCell$Hap) == -1)] = 0
    AlleleFrac_genes_cell = round(AlleleFrac_genes[HapCallGenes,cell])
    HapMatch = AlleleFrac_genes_cell == HapCell$Hap #for each gene, does AlleleFrac match expected haplotype?
    HapFrac = length(which(AlleleFrac_genes_cell == HapCell$Hap))/length(HapCallGenes)
    return(list(HapMatch = HapMatch, HapFrac = HapFrac, Hap = Hap))
}


cells = colnames(AlleleFrac2)

HapMatch_all = matrix(data = NA, nrow = nrow(D), ncol = ncol(AlleleFrac2))
rownames(HapMatch_all) = rownames(D)
colnames(HapMatch_all) = colnames(AlleleFrac2)

HapFrac_all = matrix(data = NA, nrow = 1, ncol = ncol(AlleleFrac2))
colnames(HapFrac_all) = colnames(AlleleFrac2)

a1 = proc.time()
for (cell in cells){
    HapMatch_cell = HapCallCell(cell)
    HapMatch_all[match(names(HapMatch_cell$HapMatch), rownames(HapMatch_all)),cell] = HapMatch_cell$HapMatch
    HapFrac_all[cell] = HapMatch_cell$HapFrac
}
a2 = proc.time()
a2-a1
#    user   system  elapsed
#11765.23   153.72 11969.48  #It took ~200 mins (3 hr 20 min) for this to run
###Something odd happened with HapFrac_all adding on new values to the end instead of replacing existing NA values
#This can be fixed by subsetting to keep the second half of the vector and matching the order to AlleleFrac2

HapFrac_all_fixed = HapFrac_all[495:988]
HapFrac_all_fixed = HapFrac_all_fixed[colnames(AlleleFrac2)]
#saving this to avoid having to rerun the 3 hr 20 min code to get HapMatch_all and HapFrac_all
save(D, g1, g2, genes, HapMatch_all, HapFrac_all_fixed, No_cell, plotChr, plotScaleBar, plotCell2, BIN2, g1_bin, 
    g2_bin, AlleleFrac, AlleleFrac2, HapCallV4, AlleleFrac_genes, HapCallCell, file = "At_data_9_2025.RData")


#I tested and found the problem with HapFrac_all can be fixed by creating a vector instead of a matrix. 
#I will use the below code in the future:

cells = colnames(AlleleFrac2)

HapMatch_all = matrix(data = NA, nrow = nrow(D), ncol = ncol(AlleleFrac2))
rownames(HapMatch_all) = rownames(D)
colnames(HapMatch_all) = colnames(AlleleFrac2)

HapFrac_all = vector(length = ncol(AlleleFrac2))
names(HapFrac_all) = colnames(AlleleFrac2)


a1 = proc.time()
for (cell in cells){
    HapMatch_cell = HapCallCell(cell)
    HapMatch_all[match(names(HapMatch_cell$HapMatch), rownames(HapMatch_all)),cell] = HapMatch_cell$HapMatch
    HapFrac_all[cell] = HapMatch_cell$HapFrac
}
a2 = proc.time()
a2-a1




HapMatch_all = matrix(data = NA, nrow = nrow(D), ncol = ncol(AlleleFrac2))
rownames(HapMatch_all) = rownames(D)
colnames(HapMatch_all) = colnames(AlleleFrac2)

HapFrac_all = vector(length = ncol(AlleleFrac2))
names(HapFrac_all) = colnames(AlleleFrac2)

Hap_all = HapMatch_all

a1 = proc.time()
for (cell in cells[c(1,10,20,30,100,400)]){
    HapMatch_cell = HapCallCell(cell)
    HapMatch_all[match(names(HapMatch_cell$HapMatch), rownames(HapMatch_all)),cell] = HapMatch_cell$HapMatch
    HapFrac_all[cell] = HapMatch_cell$HapFrac
    Hap_all[match(names(HapMatch_cell$HapMatch), rownames(HapMatch_all)),cell] = HapMatch_cell$Hap
}
a2 = proc.time()
a2-a1


###
load("At_data_9_2025.RData")

library('ComplexHeatmap')

FracMono = 100*colMeans(abs(AlleleFrac2 - .5) >= .3, na.rm=T)
plot(FracMono, HapFrac_all_fixed, cex = 2, pch=20)

library(circlize)
FracMono_col = colorRamp2(c(0, 100), c("white", "purple4"))
UMI_col = colorRamp2(c(4, 5.3), c("white", "forestgreen"))
HapFrac_col = colorRamp2(c(0.47, 0.97), c("white", "purple4"))

Heatmap(AlleleFrac2, cluster_rows=F, cluster_columns=T)                   
Heatmap(AlleleFrac2, name = 'At AlleleFrac',
    top_annotation = HeatmapAnnotation("UMIcounts" = log(colSums(D[,colnames(AlleleFrac2)]),10),
    "No cell" = At_meta$No_cell_well[which(colnames(D) %in% colnames(AlleleFrac2))], "FracMono" = FracMono,
    "Stage" = At_meta$Stage[which(colnames(D) %in% colnames(AlleleFrac2))],
    "HapFrac" = HapFrac_all_fixed,
    col = list(FracMono = FracMono_col, UMIcounts = UMI_col, Stage = c("tetrad" = "#eeeeee", "UM" = "#cccccc", "UM/BM" = "#aaaaaa", "BM" = "#777777", "BM/Tri" = "#444444", "Tri" = "#111111"), 
    "No cell" = c("Y" = "#111111", "N" = "#eeeeee"), "HapFrac" = HapFrac_col)),
    col = colorRampPalette(c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), 
    cluster_rows=F, cluster_columns=F, show_row_names = F, show_column_names = F)



#Used to check how HapCallCell$Hap looked
plot(HapCallCell$Hap, cex=1.5, pch=19, ylim = c(-1,1))
    abline(v=max(grep("AT1", HapCallCell$Genes)))
    abline(v=max(grep("AT2", HapCallCell$Genes)))
    abline(v=max(grep("AT3", HapCallCell$Genes)))
    abline(v=max(grep("AT4", HapCallCell$Genes)))
	abline(h=1-0.025)
	abline(h=0.025-1)


plot(HapCallChr5$Scores, cex=1.5, pch=19, type='l')
abline(v=4000, col="red2")


c(1,10,20,30,100,400)
plot(genes[1:max(grep("AT1", rownames(D))),2], Hap_all[1:max(grep("AT1", rownames(D))),1], cex=1.5, pch=19, ylim = c(-1,1), col = "red2")
    abline(v=max(grep("AT1", rownames(D))))
    abline(v=max(grep("AT2", rownames(D))))
    abline(v=max(grep("AT3", rownames(D))))
    abline(v=max(grep("AT4", rownames(D))))
	abline(h=1-0.025)
	abline(h=0.025-1)
    points(genes[names(which((g1[,cells[1]] + g2[,cells[1]]) > 0)),2], 2*(AlleleFrac_genes[names(which((g1[,cells[1]] + g2[,cells[1]]) > 0)),1])-1, cex=1.5, pch=19)

plotChr_Hap = function (cell, chr = 1, pad = 3, chro = "AT1") 
{
    Cdat = data.frame(f_col0 = AlleleFrac[, cell], Chr = floor(as.numeric(rownames(AlleleFrac))/10^6), 
        Position = (as.numeric(rownames(AlleleFrac))%%10^6)*10^6)
    Cdat = Cdat[Cdat$Chr == chr, ]
    HapDat = data.frame(Hap = (Hap_all[min(grep(chro, rownames(D))):max(grep(chro, rownames(D))),cell]+1)/2, 
    loc = genes[min(grep(chro, rownames(D))):max(grep(chro, rownames(D))),2])
    ggplot(Cdat) + geom_rect(data = data.frame(xmin = -pad, xmax = max(genes[,2])/10^6 + 
        pad, ymin = 0, ymax = 1), aes(xmin = xmin, xmax = xmax, 
        ymin = ymin, ymax = ymax), fill = "#EEEEEE") + geom_point(data = Cdat, aes(y = f_col0, 
        x = Position), cex = 3) + geom_hline(yintercept = 0.5, 
        linetype = "dashed") + theme(panel.background = element_blank(), 
        axis.title = element_blank(), panel.border = element_blank(), 
        panel.grid = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) + scale_y_continuous(breaks = seq(0, 
        1, 0.25), labels = c("0%", "", "50%", "", "100%"), limits = c(-0.4, 
        1.05)) + scale_x_continuous(expand = c(0, 0)) + annotate("segment", 
        x = -pad, xend = -pad, y = 0, yend = 1) + theme(plot.margin = margin(0, 
        0, 0.15, 0, "cm")) + geom_point(data = HapDat, aes(y = Hap, x = loc), cex = 3, col = "red2")
}

plotScaleBar = ggplot() + scale_x_continuous(expand=c(0,0), limits = c(-3, max(genes[,2])/10^6 + 3), breaks = seq(0,300,5)) + theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.x=element_line(), plot.margin = margin(0,0,0,0,'cm')) + xlab('Chromosome position (Mb)')

plotCell_Hap = function (cell) 
{
    annotate_figure(ggarrange(plotChr(cell, chr = 1, chro = "AT1"), 
        plotChr(cell, chr = 2, chro = "AT2"), plotChr(cell, chr = 3, chro = "AT3"), 
        plotChr(cell, chr = 4, chro = "AT4"), plotChr(cell, chr = 5, chro = "AT5"), 
        plotScaleBar, ncol = 1, nrow = 6, align = "v", heights = c(rep(1, 5), 0.4)), 
        left = text_grob("          % Transcripts from Col-0 allele", 
        rot = 90, size = 10), top = cell)
}

