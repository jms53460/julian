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
    HapCell$Hap = round(HapCell$Hap)
    HapCell$Hap[which(round(HapCell$Hap) == 0)] = NA
    HapCell$Hap[which(round(HapCell$Hap) == -1)] = 0
    AlleleFrac_genes_cell = round(AlleleFrac_genes[HapCallGenes,cell])
    HapMatch = AlleleFrac_genes_cell == HapCell$Hap
    HapFrac = length(which(AlleleFrac_genes_cell == HapCell$Hap))/length(HapCallGenes)
    return(list(HapMatch = HapMatch, HapFrac = HapFrac))
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

a1 = proc.time()
HapMatch_cell = HapCallCell(cell)
    HapMatch_all[match(names(HapMatch_cell$HapMatch), rownames(HapMatch_all)),cell] = HapMatch_cell$HapMatch
    HapFrac_all[cell] = HapMatch_cell$HapFrac
a2 = proc.time()
a2-a1


plot(HapCallCell$Hap, cex=1.5, pch=19, ylim = c(-1,1))
    abline(v=max(grep("AT1", HapCallCell$Genes)))
    abline(v=max(grep("AT2", HapCallCell$Genes)))
    abline(v=max(grep("AT3", HapCallCell$Genes)))
    abline(v=max(grep("AT4", HapCallCell$Genes)))
	abline(h=1-0.025)
	abline(h=0.025-1)

plot(HapCallChr5$Scores, cex=1.5, pch=19, type='l')
abline(v=4000, col="red2")






plotCellChr = function(cell, chr, maxs = 20) {
	keeps = ((g1[,cell] + g2[,cell]) > 0) & (genes$Chr == chr)
    keeps2 = which(round(as.numeric(rownames(AlleleFrac))/(10^6)) == chr)
    keeps3 = which(!is.na(AlleleFrac[keeps2,cell]))

    g1a = g1[keeps,cell]
	g1a[g1a > maxs] = maxs
    g1b = g1[,cell]
    g1b[g1b > maxs] = maxs
	g2a = g2[keeps,cell]
	g2a[g2a > maxs] = maxs
    g2b = g2[,cell]
    g2b[g2b > maxs] = maxs

    par(mar = c(5.1, 5.1, 4.1, 4.1))
	plot(genes$Position[keeps], g1a, pch = 19, cex=0.5, col = 'red', ylab=NA, xlab = 'chromosome position (Mb)', ylim = c(-maxs,maxs), xaxt = 'n', yaxt = 'n')
	abline(h=0)
	for (ii in which(keeps)) { lines(rep(genes$Position[ii], times=2), c(g1b[ii], -g2b[ii]), col='#777777') }
	points(genes$Position[keeps], -g2a, pch = 19, cex=1.5, col = 'blue')
	points(genes$Position[keeps], g1a, pch = 19, cex=1.5, col = 'red')
    lines(((1:length(keeps2) - 0.5)*10^6)[keeps3], (AlleleFrac[names(keeps3),cell]*40)-20, col='black', lwd = 6, lty='dashed')

    points(HapCallEx$Pos, HapCallEx$HapOut*20, cex=1.5, pch=19, ylim = c(-1,1), col="purple2")

    axis(1, at=seq(0,max(genes$Position[keeps]),(5*10^6)), labels=seq(0,max(genes$Position[keeps])/(10^6),5))
    axis(2, at=seq(-20,0,10), labels=c(20,10,0), col = 'blue', col.axis = 'blue', las=2)
    axis(2, at=seq(20,0,-10), labels=c(20,10,0), col = 'red', col.axis = 'red', las=2)
    axis(2, at=0, labels=0, las=2)
    axis(4, at=seq(-20,20,10), labels=seq(0,100,25), col = 'black', col.axis = 'black', las=2)
    mtext("% of transcripts from Col-0 allele", side=4, line=3, col="black")
    mtext("transcript counts", side=2, line=4)
    mtext("Col-0", side=2, at=10, line=2.5, col="red")
    mtext("Ler-0", side=2, at=-10, line=2.5, col="blue")
}

plotCellChr(cell, chr)
