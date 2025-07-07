avgs = NULL

for (j in 1:100) {


##################################################

niter = 10^4
burnin = 100
bias = 0.1
ginf = 10
gA = c(rbinom(30, ginf, .5 + bias), rbinom(15, ginf, .5 - bias), rbinom(15, ginf, .5 + bias))
gB = ginf-gA

ScoreAB = pbinom(gB,gA+gB,0.5,log.p=T) - pbinom(gA,gA+gB,0.5,log.p=T)
k = -log(1.5/length(ScoreAB))
ScoreAB[abs(ScoreAB) > k/1.25] = sign(ScoreAB[abs(ScoreAB) > k/1.25])*k/1.25
blockScore = ScoreAB

	Hap = rep(-sign(sum(blockScore)), length(blockScore))
	Scores = c(-sum(blockScore*Hap) - sum(abs(diff(Hap)/2))*k, rep(NA,niter))

	alpha = 2*blockScore*Hap
	beta0 = -c(k, rep(2*k,length(Hap)-2), k)

	l = length(Hap)
	XOs = abs(diff(Hap))
	deltaK = (XOs-1)*k*2
	beta = beta0 + (c(0, XOs) + c(XOs,0))*k
    	HapOut = Hap*0
    
	for (i in 1:niter) {
		ar = alpha + beta - log(runif(l))
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
			Scores[i] = -sum(blockScore*Hap) - sum(XOs)*k/2
			if (i > burnin) { HapOut = HapOut + Hap }
		}
		
	}
	HapOut = (HapOut/sum(!is.na(Scores[-c(1:burnin)])) + 1)/2
	


plot(HapOut, type = 'l', lwd=2, yaxs = 'i', ylim = c(0,1))
points(gA/(gA+gB), pch = 19, col='red')
abline(v=30.5,lty=2, lwd=2)
abline(v=45.5,lty=2, lwd=2)

##################################################

avgs = cbind(avgs, HapOut)
}

plot(rowMeans(avgs), type='l', lwd = 2, yaxs='i', ylim = c(0,1))
abline(v=30.5,lty=2, lwd=2)
abline(v=45.5,lty=2, lwd=2)
abline(h=.5)

avgs3 = avgs
avgs3[abs(avgs3-.5) <= (.5-.05/2)] = (avgs3[abs(avgs3-.5) <= (.5-.05/2)]-.5)/2+.5
Heatmap(t(avgs3), cluster_columns=F)

avgs2 = avgs
avgs2[31:45,] = 1 - avgs[31:45,]
hist(avgs2, breaks = 20)
abline(h = sum(avgs2 <= .5)/10, lty = 2, lwd = 2, col = 'red')


HapCallV3 = function(cell, chr, niter = 10000, burnin = 2000) {
	keeps = ((g1[,cell] + g2[,cell]) > 0) & (genes$Chr == chr)
	gA = g1[keeps,cell]
	gB = g2[keeps,cell]
	ScoreAB = pbinom(gB,gA+gB,0.5,log.p=T) - pbinom(gA,gA+gB,0.5,log.p=T)
	k = -log(1.5/length(ScoreAB))
	ScoreAB[abs(ScoreAB) > k/1.25/2] = sign(ScoreAB[abs(ScoreAB) > k/1.25/2])*k/1.25/2

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
		ar = alpha + beta - log(runif(l))
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
    return(list(HapOut=HapOut, Scores=Scores, Pos=Position))
}



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
    return(list(HapOut=HapOut, Scores=Scores, Pos=Position))
}

#HapCallV3b = function(cell, chr, recs=1.5, mx=2.5, niter = 10000, burnin = 2000) {
#	keeps0 = names(which(abs(rowMeans(g1/(g1+g2), na.rm=T) - .5) < .4))
#   keeps = names(which(((g1[,cell] + g2[,cell]) > 0) & (genes$Chr == chr)))
    
#	gA = g1[keeps[keeps %in% keeps0],cell]
#	gB = g2[keeps[keeps %in% keeps0],cell]
#	ScoreAB = pbinom(gB,gA+gB,0.5,log.p=T) - pbinom(gA,gA+gB,0.5,log.p=T)
#	k = -log(recs/length(ScoreAB))
#	ScoreAB[abs(ScoreAB) > k/mx] = sign(ScoreAB[abs(ScoreAB) > k/mx])*k/mx

#	Hap = rep(-sign(sum(ScoreAB)), length(ScoreAB))
#	Scores = c(-sum(ScoreAB*Hap) - sum(abs(diff(Hap)/2))*k, rep(NA,niter))

#	alpha = 2*ScoreAB*Hap
#	beta0 = -c(k, rep(2*k,length(Hap)-2), k)

#	l = length(Hap)
#	XOs = abs(diff(Hap))
#	deltaK = (XOs-1)*k*2
#	beta = beta0 + (c(0, XOs) + c(XOs,0))*k
 #   	HapOut = Hap*0
    
#	for (i in 1:niter) {
#		ar = alpha + beta - log(runif(l))
#		swaps = ar >= 0
#		if (any(swaps)) {
#			for (j in 2:l) {
#				swaps[j] = ar[j] >= swaps[j-1]*deltaK[j-1]
#			}
#			alpha[swaps] = -alpha[swaps]
#			Hap[swaps] = -Hap[swaps]
#			XOs = abs(diff(Hap))
#			deltaK = (XOs-1)*k*2
#			beta = beta0 + (c(0, XOs) + c(XOs,0))*k
#			Scores[i] = -sum(ScoreAB*Hap) - sum(XOs)*k/2
#		}
#		if (i > burnin) { HapOut = HapOut + Hap }		
#	}
#	HapOut = HapOut/(niter - burnin)
#	Position = genes[names(ScoreAB),2]

 #   return(list(HapOut=HapOut, Scores=Scores, Pos=Position))
#}

plotCell2("A230-241_79s")
cell = "A230-241_79s"
chr = 5

plotCell2("A114-116_42s")
cell = "A114-116_42s"
chr = 5

cell = "A288-301_75s"

#for genes
set.seed(1)
burnin=2000
#a1 = proc.time()
HapCallEx = HapCallV3(cell,chr, niter=10000, burnin=2000)
#a2 = proc.time()
#a2-a1
plot(HapCallEx$Pos, HapCallEx$HapOut, cex=1.5, pch=19, ylim = c(-1,1))
	keeps2 = which(round(as.numeric(rownames(AlleleFrac))/(10^6)) == chr)
	keeps3 = which(!is.na(AlleleFrac[keeps2,cell]))
	lines(((1:length(keeps2) - 0.5)*10^6)[keeps3], (AlleleFrac[names(keeps3),cell]*2)-1, col='purple2', lwd = 2, lty='dashed')
	abline(h=1-0.025)
	abline(h=0.025-1)




set.seed(1)
burnin=4000
#a1 = proc.time()
HapCallEx = HapCallV4(cell,chr,recs=1.5, mx=(2.5), ss=1, niter=10000, burnin=4000)
#a2 = proc.time()
#a2-a1
plot(HapCallEx$Pos, HapCallEx$HapOut, cex=1.5, pch=19, ylim = c(-1,1))
	keeps2 = which(round(as.numeric(rownames(AlleleFrac))/(10^6)) == chr)
	keeps3 = which(!is.na(AlleleFrac[keeps2,cell]))
	lines(((1:length(keeps2) - 0.5)*10^6)[keeps3], (AlleleFrac[names(keeps3),cell]*2)-1, col='purple2', lwd = 2, lty='dashed')
	abline(h=1-0.025)
	abline(h=0.025-1)



plot(HapCallEx$Scores, cex=1.5, pch=19, type='l')
abline(v=burnin, col="red2")






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


svg('At_UMI_counts_plot2.svg', width=5, height=4)
plotCellChr(cell = 'A230-241_79s', 5)
dev.off()


svg('At_UMI_counts_plot3.svg', width=5, height=4)
plotCellChr(cell = 'A288-301_75s', 5)
dev.off()



cell1 = "A230-241_79s"
c1_chr = 5
cell2 = "A114-116_42s"
c2_chr = 5
cellS = "A288-301_75s"
cS_chr = 5

rowMeans(g1[,c(14,16,30,31)])


###simulate some cells (cellS)


metaFn = function(cell1, c1_chr, cell2, c2_chr, cellS, cS_chr, recs=1.5, mx = (2.5), ss=1, niter = 10000, burnin = 2000) {
par(mfrow = c(3,2))
HapCallEx = HapCallV4(cell1, c1_chr, mx = mx, recs = recs)
plotCellChr(cell=cell1, chr=c1_chr)
plot(HapCallEx$Scores, cex=1.5, pch=19, type='l')
abline(v=burnin, col="red2")

HapCallEx = HapCallV4(cell2, c2_chr, mx = mx, recs = recs)
plotCellChr(cell=cell2, chr=c2_chr)
plot(HapCallEx$Scores, cex=1.5, pch=19, type='l')
abline(v=burnin, col="red2")

HapCallEx = HapCallV4(cellS, cS_chr, mx = mx, recs = recs)
plotCellChr(cell=cellS, chr=cS_chr)
plot(HapCallEx$Scores, cex=1.5, pch=19, type='l')
abline(v=burnin, col="red2")
}

metaFn(cell1, c1_chr, cell2, c2_chr, cellS, cS_chr)

metaFn = function(mx, recs, etc) {
par(mfrow = c(3,2))
result = HapCallV4(cell1, mx = mx, recs = recs, etc)
plotA(result)
plotB(result)

result = HapCallV4(cell2, mx = mx, recs = recs, etc)
plotA(result)
plotB(result)

result = HapCallV4(cellS, mx = mx, recs = recs, etc)
plotA(result)
plotB(result)
}
