avgs = NULL

for (j in 1:100) {



##################################################
bias = 0.2
g1s = c(rbinom(30, 10, .5 + bias),rbinom(15, 10,.5 - bias),rbinom(15, 10,.5 + bias))
g2s = 10-g1s

r = 6/length(g1s)
g1p = pbinom(g1s,g1s+g2s,.5) + .05
g2p = pbinom(g2s,g1s+g2s,.5) + .05
g1p = g1p/(g1p+g2p)

mcmc = function(i) {
	mc = rep(NA, length(g1s))
	mc[1] = as.numeric(runif(1) < g1p[1])
	for (i in 2:length(g1s)) {
		pr = abs(mc[i-1] - r)
		p1 = pr * g1p[i]/(pr*g1p[i] + (1-pr)*(1-g1p[i]))
		mc[i] = as.numeric(runif(1) < p1)
	}
	return(mc)
}

res = sapply(1:100, mcmc)
plot(rowMeans(res), type = 'l', col='#999999', yaxs = 'i', ylim = c(0,1))
lines(rowMeans(sapply(1:100, mcmc)), col='#999999')
lines(rowMeans(sapply(1:100, mcmc)), col='#999999')
lines(rowMeans(sapply(1:100, mcmc)), col='#999999')
lines(rowMeans(sapply(1:100, mcmc)), col='#999999')

g1p=g1p[length(g2p):1]
lines(rowMeans(sapply(1:100, mcmc))[length(g2p):1], col='red')
lines(rowMeans(sapply(1:100, mcmc))[length(g2p):1], col='red')
lines(rowMeans(sapply(1:100, mcmc))[length(g2p):1], col='red')
lines(rowMeans(sapply(1:100, mcmc))[length(g2p):1], col='red')
lines(rowMeans(sapply(1:100, mcmc))[length(g2p):1], col='red')

avg = (rowMeans(sapply(1:100, mcmc))[length(g2p):1] + rowMeans(res))/2
lines(avg, lwd=2)

points(g1s/(g1s+g2s), pch = 19)
abline(v=30.5,lty=2, lwd=2)
abline(v=45.5,lty=2, lwd=2)




### 
#recomb chance per gene = 4% chance of recombination for each Mb / 38312 genes over 124 Mb bins (https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0020144#s2)
#1.5 / chr
r = 0.1

MC_chr = function(cell, chr){
gene_use = rownames(g1[which(genes[,1] == chr),])
gene_use = names(which(g1[gene_use,cell]+g2[gene_use,cell] > 2))
g1p = pbinom(g1[gene_use,cell],g1[gene_use,cell]+g2[gene_use,cell],.5) + .05
g2p = pbinom(g2[gene_use,cell],g1[gene_use,cell]+g2[gene_use,cell],.5) + .05
g1p = g1p/(g1p+g2p)


mcmc = function(i) {
	mc = rep(NA, length(g1[gene_use,cell]))
	mc[1] = as.numeric(runif(1) < g1p[1])
	for (i in 2:length(g1[gene_use,cell])) {
		pr = abs(mc[i-1] - r)
		p1 = pr * g1p[i]/(pr*g1p[i] + (1-pr)*(1-g1p[i]))
		mc[i] = as.numeric(runif(1) < p1)
	}
	return(mc)
}

res = sapply(1:100, mcmc)
plot(rowMeans(res), type = 'l', col='#999999', yaxs = 'i', ylim = c(0,1))

g1p=g1p[length(g2p):1]
lines(rowMeans(sapply(1:100, mcmc))[length(g2p):1], col='red')

avg = (rowMeans(sapply(1:100, mcmc))[length(g2p):1] + rowMeans(res))/2
lines(avg, lwd=2)

points(g1[gene_use,cell]/(g1[gene_use,cell]+g2[gene_use,cell]), pch = 19)
}


MC_chr("A288-301_84s", chr=1)
plotCell2('A288-301_84s')

MC_chr("A230-241_1s", chr=1)
plotCell2('A230-241_1s')

calc -log p val per gene (treat as binom distri)
if more likely allele a -> make positive


##################################################
cell = "A194-205_55s"
chr = 4

cell = "A288-301_84s"
chr = 1


HapCall = function(cell, chr, niter = 1000000) {
	keeps = ((g1[,cell] + g2[,cell]) > 0) & (genes$Chr == chr)
	gA = g1[keeps,cell]
	gB = g2[keeps,cell]
	ScoreAB = pbinom(gB,gA+gB,0.5,log.p=T) - pbinom(gA,gA+gB,0.5,log.p=T)
	k = -2*log(1.5/length(ScoreAB))
	ScoreAB[abs(ScoreAB) > k/1.25/2] = sign(ScoreAB[abs(ScoreAB) > k/1.25/2])*k/1.25/2

	blocks = sign(ScoreAB)
	blocks[1] = 1
	for (i in 2:length(ScoreAB)){
		if (sign(ScoreAB[i-1]) == sign(ScoreAB[i])){
			blocks[i] = blocks[i-1]
		} else {
			blocks[i] = blocks[i-1] + 1
		}
	}
	checkBlocks = data.frame(sign_ScoreAB = sign(ScoreAB), blocks = blocks, Position = genes[names(ScoreAB),2])

	blockScore = by(ScoreAB, blocks, sum)
	#blockPos = by(checkBlocks$Position, checkBlocks$blocks, mean)
	#blockLengths = by(ScoreAB, blocks, length)

	#Hap = as.numeric(-sign(blockScore))
	#Hap[which(Hap==0)] = c(1,Hap)[which(Hap==0)]
	Hap = rep(-sign(sum(blockScore)), length(blockScore))

	Scores = c(-sum(blockScore*Hap) - sum(abs(diff(Hap)/2))*k/2, rep(NA,niter))

	for (i in 1:niter) {
    	candidate = sample(1:length(blockScore), 1)
    	Haps = na.omit(Hap[candidate + -1:1])
    	if (length(Haps) == 3) {
    	    recCost = sum(abs(diff(Haps)) - 1)*k
    	} else {
    	    recCost = sum(abs(diff(Haps)) - .5)*k
    	}
    	alpha = blockScore[candidate] * Hap[candidate] + recCost
    	if (alpha >= 0) {
    	    Hap[candidate] = -Hap[candidate]
    	    Scores[i+1] = Scores[i] + alpha
    	} else if (exp(alpha) >= runif(1)) {
    	    Hap[candidate] = -Hap[candidate]
    	    Scores[i+1] = Scores[i] + alpha
    	} else { Scores[i+1] = Scores[i] }
	}

	for (i in 1:length(Hap)){
		checkBlocks[which(checkBlocks$blocks == i),4] = Hap[i]
	}
	colnames(checkBlocks) = c("sign_ScoreAB", "blocks", "Position", "Hap")
	plot(checkBlocks$Position, checkBlocks$Hap, cex=1.5, pch=19)
	keeps2 = which(round(as.numeric(rownames(AlleleFrac))/(10^6)) == chr)
	keeps3 = which(!is.na(AlleleFrac[keeps2,cell]))
	lines(((1:length(keeps2) - 0.5)*10^6)[keeps3], (AlleleFrac[names(keeps3),cell]*2)-1, col='purple2', lwd = 2, lty='dashed')
}


plot(Scores)
plot(Hap, cex=1.5, pch=19)








plotCellChr2 = function(cell, chr, maxs = 20) {
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
	polygon(c(rep(min(checkBlocks$Position[which(checkBlocks$Hap == 1)]), times=2), rep(max(checkBlocks$Position[which(checkBlocks$Hap == 1)]), times=2)), c(-20,20,20,-20), col='mistyrose', border=NA)
	polygon(c(rep(min(checkBlocks$Position[which(checkBlocks$Hap == -1)]), times=2), rep(max(checkBlocks$Position[which(checkBlocks$Hap == -1)]), times=2)), c(-20,20,20,-20), col='lightsteelblue1', border=NA)

	abline(h=0)
	for (ii in which(keeps)) { lines(rep(genes$Position[ii], times=2), c(g1b[ii], -g2b[ii])) }
	points(genes$Position[keeps], -g2a, pch = 19, cex=0.5, col = 'blue')
	points(genes$Position[keeps], g1a, pch = 19, cex=0.5, col = 'red')
    lines(((1:length(keeps2) - 0.5)*10^6)[keeps3], (AlleleFrac[names(keeps3),cell]*40)-20, col='purple2', lwd = 2, lty='dashed')

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
plotCellChr2(cell = 'A194-205_55s', 4)




###original by Brad
niter = 1000000
Hap = -sign(ScoreAB)
Scores = c(-sum(ScoreAB*Hap), rep(NA,niter))

for (i in 1:niter) {
    candidate = sample(1:length(ScoreAB), 1)
    Haps = na.omit(Hap[candidate + -1:1])
    if (length(Haps) == 3) {
        recCost = (sum(abs(diff(Haps))) - 1)*k
    } else {
        recCost = (sum(abs(diff(Haps))) - .5)*k
    }
    alpha = ScoreAB[candidate] * Hap[candidate] + recCost
    if (alpha >= 0) {
        Hap[candidate] = -Hap[candidate]
        Scores[i+1] = Scores[i] + alpha
    } else if (exp(alpha) >= runif(1)) {
        Hap[candidate] = -Hap[candidate]
        Scores[i+1] = Scores[i] + alpha
    } else { Scores[i+1] = Scores[i] }
}
###original by Brad


