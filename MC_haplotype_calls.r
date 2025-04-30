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




avgs = cbind(avgs, avg)
}

plot(rowMeans(avgs),type='l', lwd = 2, yaxs='i')
abline(v=30.5,lty=2, lwd=2)
abline(v=45.5,lty=2, lwd=2)
abline(h=.5)




H1 = sample(0:1)
for (i in 2:100) {
	H1 = c(H1, 






#####################################