#SGT Figures
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

#Make Arabidopsis pollen chromosome 1 plot
#"At_data_2_2026.RData" from "At_haplotype_calls.r" has D, g1, g2, genes, At_Stages, At_meta, HapMatch_all, HapFrac_all, Hap_all, HapFrac_genes_all, No_cell, plotChr, plotScaleBar, plotCell2, HapCallV4, AlleleFrac_genes, HapCallCell
load("At_data_2_2026.RData")
ratio = g1/(g1+g2)
ratio[(g1+g2) < 10] = NA  # Remove measurements with under 10 genoinformative transcripts
geneUse = which(abs(rowMeans(g1/(g1+g2), na.rm=T) - .5) < .4)  # Exclude genes with >90% of all transcripts mapping to the same allele across all samples
ratio[-geneUse,] = NA

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

plotCell1 = function (cell) 
{
    annotate_figure(ggarrange(plotChr(cell, chr = 1), plotScaleBar, ncol = 1, nrow = 2, align = "v",
    heights = c(1, 0.25)), left = text_grob("          % Transcripts from 
    Col-0 allele", 
        rot = 90, size = 10), top = "Arabidopsis pollen")
}

svg('Arabidopsis_plotCell1_MonoTP.svg', width=3, height=2)
plotCell1('A242-253_37s')
dev.off()



#Make Tomato pollen chromosome 1 plot
#"TomatoRemap_11_2024-7_2025.RData" from TomatoRemapping has D2,g1_2,g2_2,Sl_genes2,reads_UMIs2,Sl_Stages,Sl_meta
load("TomatoRemap_11_2024-7_2025_with_meta.RData")
g1 = g1_2
g2 = g2_2
genes = Sl_genes2
ratio = g1/(g1+g2)
ratio[(g1+g2) < 10] = NA  # Remove measurements with under 10 genoinformative transcripts
geneUse = which(abs(rowMeans(g1/(g1+g2), na.rm=T) - .5) < .4)  # Exclude genes with >90% of all transcripts mapping to the same allele across all samples
ratio[-geneUse,] = NA

#plotScaleBar = ggplot() + scale_x_continuous(expand=c(0,0), limits = c(-3, max(genes[,2])/10^6 + 3), breaks = seq(0,30,5)) + theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.x=element_line(), plot.margin = margin(0,0,0,0,'cm')) + xlab('Chromosome position (Mb)')
plotScaleBar = ggplot() + scale_x_continuous(expand=c(0,0), limits = c(-3, max(genes[,2])/10^6 + 3), breaks = seq(0,91,15)) + theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.x=element_line(), plot.margin = margin(0,0,0,0,'cm')) + xlab('Chromosome position (Mb)')

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

plotCell1 = function (cell) 
{
    annotate_figure(ggarrange(plotChr(cell, chr = 1), plotScaleBar, ncol = 1, nrow = 2, align = "v",
    heights = c(1, 0.25)), left = text_grob("          % Transcripts from 
    lycopersicum allele", 
        rot = 90, size = 10), top = "Tomato pollen")
}

svg('Tomato_plotCell1_BM_shed.svg', width=3, height=2)
plotCell1('S80-91_57s')
dev.off()



#Make Rice pollen chromosome 1 plot
#"RiceRemap_7-8_2025.RData" from Rice Remapping has D2,g1_2,g2_2,Rice_genes2,reads_UMIs2
#"Rice.RData" from Prep_for_ortholog_comp2 has Os_meanExp,Os_meanExpL, Os_UMIs, Os_g1, Os_g2, Os_Meta, Os_stages, Os_AlleleFracBin, Os_FracMono, Os_TPM

load("RiceRemap_7-8_2025.RData")
load("Rice.RData")
g1 = g1_2
g2 = g2_2
genes = Rice_genes2
ratio = g1/(g1+g2)
ratio[(g1+g2) < 10] = NA  # Remove measurements with under 10 genoinformative transcripts
geneUse = which(abs(rowMeans(g1/(g1+g2), na.rm=T) - .5) < .4)  # Exclude genes with >90% of all transcripts mapping to the same allele across all samples
ratio[-geneUse,] = NA

#plotScaleBar = ggplot() + scale_x_continuous(expand=c(0,0), limits = c(-3, max(genes[,2])/10^6 + 3), breaks = seq(0,30,5)) + theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.x=element_line(), plot.margin = margin(0,0,0,0,'cm')) + xlab('Chromosome position (Mb)')
plotScaleBar = ggplot() + scale_x_continuous(expand=c(0,0), limits = c(-3, max(genes[,2])/10^6 + 3), breaks = seq(0,44,10)) + theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.x=element_line(), plot.margin = margin(0,0,0,0,'cm')) + xlab('Chromosome position (Mb)')

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

plotCell1 = function (cell) 
{
    annotate_figure(ggarrange(plotChr(cell, chr = 1), plotScaleBar, ncol = 1, nrow = 2, align = "v",
    heights = c(1, 0.25)), left = text_grob("          % Transcripts matching
    reference genome", 
        rot = 90, size = 10), top = "Rice pollen")
}

svg('Rice_plotCell1_TP.svg', width=3, height=2)
plotCell1('R37-90_76s')
dev.off()



#Calculating % of genoinformative transcripts and # of genoinformative genes
mean((colSums(Zm_g1)+colSums(Zm_g2))/colSums(Zm_UMIs))*100 #37.16367% genoinformative transcripts
mean((colSums(Zm2_g1)+colSums(Zm2_g2))/colSums(Zm2_UMIs))*100 #11.36436% genoinformative transcripts
mean((colSums(At_g1)+colSums(At_g2))/colSums(At_UMIs)) #24.38867% genoinformative transcripts
mean((colSums(Os_g1)+colSums(Os_g2))/colSums(Os_UMIs)) #8.955501% genoinformative transcripts
mean((colSums(Sl_g1)+colSums(Sl_g2))/colSums(Sl_UMIs)) #17.61042% genoinformative transcripts

length(which(rowSums(Zm_g1+Zm_g2)>0)) #15814 genoinformative genes
length(which(rowSums(Zm_g1+Zm_g2)>0))/length(rownames(Zm_g1))*100 #66.43421% genoinformative genes
length(which(rowSums(Zm2_g1+Zm2_g2)>0)) #5455
length(which(rowSums(Zm2_g1+Zm2_g2)>0))/length(rownames(Zm2_g1))*100 #22.91632% genoinformative genes
length(which(rowSums(At_g1+At_g2)>0)) #17928
length(which(rowSums(At_g1+At_g2)>0))/length(rownames(At_g1))*100 #46.79474% genoinformative genes
length(which(rowSums(Os_g1+Os_g2)>0)) #6706
length(which(rowSums(Os_g1+Os_g2)>0))/length(rownames(Os_g1))*100 #12.01771% genoinformative genes
length(which(rowSums(Sl_g1+Sl_g2)>0)) #12763
length(which(rowSums(Sl_g1+Sl_g2)>0))/length(rownames(Sl_g1))*100 #38.02813% genoinformative genes




load("Pseudotime_comparison.RData")

library(scales)
pT_summary = data.frame(pT_scaled = c(rescale(Zm_pT_mean, to=c(0,100)), rescale(Os_pT_mean, to=c(0,100)), 
        rescale(Sl_pT_mean, to=c(0,100)), rescale(At_pT_mean, to=c(0,100))), 
    stages = c(names(Zm_pT_mean), names(Os_pT_mean), names(Sl_pT_mean), names(At_pT_mean)), 
    species = c(rep("Zm", length(Zm_pT_mean)), rep("Os", length(Os_pT_mean)), rep("Sl", length(Sl_pT_mean)), rep("At", length(At_pT_mean))))
pT_summary$stages = factor(pT_summary$stages, levels=c("Tetrad", "UM", "Early", "Middle", "Late", "UM/BM", "BM", "BM/TP", "BM_shed", "TP"))

svg('Pseudotime_comparison.svg', width=4, height=3)
ggplot(pT_summary, aes(x = pT_scaled, y = species, color = stages)) +
  geom_point(size = 5, alpha = 0.8) +
  theme_minimal() +
  labs(
    x = "Scaled Pseudotime",
    y = "Species") +
    scale_color_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#F9A729", "Early" = "#F9A729", 
    "Middle" = "#F9A729", "Late" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", 
     "BM_shed" = "#A4A4D5"))
dev.off()

pT_summary2 = data.frame(pT_scaled = c(rescale(Zm_pT_mean[-1], to=c(0,100)), rescale(Os_pT_mean[-1], to=c(0,100)), 
        rescale(Sl_pT_mean[-1], to=c(0,100)), rescale(At_pT_mean[-1], to=c(0,100))), 
    stages = c(names(Zm_pT_mean[-1]), names(Os_pT_mean[-1]), names(Sl_pT_mean[-1]), names(At_pT_mean[-1])), 
    species = c(rep("Zm", length(Zm_pT_mean[-1])), rep("Os", length(Os_pT_mean[-1])), rep("Sl", length(Sl_pT_mean[-1])), rep("At", length(At_pT_mean[-1]))))
pT_summary2$stages = factor(pT_summary$stages, levels=c("UM", "Early", "Middle", "Late", "UM/BM", "BM", "BM/TP", "BM_shed", "TP"))


svg('Pseudotime_comparison2.svg', width=4, height=3)
ggplot(pT_summary2, aes(x = pT_scaled, y = species, color = stages)) +
  geom_point(size = 5, alpha = 0.8) +
  theme_minimal() +
  labs(
    x = "Scaled Pseudotime",
    y = "Species") +
    scale_color_manual(values = c("UM" = "#F9A729", "Early" = "#F9A729", 
    "Middle" = "#F9A729", "Late" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", 
     "BM_shed" = "#A4A4D5"))
dev.off()



library(scales)
pT_summary3 = data.frame(pT_scaled = c(rescale(Zm_pT_mean2, to=c(0,100)), rescale(Os_pT_mean, to=c(0,100)), 
        rescale(Sl_pT_mean, to=c(0,100)), rescale(At_pT_mean2, to=c(0,100))), 
    stages = c(names(Zm_pT_mean), names(Os_pT_mean), names(Sl_pT_mean), names(At_pT_mean2)), 
    species = c(rep("Zm", length(Zm_pT_mean)), rep("Os", length(Os_pT_mean)), rep("Sl", length(Sl_pT_mean)), rep("At", length(At_pT_mean2))))
pT_summary3$stages = factor(pT_summary3$stages, levels=c("Tetrad", "UM", "Early", "Middle", "Late", "UM/BM", "BM", "BM/TP", "BM_shed", "TP"))
pT_summary3$species = factor(pT_summary3$species, levels=c("At", "Sl", "Os", "Zm"))

svg('Pseudotime_comparison3.svg', width=4, height=3)
ggplot(pT_summary3, aes(x = pT_scaled, y = species, color = stages)) +
  geom_point(size = 5, alpha = 0.8) +
  theme_minimal() +
  labs(
    x = "Scaled Pseudotime",
    y = "Species") +
    scale_color_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#F9A729", "Early" = "#F9A729", 
    "Middle" = "#F9A729", "Late" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", 
     "BM_shed" = "#A4A4D5"))
dev.off()




library(scales)
pT_summary4 = data.frame(pT_scaled = c(rescale(Zm_pT_stage2$pT, to=c(0,100)), rescale(Os_pT_stage$pT, to=c(0,100)), 
        rescale(Sl_pT_stage$pT, to=c(0,100)), rescale(At_pT_stage2$pT, to=c(0,100))), 
    stages = c(Zm_pT_stage2$stages, Os_pT_stage$stages, Sl_pT_stage$stages, At_pT_stage2$stages), 
    species = c(rep("Zm", length(Zm_pT_stage2$stages)), rep("Os", length(Os_pT_stage$stages)), rep("Sl", length(Sl_pT_stage$stages)), rep("At", length(At_pT_stage2$stages))))
pT_summary4$stages = factor(pT_summary4$stages, levels=c("Tetrad", "UM", "UM Bi", "UM Mono", "Early", "Middle", "Late", "UM/BM", "BM", "BM/TP", "BM_shed", "TP", "TP_shed"))
pT_summary4$species = factor(pT_summary4$species, levels=c("At", "Sl", "Os", "Zm"))

svg('Pseudotime_comparison4.svg', width=4, height=3)
ggplot(pT_summary4, aes(x = pT_scaled, y = species, color = stages)) +
  geom_jitter(size = 0.5, alpha = 0.8, height=0.25) +
  theme_minimal() + #geom_jitter() + 
  labs(
    x = "Scaled Pseudotime",
    y = "Species") +
    scale_color_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#F9A729", "UM Bi" = "#F9A729", "UM Mono" = "#F9A729", 
    "Early" = "#F9A729", "Middle" = "#F9A729", "Late" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "TP" = "#A4A4D5", 
     "BM_shed" = "#A4A4D5", "TP_shed" = "#A4A4D5"))
dev.off()



#pseudo plot, merged Zm, no alpha, new At cutoff
load("Arabidopsis.RData")
load("Rice.RData")
load("Maize_merge.RData")
load("Tomato.RData") 
#load("Pseudotime_comparison2.RData") #Os_pT_stage, Zm_merge_pT_stage, Sl_pT_stage, At_pT_stage2
load("Pseudotime_comparison2_At10k.RData") #Os_pT_stage, Zm_merge_pT_stage, Sl_pT_stage, At_pT_stage2

library(scales)
pT_summary5 = data.frame(pT_scaled = c(rescale(Zm_merge_pT_stage$pT, to=c(0,100)), rescale(Os_pT_stage$pT, to=c(0,100)), 
        rescale(Sl_pT_stage$pT, to=c(0,100)), rescale(At_pT_stage2$pT, to=c(0,100))), 
    stages = c(Zm_merge_pT_stage$stages, Os_pT_stage$stages, Sl_pT_stage$stages, At_pT_stage2$stages), 
    species = c(rep("Zm", length(Zm_merge_pT_stage$stages)), rep("Os", length(Os_pT_stage$stages)), rep("Sl", length(Sl_pT_stage$stages)), rep("At", length(At_pT_stage2$stages))))
pT_summary5$stages[grep("UM ", pT_summary5$stages)] = "UM"
pT_summary5$stages[which(pT_summary5$stages == "TP" | pT_summary5$stages == "BM_shed" | pT_summary5$stages == "TP_shed")] = "Pollen"
pT_summary5$stages = factor(pT_summary5$stages, levels=c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "Pollen"))
pT_summary5$species = factor(pT_summary5$species, levels=c("At", "Sl", "Os", "Zm"))                                
rownames(pT_summary5) = c(rownames(Zm_merge_pT_stage), rownames(Os_pT_stage), rownames(Sl_pT_stage), rownames(At_pT_stage2))
pT_summary5$hapfrac = c(Zm_merge_FracMono, Os_FracMono, Sl_FracMono, At_HapFrac)[rownames(pT_summary5)]

svg('Pseudotime_comparison5.svg', width=4, height=3)
ggplot(pT_summary5, aes(x = pT_scaled, y = species, color = stages)) +
  geom_jitter(size = 0.5, height=0.25) +
  theme_minimal() + #geom_jitter() + 
  labs(
    x = "Scaled Pseudotime",
    y = "Species") +
    scale_color_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#F9A729", 
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "Pollen" = "#A4A4D5"))
dev.off()

svg('Pseudotime_comparison5_.svg', width=4, height=3)
ggplot(pT_summary5, aes(x = pT_scaled, y = species, color = stages)) +
  geom_jitter(size = 0.5, height=0.25) +
  theme_minimal() + #geom_jitter() + 
  labs(
    x = "Scaled Pseudotime",
    y = "Species") +
    scale_color_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#F9A729", 
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "Pollen" = "#A4A4D5"))
dev.off()

FracMono_col = colorRamp2(c(0, 100), c("#C5BFBE", "purple4"))
FracMono_col = colorRamp2(c(0, 100), c("#C5BFBE", "#8074A8"))

svg('Pseudotime_comparison5_hapfrac.svg', width=4, height=3)
ggplot(pT_summary5, aes(x = pT_scaled, y = species, color = hapfrac)) +
  geom_jitter(size = 0.5, height=0.25) +
  theme_minimal() + #geom_jitter() + 
  labs(
    x = "Scaled Pseudotime",
    y = "Species") + scale_color_gradient(low = "#999", high = "#B10318")
dev.off()
#B10318
pT_summary5$hapfrac = c(rescale(Zm_merge_FracMono, to=c(0,100)), rescale(Os_FracMono, to=c(0,100)), 
    rescale(Sl_FracMono, to=c(0,100)), rescale(At_HapFrac, to=c(0,100)))[rownames(pT_summary5)]

svg('Pseudotime_comparison5_hapfrac_rescaled.svg', width=4, height=3)
ggplot(pT_summary5, aes(x = pT_scaled, y = species, color = hapfrac)) +
  geom_jitter(size = 0.5, height=0.25) +
  theme_minimal() + #geom_jitter() + 
  labs(
    x = "Scaled Pseudotime",
    y = "Species") + scale_color_gradient(low = "#C5BFBE", high = "purple4")
dev.off()


#Using At 10k cutoff
load("Arabidopsis.RData")
load("Rice.RData")
load("Maize_merge.RData")
load("Tomato.RData") 
load("Pseudotime_comparison2_At10k.RData") #Os_pT_stage, Zm_merge_pT_stage, Sl_pT_stage, At_pT_stage2

library(scales)
pT_summary5 = data.frame(pT_scaled = c(rescale(Zm_merge_pT_stage$pT, to=c(0,100)), rescale(Os_pT_stage$pT, to=c(0,100)), 
        rescale(Sl_pT_stage$pT, to=c(0,100)), rescale(At_pT_stage2$pT, to=c(0,100))), 
    stages = c(Zm_merge_pT_stage$stages, Os_pT_stage$stages, Sl_pT_stage$stages, At_pT_stage2$stages), 
    species = c(rep("Zm", length(Zm_merge_pT_stage$stages)), rep("Os", length(Os_pT_stage$stages)), rep("Sl", length(Sl_pT_stage$stages)), rep("At", length(At_pT_stage2$stages))))
pT_summary5$stages[grep("UM ", pT_summary5$stages)] = "UM"
pT_summary5$stages[which(pT_summary5$stages == "TP" | pT_summary5$stages == "BM_shed" | pT_summary5$stages == "TP_shed")] = "Pollen"
pT_summary5$stages = factor(pT_summary5$stages, levels=c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "Pollen"))
pT_summary5$species = factor(pT_summary5$species, levels=c("At", "Sl", "Os", "Zm"))                                
rownames(pT_summary5) = c(rownames(Zm_merge_pT_stage), rownames(Os_pT_stage), rownames(Sl_pT_stage), rownames(At_pT_stage2))
pT_summary5$hapfrac = c(Zm_merge_FracMono, Os_FracMono, Sl_FracMono, At_HapFrac)[rownames(pT_summary5)]

randomSamp = sample(length(pT_summary5$stages))


set.seed(7)
svg('Pseudotime_comparison5_At10k.svg', width=4, height=3)
ggplot(pT_summary5[randomSamp,], aes(x = pT_scaled, y = species, color = stages)) +
  geom_jitter(size = 0.5, height=0.25) +
  theme_minimal() + 
  labs(
    x = "Scaled Pseudotime",
    y = "Species") +
    scale_color_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#F9A729", 
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "Pollen" = "#A4A4D5"))
dev.off()

set.seed(7)
svg('Pseudotime_comparison5_hapfrac_At10k.svg', width=4, height=3)
ggplot(pT_summary5[randomSamp,], aes(x = pT_scaled, y = species, color = hapfrac)) +
  geom_jitter(size = 0.5, height=0.25) +
  theme_minimal() + 
  labs(
    x = "Scaled Pseudotime",
    y = "Species") + scale_color_gradient(low = "#999", high = "#B10318")
dev.off()



#Checking whether stage or hapfrac separate Arabidopsis UM and BM samples more
pT_summary5[which(pT_summary5$stages == "UM" & pT_summary5$species == "At" | pT_summary5$stages == "BM" & pT_summary5$species == "At"),]
AtUMBM = pT_summary5[which(pT_summary5$stages == "UM" & pT_summary5$species == "At" | pT_summary5$stages == "BM" & pT_summary5$species == "At"),]
randomSamp2 = sample(length(AtUMBM$stages))
AtUMBM$stages = factor(AtUMBM$stages, levels = c("UM", "BM"))

set.seed(7)
svg('Pseudotime_comparison5_At10k_UMBM.svg', width=4, height=3)
ggplot(AtUMBM[randomSamp2,], aes(x = pT_scaled, y = species, color = stages)) +
  geom_jitter(size = 0.5, height=0.25) +
  theme_minimal() + 
  labs(
    x = "Scaled Pseudotime",
    y = "Species") +
    scale_color_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#F9A729", 
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "Pollen" = "#A4A4D5"))
dev.off()

set.seed(7)
svg('Pseudotime_comparison5_hapfrac_At10k_UMBM.svg', width=4, height=3)
ggplot(AtUMBM[randomSamp2,], aes(x = pT_scaled, y = species, color = hapfrac)) +
  geom_jitter(size = 0.5, height=0.25) +
  theme_minimal() + 
  labs(
    x = "Scaled Pseudotime",
    y = "Species") + scale_color_gradient(low = "#999", high = "#B10318")
dev.off()

set.seed(7)
svg('Pseudotime_comparison5_hapfrac_At10k_UM.svg', width=4, height=3)
ggplot(AtUMBM[AtUMBM$stages == "UM",], aes(x = pT_scaled, y = species, color = hapfrac)) +
  geom_jitter(size = 0.5, height=0.25) +
  theme_minimal() + xlim(0, 90) +
  labs(
    x = "Scaled Pseudotime",
    y = "Species") + scale_color_gradient(low = "#999", high = "#B10318")
dev.off()


set.seed(7)
svg('Pseudotime_comparison5_hapfrac_At10k_BM.svg', width=4, height=3)
ggplot(AtUMBM[AtUMBM$stages == "BM",], aes(x = pT_scaled, y = species, color = hapfrac)) +
  geom_jitter(size = 0.5, height=0.25) +
  theme_minimal() + xlim(0, 90) +
  labs(
    x = "Scaled Pseudotime",
    y = "Species") + scale_color_gradient(low = "#999", high = "#B10318")
dev.off()




#####Making two color stacked barplots with % haploid vs. diploid expressing samples at each stage in each species
load('Rice.RData')
load("Tomato.RData")
load("Arabidopsis.RData")
load("Maize.RData")

hist(Zm_FracMono) #threshold of 40
hist(Zm2_FracMono) #threshold of 40
hist(Os_FracMono) #threshold of 40
hist(Sl_FracMono) #threshold of 40
hist(At_HapFrac) #threshold of 40

Zm_merge_stages = c(as.character(Zm_stages), as.character(Zm2_stages))
names(Zm_merge_stages) = c(names(Zm_stages), names(Zm2_stages))
Zm_merge_stages[rownames(Zm_Meta[grep("mix", Zm_Meta[,5]),])] = "UM/BM"
Zm_merge_stages[which(Zm_merge_stages == "Early" | Zm_merge_stages == "Middle" | Zm_merge_stages == "Late")] = "UM"
Zm_merge_stages[which(Zm_merge_stages == "TP_shed" | Zm_merge_stages == "TP")] = "Pollen"
Zm_merge_stages = factor(Zm_merge_stages, levels = c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "Pollen"))

stages = c(Zm_merge_stages, Os_stages, Sl_stages, At_stages)
stages = sub("Early", "UM", sub("Middle", "UM", sub("Late", "UM", sub("UM .*", "UM", sub(".*_shed", "Shed", stages)))))
stages[which(stages == "TP" | stages == "Shed")] = "Pollen"
stages = factor(stages, levels = c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "Pollen"))
names(stages) = c(names(Zm_merge_stages), names(Os_stages), names(Sl_stages), names(At_stages))
species = sub("L", "Zm", sub("M", "Zm", sub("R", "Os", sub("S", "Sl", sub("A", "At", substr(c(colnames(Zm_TPM), colnames(Zm2_TPM), colnames(Os_TPM), colnames(Sl_TPM), colnames(At_TPM)), 1,1))))))

FracMono = c(Zm_FracMono, Zm2_FracMono, Os_FracMono, Sl_FracMono, At_HapFrac)

FracMono_data = data.frame(stages = stages, species = species, FracMono = FracMono, TranscriptOrigin = FracMono > 40)
FracMono_data$TranscriptOrigin[FracMono_data$TranscriptOrigin == "TRUE"] = "Haploid"
FracMono_data$TranscriptOrigin[FracMono_data$TranscriptOrigin == "FALSE"] = "Diploid"

FracMono_count <- FracMono_data %>%
  count(species, stages, TranscriptOrigin, name = "count")
FracMono_count = FracMono_count[-grep("BM/TP", FracMono_count$stages),]

FracPlot = function(sp){ggplot(FracMono_count[FracMono_count$species == sp,], aes(fill=TranscriptOrigin, y=count, x=stages)) + 
    geom_bar(position="fill", stat="identity") + theme_classic() + 
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25)) +
    labs(title = sp, x = "Stage", y = "% of pollen percursors 
    or pollen") + 
    scale_fill_manual(values = c("Haploid" = "#8074A8", "Diploid" = "#C5BFBE")) 
}


#FracPlot = function(sp){ggplot(FracMono_data[grep(sp, FracMono_data$species),], aes(fill=Transcripts, y=1, x=stages)) + 
#    geom_bar(position="fill", stat="identity") + theme_classic() + labs(title = "Maize", x = "Stage", 
#    y = "% of pollen precursors or pollen") + scale_y_continuous(limits = c(0, 1)) + geom_col(color = NA) +
#    scale_fill_manual(values = c("Haploid" = "#8074A8", "Diploid" = "#C5BFBE")) 
#}

ggsave("FracMono_barZm.svg", FracPlot("Zm"), width = 4, height = 1.8, units = "in")
ggsave("FracMono_barOs.svg", FracPlot("Os"), width = 4, height = 1.8, units = "in")
ggsave("FracMono_barSl.svg", FracPlot("Sl"), width = 4, height = 1.8, units = "in")
ggsave("FracMono_barAt.svg", FracPlot("At"), width = 4, height = 1.8, units = "in") #need to rerun, accidentally overwrote


#####Remaking stacked barplots with % haploid expressing samples at each stage in each species w/o UM/BM and w/ more bins

Zm_merge_stages = c(as.character(Zm_stages), as.character(Zm2_stages))
names(Zm_merge_stages) = c(names(Zm_stages), names(Zm2_stages))
#Zm_merge_stages[rownames(Zm_Meta[grep("mix", Zm_Meta[,5]),])] = "UM/BM"
Zm_merge_stages[which(Zm_merge_stages == "Early" | Zm_merge_stages == "Middle" | Zm_merge_stages == "Late")] = "UM"
Zm_merge_stages[which(Zm_merge_stages == "TP_shed" | Zm_merge_stages == "TP")] = "Pollen"
Zm_merge_stages = factor(Zm_merge_stages, levels = c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "Pollen"))

stages = c(Zm_merge_stages, Os_stages, Sl_stages, At_stages)
stages = sub("Early", "UM", sub("Middle", "UM", sub("Late", "UM", sub("UM .*", "UM", sub(".*_shed", "Shed", stages)))))
stages[which(stages == "TP" | stages == "Shed")] = "Pollen"
stages = factor(stages, levels = c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "Pollen"))
names(stages) = c(names(Zm_merge_stages), names(Os_stages), names(Sl_stages), names(At_stages))
species = sub("L", "Zm", sub("M", "Zm", sub("R", "Os", sub("S", "Sl", sub("A", "At", substr(c(colnames(Zm_TPM), colnames(Zm2_TPM), colnames(Os_TPM), colnames(Sl_TPM), colnames(At_TPM)), 1,1))))))

FracMono = c(Zm_FracMono, Zm2_FracMono, Os_FracMono, Sl_FracMono, At_HapFrac)


# Define bins
bin_FracMono <- function(x) {
  case_when(
    is.na(x)                 ~ NA_character_,
    x < 20                   ~ "0-20",
    x >= 20 & x < 40             ~ "20-40",
    x >= 40 & x < 60           ~ "40-60",
    x >= 60 & x < 80          ~ "60-80",
    x >= 80                 ~ "80-100",
  )
}

FracMono_bins = factor(bin_FracMono(FracMono), levels=c("0-20", "20-40", "40-60", "60-80", "80-100"))

FracMono_data = data.frame(stages = stages, species = species, FracMono = FracMono, FracMono_bins=FracMono_bins)

FracMono_count <- FracMono_data %>%
  count(species, stages, FracMono_bins, name = "count")
FracMono_count = FracMono_count[-grep("BM/TP", FracMono_count$stages),]
FracMono_count = FracMono_count[-grep("UM/BM", FracMono_count$stages),]


c("0-20" = "#CDCBDB", "20-40" = "#B4B1CC", "40-60" = "#9995BC", "60-80" = "#7F78AC", "80-100" = "#645A9D")
"#CDCBDB", 
"#B4B1CC",
"#9995BC",
"#7F78AC",
"#645A9D"
FracPlot = function(sp){ggplot(FracMono_count[FracMono_count$species == sp,], aes(fill=FracMono_bins, y=count, x=stages)) + 
    geom_bar(position="fill", stat="identity") + theme_classic() + 
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25)) +
    labs(title = sp, x = "Stage", y = "% of pollen percursors 
    or pollen") + 
    scale_fill_manual(values = c("0-20" = "#CDCBDB", "20-40" = "#B4B1CC", "40-60" = "#9995BC", 
    "60-80" = "#7F78AC", "80-100" = "#645A9D")) 
}


#FracPlot = function(sp){ggplot(FracMono_data[grep(sp, FracMono_data$species),], aes(fill=Transcripts, y=1, x=stages)) + 
#    geom_bar(position="fill", stat="identity") + theme_classic() + labs(title = "Maize", x = "Stage", 
#    y = "% of pollen precursors or pollen") + scale_y_continuous(limits = c(0, 1)) + geom_col(color = NA) +
#    scale_fill_manual(values = c("Haploid" = "#8074A8", "Diploid" = "#C5BFBE")) 
#}

ggsave("FracMono_barZm_bins.svg", FracPlot("Zm"), width = 4, height = 1.8, units = "in")
ggsave("FracMono_barOs_bins.svg", FracPlot("Os"), width = 4, height = 1.8, units = "in")
ggsave("FracMono_barSl_bins.svg", FracPlot("Sl"), width = 4, height = 1.8, units = "in")
ggsave("FracMono_barAt_bins.svg", FracPlot("At"), width = 4, height = 1.8, units = "in") #need to rerun, accidentally overwrote




#Same type of stacked barplots but for microscopy data
AtDsRed_OG <- read_excel("C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R/AtDsRed_9-24-26.xlsx")
SlGFP_OG <- read_excel("C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R/SlGFP_9-24-26.xlsx")

AtDsRed = as.data.frame(AtDsRed_OG[,c(9,18,20)])
colnames(AtDsRed) = c("stages", "haploid", "total")
AtDsRed_count = data.frame(haploid = do.call(rbind, as.list(by(AtDsRed$haploid, AtDsRed$stages, sum))),
    total = do.call(rbind, as.list(by(AtDsRed$total, AtDsRed$stages, sum))))
AtDsRed_count$stages = rownames(AtDsRed_count)
AtDsRed_count = AtDsRed_count[-c(grep("BM/TP", AtDsRed_count$stages), grep("UM/BM", AtDsRed_count$stages)),]
AtDsRed_count$stages = factor(AtDsRed_count$stages, levels = c("Tetrad", "UM", "BM", "Pollen"))
AtDsRed_count$diploid = AtDsRed_count$total - AtDsRed_count$haploid
AtDsRed_count = AtDsRed_count[levels(AtDsRed_count$stages),]
AtDsRed_data = data.frame(TranscriptOrigin = c(rep("Diploid", times = length(AtDsRed_count$stages)), 
    rep("Haploid", times = length(AtDsRed_count$stages))), 
    count = c(AtDsRed_count$diploid, AtDsRed_count$haploid), stages = rep(AtDsRed_count$stages, times=2))

SlGFP = as.data.frame(SlGFP_OG[,c(4,5,6)])
colnames(SlGFP) = c("stages", "haploid", "total")
SlGFP_count = data.frame(haploid = do.call(rbind, as.list(by(SlGFP$haploid, SlGFP$stages, sum))),
    total = do.call(rbind, as.list(by(SlGFP$total, SlGFP$stages, sum))))
SlGFP_count$stages = rownames(SlGFP_count)
SlGFP_count$stages = factor(SlGFP_count$stages, levels = c("Tetrad", "UM", "BM", "Pollen"))
SlGFP_count$diploid = SlGFP_count$total - SlGFP_count$haploid
SlGFP_count = SlGFP_count[levels(SlGFP_count$stages),]
SlGFP_data = data.frame(TranscriptOrigin = c(rep("Diploid", times = length(SlGFP_count$stages)), 
    rep("Haploid", times = length(SlGFP_count$stages))), 
    count = c(SlGFP_count$diploid, SlGFP_count$haploid), stages = rep(SlGFP_count$stages, times=2))

FracPlotLAT52 = function(data, sp){ggplot(data, aes(fill=TranscriptOrigin, y=count, x=stages)) + 
    geom_bar(position="fill", stat="identity") + theme_classic() + 
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25)) +
    labs(title = sp, x = "Stage", y = "% of pollen percursors 
    or pollen") + 
    scale_fill_manual(values = c("Haploid" = "#8074A8", "Diploid" = "#C5BFBE")) 
}

ggsave("FracMonoLAT52_barAt.svg", FracPlotLAT52(AtDsRed_data, "Arabidopsis LAT52:DsRed qrt1"), width = 4, height = 1.8, units = "in")

ggsave("FracMonoLAT52_barSl.svg", FracPlotLAT52(SlGFP_data, "Tomato LAT52:GFP"), width = 4, height = 1.8, units = "in")





GN = list(Diploid = names(which(rowMeans(Zm_merge_TPM[,names(which(Zm_merge_FracMono <= 40))] > 0) > 0)), 
    Haploid = names(which(rowMeans(Zm_merge_TPM[,names(which(Zm_merge_FracMono > 40))] > 0) > 0)))
expC = matrix('-', ncol = 2, nrow = length(unique(unlist(GN))))
rownames(expC) = unique(unlist(GN))
colnames(expC) = names(GN)
for (i in 1:length(GN)) {
	expC[GN[[i]],i] = '+'
}
expC = apply(expC,1,paste,collapse=':')
expC[expC == '-:+:-'] = '+:+:-'

ggplot(data.frame(x=table(expC)), aes(x="", y=x.Freq, fill=x.expC)) + 
geom_bar(stat="identity", width=1, color="white") + coord_polar("y", start=0) + theme_void() + 
scale_fill_manual(values=c('#d4aa00','#33a02c','#1f78b4'))#,'#144d74'))
Zm_expC = data.frame(x=table(expC), species = rep("Zm", times=3))

length(GN$Haploid)/length(unique(unlist(GN)))*100 #89.99826% haploid expressed



GN = list(Diploid = names(which(rowMeans(Os_TPM[,names(which(Os_FracMono <= 40))] > 0) > 0)), 
    Haploid = names(which(rowMeans(Os_TPM[,names(which(Os_FracMono > 40))] > 0) > 0)))
expC = matrix('-', ncol = 2, nrow = length(unique(unlist(GN))))
rownames(expC) = unique(unlist(GN))
colnames(expC) = names(GN)
for (i in 1:length(GN)) {
	expC[GN[[i]],i] = '+'
}
expC = apply(expC,1,paste,collapse=':')
expC[expC == '-:+:-'] = '+:+:-'

ggplot(data.frame(x=table(expC)), aes(x="", y=x.Freq, fill=x.expC)) + 
geom_bar(stat="identity", width=1, color="white") + coord_polar("y", start=0) + theme_void() + 
scale_fill_manual(values=c('#d4aa00','#33a02c','#1f78b4'))#,'#144d74'))
Os_expC = data.frame(x=table(expC), species = rep("Os", times=3))

length(GN$Haploid)/length(unique(unlist(GN)))*100 #91.3691% haploid expressed



GN = list(Diploid = names(which(rowMeans(Sl_TPM[,names(which(Sl_FracMono <= 40))] > 0) > 0)), 
    Haploid = names(which(rowMeans(Sl_TPM[,names(which(Sl_FracMono > 40))] > 0) > 0)))
expC = matrix('-', ncol = 2, nrow = length(unique(unlist(GN))))
rownames(expC) = unique(unlist(GN))
colnames(expC) = names(GN)
for (i in 1:length(GN)) {
	expC[GN[[i]],i] = '+'
}
expC = apply(expC,1,paste,collapse=':')
expC[expC == '-:+:-'] = '+:+:-'

ggplot(data.frame(x=table(expC)), aes(x="", y=x.Freq, fill=x.expC)) + 
geom_bar(stat="identity", width=1, color="white") + coord_polar("y", start=0) + theme_void() + 
scale_fill_manual(values=c('#d4aa00','#33a02c','#1f78b4'))#,'#144d74'))
Sl_expC = data.frame(x=table(expC), species = rep("Sl", times=3))


length(GN$Haploid)/length(unique(unlist(GN)))*100 #94.28736% haploid expressed



GN = list(Diploid = names(which(rowMeans(At_TPM[,names(which(At_HapFrac <= 40))] > 0) > 0)), 
    Haploid = names(which(rowMeans(At_TPM[,names(which(At_HapFrac > 40))] > 0) > 0)))
expC = matrix('-', ncol = 2, nrow = length(unique(unlist(GN))))
rownames(expC) = unique(unlist(GN))
colnames(expC) = names(GN)
for (i in 1:length(GN)) {
	expC[GN[[i]],i] = '+'
}
expC = apply(expC,1,paste,collapse=':')
expC[expC == '-:+:-'] = '+:+:-'

ggplot(data.frame(x=table(expC)), aes(x="", y=x.Freq, fill=x.expC)) + 
geom_bar(stat="identity", width=1, color="white") + coord_polar("y", start=0) + theme_void() + 
scale_fill_manual(values=c('#d4aa00','#33a02c','#1f78b4'))#,'#144d74'))
At_expC = data.frame(x=table(expC), species = rep("At", times=3))

length(GN$Haploid)/length(unique(unlist(GN)))*100 #91.23552% haploid expressed


expC = rbind(Zm_expC, Os_expC, Sl_expC, At_expC)


HapExp = ggplot(expC, aes(fill=x.expC, y=x.Freq, x=species)) + 
    geom_bar(position="fill", stat="identity") + theme_classic() + 
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = seq(0, 100, by = 25)) +
    labs(x = "Species", y = "% of genes
    expressed") + scale_fill_manual(values=c('#d4aa00','#33a02c','#1f78b4'))


ggsave("HapExpBar.svg", HapExp, width = 4, height = 1.8, units = "in")



###Making example pollen precursor and pollen grain histograms

CalcGeneFrac = function(g1, g2){
    ratio = g1/(g1+g2)
    ratio[(g1+g2) < 10] = NA  # Remove measurements with under 10 genoinformative transcripts
    geneUse = which(abs(rowMeans(g1/(g1+g2), na.rm=T) - .5) < .4)  # Exclude genes with >90% of all transcripts mapping to the same allele across all samples
    ratio[-geneUse,] = NA
    geneFrac = ratio*100
    return(geneFrac)
}
Zm_merge_geneFrac = CalcGeneFrac(Zm_merge_g1, Zm_merge_g2)
Os_geneFrac = CalcGeneFrac(Os_g1, Os_g2)
Sl_geneFrac = CalcGeneFrac(Sl_g1, Sl_g2)
At_geneFrac = CalcGeneFrac(At_g1, At_g2)

ratioHist = function(geneFrac, cell, allele){
    gghistogram(data.frame(geneFrac = geneFrac[,cell]), x = 'geneFrac', fill = 'black', bins = 15) + 
    theme(strip.background = element_blank(), strip.text.x = element_blank(), legend.position = 'none') + 
    scale_y_continuous(expand = c(0,0)) + xlab(paste('% Transcripts matching
    ', allele, sep = " ")) + ylab('Genes') + 
    scale_x_continuous(breaks = seq(0,100,10), labels=c('0', rep('',4), '50', rep('',4), '100')) +
    coord_cartesian(xlim = c(0, 100))
}


svg('Maize_ratioHist_tetrad.svg', width=2.5, height=2)
ratioHist(Zm_merge_geneFrac, 'Lane1_RPI3_11s', 'B73 allele') #use to pick tetrad (sep): colSums(Zm_merge_UMIs[,names(Zm_merge_stages[Zm_merge_stages == "Tetrad"])])
dev.off()

svg('Maize_ratioHist_pollen.svg', width=2.5, height=2)
ratioHist(Zm_merge_geneFrac, 'Lane1_RPI6_77s', 'B73 allele') #use to pick pollen: colSums(Zm_merge_UMIs[,names(Zm_merge_stages[Zm_merge_stages == "Pollen"])])
dev.off()

svg('Rice_ratioHist_intact_tetrad.svg', width=2.5, height=2)
ratioHist(Os_geneFrac, 'R119-136_75s', 'reference') #use to pick intact tetrad: colSums(Os_UMIs[,names(Os_stages[which(Os_Meta[,12] == "tetrads")])])
dev.off()

svg('Rice_ratioHist_pollen.svg', width=2.5, height=2)
ratioHist(Os_geneFrac, 'R125-142_75s', 'reference') #use to pick pollen: colSums(Os_UMIs[,names(Os_stages[Os_stages == "TP"])])
dev.off()

ratioHist(Sl_geneFrac, 'S80-91_1s', 'lycopersicum allele')
ratioHist(Sl_geneFrac, 'S80-91_3s', 'lycopersicum allele')
ratioHist(Sl_geneFrac, 'S92-103_11s', 'lycopersicum allele')
ratioHist(Sl_geneFrac, 'S92-103_14s', 'lycopersicum allele')

svg('Tomato_ratioHist_intact_tetrad.svg', width=2.5, height=2)
ratioHist(Sl_geneFrac, 'S80-91_1s', 'lycopersicum allele') #use to pick intact tetrad: colSums(Sl_UMIs[,names(Sl_stages[which(Sl_Meta[,8] == "tetrad")])])
dev.off()

svg('Tomato_ratioHist_pollen.svg', width=2.5, height=2)
ratioHist(Sl_geneFrac, 'S26-37_27s', 'lycopersicum allele') #use to pick pollen: colSums(Sl_UMIs[,names(Sl_stages[Sl_stages == "BM_shed"])])
dev.off()


svg('Arabidopsis_ratioHist_intact_tetrad.svg', width=2.5, height=2)
ratioHist(At_geneFrac, 'A268-287_91s', 'Col-0 allele') #use to pick intact tetrad: colSums(At_UMIs[,names(At_stages[At_stages == "Tetrad"])])
dev.off()

svg('Arabidopsis_ratioHist_pollen.svg', width=2.5, height=2)
ratioHist(At_geneFrac, 'A194-205_76s', 'Col-0 allele') #use to pick pollen: colSums(At_UMIs[,names(At_stages[At_stages == "TP"])])
dev.off()








> table(Os_stages)
Os_stages
 Tetrad   UM Bi UM Mono   UM/BM      BM   BM/TP      TP TP_shed 
     33      70      33      10      24      10     124       3
mean(colSums(Os_UMIs[,names(Os_stages[which(Os_Meta[,12] == "tetrad_sep")])])) #16 sep tetrads
[1] 22056.06
mean(colSums(Os_UMIs[, names(Os_stages[which(Os_stages == "UM Bi")])]))
[1] 14281.44
> mean(colSums(Os_UMIs[, names(Os_stages[which(Os_stages == "UM Mono")])]))
[1] 72989.42
> mean(colSums(Os_UMIs[, names(Os_stages[which(Os_stages == "UM/BM")])]))  
[1] 87327.3
> mean(colSums(Os_UMIs[, names(Os_stages[which(Os_stages == "BM")])]))   
[1] 122038.5
> mean(colSums(Os_UMIs[, names(Os_stages[which(Os_stages == "BM/TP")])]))
[1] 73317.7
> mean(colSums(Os_UMIs[, names(Os_stages[which(Os_stages == "TP")])]))   
[1] 91745.29
> mean(colSums(Os_UMIs[, names(Os_stages[which(Os_stages == "TP_shed")])]))
[1] 16060.33

calcFC = function(UMIs, stages, pre, post){
out = log(UMIs + 0.5, 2)
out = rowMeans(out[,stages == post]) - rowMeans(out[,stages == pre])
out = data.frame(FC = out, Cat = c('Up', 'NC', 'Down')[(out <= -1) + (out < 1) + 1])
out = out[rowMeans(UMIs[,stages %in% c(post, pre)]) >= 1,]
ggplot(out, aes(x=FC, fill=Cat)) + geom_histogram(binwidth=2/15) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab('Genes') + xlab('log2 fold change') + scale_x_continuous(lim=c(-4,6), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#888888", "#AAAAAA", "#888888")) + geom_vline(xintercept=0, linetype='dashed') + geom_vline(xintercept=1) + geom_vline(xintercept=-1)
}
calcFC(Os_UMIs, Os_stages, "UM Bi", "UM Mono")




bootFn = function(x, trim = 0.2, nboot = 2000) {
	out = sapply(1:nboot, function(i) { mean(sample(x, replace = T), trim = trim, na.rm=T) })
	return(c(mean(out), mean(out) - sd(out), mean(out) + sd(out)))
}

repFig4 = function(UMIs, stages, pre, post, g1, g2){
TranscriptsByStage = data.frame(t(matrix(unlist(by(colSums(UMIs)/10^3, stages, bootFn)), nrow = 3)))
colnames(TranscriptsByStage) = c('val', 'lower', 'upper')
TranscriptsByStage$Stage = unique(stages[order(stages)])
Fig4a = ggplot(TranscriptsByStage) + geom_bar(aes(x=Stage, y=val), stat='identity', fill='gray') + geom_errorbar(aes(x=Stage, ymin=lower, ymax=upper), width=0.3, size=.9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90),  axis.title.x = element_blank()) + scale_y_continuous(expand=c(0,0)) + ylab('Transcripts per\npollen precursor (x1000)')

## Fig 4b
out = log(UMIs + 0.5, 2)
out = rowMeans(out[,stages == post]) - rowMeans(out[,stages == pre])
out = data.frame(FC = out, Cat = c('Up', 'NC', 'Down')[(out <= -1) + (out < 1) + 1])
out = out[rowMeans(UMIs[,stages %in% c(post, pre)]) >= 1,]
Fig4b = ggplot(out, aes(x=FC, fill=Cat)) + geom_histogram(binwidth=2/15) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab('Genes') + xlab('log2 fold change') + scale_x_continuous(lim=c(-4,6), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#888888", "#AAAAAA", "#888888")) + geom_vline(xintercept=0, linetype='dashed') + geom_vline(xintercept=1) + geom_vline(xintercept=-1)

## Fig 4c
gMin = (g1 + g2 - abs(g1-g2))/2
gMax = g1+g2-gMin
gMin2=gMin
gMax2=gMax
gMin2[(g1+g2) < 10] = NA
gMax2[(g1+g2) < 10] = NA

use = rowMeans((g1+g2)[rownames(out), stages %in% c(post,pre)]) >= 1
out2 = out[use & (out$Cat == 'Up'),]
out2$FCmax = rowMeans(log(gMax[rownames(out2),stages == post] + 0.5,2), na.rm=T) -
rowMeans(log(gMax[rownames(out2),stages == pre] + 0.5,2), na.rm=T)
out2$FCmin =
rowMeans(log(gMin[rownames(out2),stages == post] + 0.5,2), na.rm=T) -
rowMeans(log(gMin[rownames(out2),stages == pre] + 0.5,2), na.rm=T)
out2 = data.frame(FC = c(out2$FCmax, out2$FCmin), Cat = rep(out2$Cat,2), Type = rep(c('Max','Min'), each = nrow(out2)))

Fig4c = ggplot(out2, aes(x=FC, fill=Type)) + geom_histogram(position='identity', alpha=.75, binwidth=2/15) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab('Genes') + xlab('log2 fold change') + scale_x_continuous(lim = c(-4,6), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

ggarrange(Fig4a, Fig4b, Fig4c, ncol = 3)
}

svg("Rice_repFig4.svg", width=7, height=2.5)
repFig4(Os_UMIs, Os_stages, "UM Bi", "UM Mono", Os_g1, Os_g2)
dev.off()

svg("Tomato_repFig4.svg", width=7, height=2.5)
repFig4(Sl_UMIs, Sl_stages, "UM Bi", "UM Mono", Sl_g1, Sl_g2)
dev.off()

svg("Tomato_repFig4_UMbi_vs_UMBM.svg", width=7, height=2.5)
repFig4(Sl_UMIs, Sl_stages, "UM Bi", "UM/BM", Sl_g1, Sl_g2)
dev.off()



###No cell controls
#mean / median / MAD / quantile(, p = c(.9,.95,.98)) for NCC, for Arabidopsis, whats mean / median / MAD / top 3
load("At_data_2_2026.RData")
At_NCC = colSums(D[,rownames(At_meta[At_meta$No_cell_well == "Y",])])
length(At_NCC)
#[1] 20
summary(At_NCC)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1028    3747    5535    8112   10516   22950
mad(At_NCC)
#[1] 4783.609
At_NCC[order(At_NCC, decreasing=TRUE)][1:3]
#A268-287_80s A288-301_64s A288-301_25s 
#       22950        22692        15699

load("TomatoRemap_11_2024-7_2025_with_meta.RData")
Sl_NCC = colSums(D2[,rownames(Sl_meta[Sl_meta$No_cell2 == "Y",])])
length(Sl_NCC)                     
#[1] 126
summary(Sl_NCC)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    232    1734    2988    3982    4744   21091
mad(Sl_NCC)
#[1] 2104.551
quantile(Sl_NCC, p = c(.9,.95,.98))
#  90%   95%   98% 
# 7239 12073 14346


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
Os_NCC = colSums(D2[,rownames(Rice_meta[Rice_meta$No_cell2 == "Y",])])
length(Os_NCC)
[1] 186
summary(Os_NCC)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    711    3068    4692    5891    7152   27917
mad(Os_NCC)
#[1] 2900.707
quantile(Os_NCC, p = c(.9,.95,.98))
#    90%     95%     98% 
#11367.5 14248.0 18339.1



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
Zm_NCC = colSums(D[,rownames(Meta[Meta$No_cell2 == "Y",])])

length(Zm_NCC)
#[1] 98
Zm_NCC_BBB = Zm_NCC[-grep("ETOH", names(Zm_NCC))]
length(Zm_NCC_BBB)
#[1] 84
Zm_NCC_ETOH = Zm_NCC[grep("ETOH", names(Zm_NCC))]
length(Zm_NCC_ETOH)
#[1] 14
summary(Zm_NCC)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    316    1005    1412    3997    2388   50931
mad(Zm_NCC)
#[1] 859.1667
quantile(Zm_NCC, p = c(.9,.95,.98))
#     90%      95%      98% 
#11190.10 16894.80 22890.16
summary(Zm_NCC_BBB)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  316.0   940.8  1271.0  1695.6  1798.5 10861.0
mad(Zm_NCC_BBB)    
#[1] 575.2488
quantile(Zm_NCC_BBB, p = c(.9,.95,.98))
#    90%     95%     98% 
#2921.70 4847.90 6445.74
summary(Zm_NCC_ETOH)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   7400   10393   14051   17802   20451   50931
mad(Zm_NCC_ETOH)   
#[1] 7834.058 
Zm_NCC_ETOH[order(Zm_NCC_ETOH, decreasing=TRUE)][1:3]
#M37-48_ETOH_51s M37-48_ETOH_32s M37-48_ETOH_69s 
#          50931           34016           22180 



#raw UMIcount plots
load("Arabidopsis.RData")
load("Rice.RData")
load("Maize_merge.RData")
load("Maize.RData")
load("Tomato.RData")

hist(Zm_FracMono) #threshold of 40
hist(Zm2_FracMono) #threshold of 40
hist(Os_FracMono) #threshold of 40
hist(Sl_FracMono) #threshold of 40
hist(At_HapFrac) #threshold of 40

Zm_merge_stages = c(as.character(Zm_stages), as.character(Zm2_stages))
names(Zm_merge_stages) = c(names(Zm_stages), names(Zm2_stages))
Zm_merge_stages[rownames(Zm_Meta[grep("mix", Zm_Meta[,5]),])] = "UM/BM"
Zm_merge_stages[which(Zm_merge_stages == "Early" | Zm_merge_stages == "Middle" | Zm_merge_stages == "Late")] = "UM"
#Zm_merge_stages[which(Zm_merge_stages == "TP_shed" | Zm_merge_stages == "TP")] = "Pollen"
Zm_merge_stages = factor(Zm_merge_stages, levels = c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "TP", "TP_shed"))


UMI_summary = data.frame(UMIs = c(colSums(Zm_merge_UMIs), colSums(Os_UMIs), colSums(Sl_UMIs), colSums(At_UMIs))/10^3,
    stages = c(Zm_merge_stages, Os_stages, Sl_stages, At_stages), 
    species = c(rep("Zm", length(Zm_merge_stages)), rep("Os", length(Os_stages)), 
    rep("Sl", length(Sl_stages)), rep("At", length(At_stages))))
UMI_summary$stages = factor(UMI_summary$stages, levels=c("Tetrad", "UM","UM Bi","UM Mono", "UM/BM", "BM", "BM/TP", "TP", "BM_shed", "TP_shed"))


svg('UMI_plot.svg', width=4, height=3)
ggplot(UMI_summary, aes(x = UMIs, y = species, color = stages)) +
  geom_jitter(size = 0.5, height=0.25) +
  theme_minimal() + #geom_jitter() + 
  labs(
    x = paste("Transcripts per pollen precursor (x1000)
    with >", cutoff, "transcripts", sep = " "),
    y = "Species") +
    scale_color_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#FD6F30", "UM Bi" = "#FD6F30", "UM Mono" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "BM_shed" = "#BBC9E5", "TP" = "#A4A4D5", "TP_shed" = "#BBC9E5"))
dev.off()



UMIplot = function(Sp, cutoff){
ggplot(UMI_summary[UMI_summary$species == Sp & UMI_summary$UMIs > cutoff,], aes(x = stages, y = UMIs, color = stages)) +
  geom_jitter(size = 0.5, width=0.35) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = Sp,
    x = "",
    y = "") + theme(legend.position="none") +
    scale_color_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#FD6F30", "UM Bi" = "#FD6F30", "UM Mono" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "BM_shed" = "#BBC9E5", "TP" = "#A4A4D5", "TP_shed" = "#BBC9E5"))
}

UMIplotl = function(Sp, cutoff){
ggplot(UMI_summary[UMI_summary$species == Sp & UMI_summary$UMIs > cutoff,], aes(x = stages, y = UMIs, color = stages)) +
  geom_jitter(size = 0.5, width=0.35) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = Sp,
    x = "",
    y = "Transcripts per pollen precursor 
    (x1000)") + theme(legend.position="none") +
    scale_color_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#FD6F30", "UM Bi" = "#FD6F30", "UM Mono" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "BM_shed" = "#BBC9E5", "TP" = "#A4A4D5", "TP_shed" = "#BBC9E5"))
}


set.seed(2)
svg('UMI_plot2.svg', width=8, height=3)
ggarrange(UMIplotl("Zm", 5), UMIplot("Os", 10), UMIplot("Sl", 10), UMIplot("At", 23), ncol = 4, nrow = 1, align = "h")
dev.off()

set.seed(2)
svg('UMI_plot3.svg', width=8, height=3)
ggarrange(UMIplotl("Zm", 5), UMIplot("Os", 15), UMIplot("Sl", 10), UMIplot("At", 23), ncol = 4, nrow = 1, align = "h")
dev.off()

UMIplot2 = function(Sp, cutoff, St){
ggplot(UMI_summary[UMI_summary$species == Sp & UMI_summary$UMIs > cutoff & UMI_summary$stages == St,], aes(x = UMIs, fill = stages)) +
  geom_histogram(binwidth=5) +
  theme_minimal() + coord_cartesian(xlim = c(0, max(UMI_summary[UMI_summary$species == Sp & UMI_summary$UMIs > cutoff,1]) +10)) +
  labs(
    title = paste(Sp, St, sep = " "),
    x = "",
    y = "") + theme_classic() + theme(legend.position="none") +
    scale_fill_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#FD6F30", "UM Bi" = "#FD6F30", "UM Mono" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "BM_shed" = "#BBC9E5", "TP" = "#A4A4D5", "TP_shed" = "#BBC9E5"))
}


UMIplot2l = function(Sp, cutoff, St){
ggplot(UMI_summary[UMI_summary$species == Sp & UMI_summary$UMIs > cutoff & UMI_summary$stages == St,], aes(x = UMIs, fill = stages)) +
  geom_histogram(binwidth=5) +
  theme_minimal() + coord_cartesian(xlim = c(0, max(UMI_summary[UMI_summary$species == Sp & UMI_summary$UMIs > cutoff,1]) +10)) +
  labs(
    title = paste(Sp, St, sep = " "),
    x = "Transcripts per pollen
    precursor (x1000)",
    y = "# of samples") + theme_classic() + theme(legend.position="none") +
    scale_fill_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#FD6F30", "UM Bi" = "#FD6F30", "UM Mono" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "BM_shed" = "#BBC9E5", "TP" = "#A4A4D5", "TP_shed" = "#BBC9E5"))
}

Zm_UMI_hist = ggarrange(UMIplot2l("Zm", 5, "Tetrad"), UMIplot2("Zm", 5, "UM"), UMIplot2("Zm", 5, "UM/BM"), UMIplot2("Zm", 5, "BM"),
    UMIplot2("Zm", 5, "BM/TP"), UMIplot2("Zm", 5, "TP"), UMIplot2("Zm", 5, "TP_shed"), ncol = 7, nrow = 1, align = "h")

Os_UMI_hist = ggarrange(UMIplot2l("Os", 10, "Tetrad"), UMIplot2("Os", 10, "UM Bi"), UMIplot2("Os", 10, "UM Mono"), UMIplot2("Os", 10, "UM/BM"),
    UMIplot2("Os", 10, "BM"), UMIplot2("Os", 10, "BM/TP"), UMIplot2("Os", 10, "TP"), UMIplot2("Os", 10, "TP_shed"), ncol = 8, nrow = 1, align = "h")

Sl_UMI_hist = ggarrange(UMIplot2l("Sl", 10, "Tetrad"), UMIplot2("Sl", 10, "UM Bi"), UMIplot2("Sl", 10, "UM Mono"), UMIplot2("Sl", 10, "UM/BM"),
    UMIplot2("Sl", 10, "BM"), UMIplot2("Sl", 10, "BM_shed"), ncol = 6, nrow = 1, align = "h")

At_UMI_hist = ggarrange(UMIplot2l("At", 23, "Tetrad"), UMIplot2("At", 23, "UM"), UMIplot2("At", 23, "UM/BM"),
    UMIplot2("At", 23, "BM"), UMIplot2("At", 23, "BM/TP"), UMIplot2("At", 23, "TP"), ncol = 6, nrow = 1, align = "h")




svg('UMI_plot4.svg', width=14, height=8)
ggarrange(Zm_UMI_hist,Os_UMI_hist,Sl_UMI_hist,At_UMI_hist, ncol = 1, nrow = 4, align = "v")
dev.off()




repFig4(Zm_UMIs, Zm_stages, "Late", "BM", Zm_g1, Zm_g2)
repFig4(At_UMIs, At_stages, "BM", "TP", At_g1, At_g2)



> table(Sl_stages)                                                       
Sl_stages
 Tetrad   UM Bi UM Mono   UM/BM      BM BM_shed 
     19       5       4      10      47      31
mean(colSums(Sl_UMIs[,names(Sl_stages[which(Sl_Meta[,8] == "tetrad_sep")])])) #only 2 sep tetrads
[1] 14409.5
> mean(colSums(Sl_UMIs[, names(Sl_stages[which(Sl_stages == "UM Bi")])]))
[1] 13846.6
> mean(colSums(Sl_UMIs[, names(Sl_stages[which(Sl_stages == "UM Mono")])]))
[1] 27165
> mean(colSums(Sl_UMIs[, names(Sl_stages[which(Sl_stages == "UM/BM")])]))  
[1] 43114.6
> mean(colSums(Sl_UMIs[, names(Sl_stages[which(Sl_stages == "BM")])]))   
[1] 111594.2
> mean(colSums(Sl_UMIs[, names(Sl_stages[which(Sl_stages == "BM_shed")])]))
[1] 111514.9


table(At_stages)
At_stages
Tetrad     UM  UM/BM     BM  BM/TP     TP 
     2     38     11     58      3     79
mean(colSums(At_UMIs[, names(At_stages[which(At_stages == "Tetrad")])])) 
[1] 32111.5
> mean(colSums(At_UMIs[, names(At_stages[which(At_stages == "UM")])]))    
[1] 37525.08
> mean(colSums(At_UMIs[, names(At_stages[which(At_stages == "UM/BM")])]))
[1] 41752.91
> mean(colSums(At_UMIs[, names(At_stages[which(At_stages == "BM")])]))   
[1] 40046.14
> mean(colSums(At_UMIs[, names(At_stages[which(At_stages == "BM/TP")])]))
[1] 29262.67
> mean(colSums(At_UMIs[, names(At_stages[which(At_stages == "TP")])]))   
[1] 65197.92



boxplot(Sl_UMIs["Solyc10g007270",] ~ Sl_stages, pch=19, xlab = "Stage", ylab = "Raw transcript count", main = "Tomato LAT52 expression")
boxplot(Sl_TPM["Solyc10g007270",] ~ Sl_stages, pch=19, xlab = "Stage", ylab = "Transcripts per million", main = "Tomato LAT52 expression")

Sl_LAT52_frac = Sl_g1["Solyc10g007270",]/(Sl_g1["Solyc10g007270",]+Sl_g2["Solyc10g007270",])
hist(Sl_LAT52_frac*100, breaks=10, xlab="% transcripts from lycopersicum allele", ylab="# of pollen precursors or pollen", main="Tomato LAT52 Allelic Bias")



LAT52_UMIs_stages = data.frame(UMIs = Sl_UMIs["Solyc10g007270",], stages = Sl_stages)

MeanLAT52plot = function(Sp, cutoff){
ggplot(LAT52_UMIs_stages, aes(x = stages, y = UMIs, fill = stages)) +
  geom_bar(position="fill") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Tomato LAT52 Expression",
    x = "Stage",
    y = "Transcripts per pollen precursor") + theme(legend.position="none") +
    scale_color_manual(values = c("Tetrad" = "#EB1E2C", "UM" = "#FD6F30", "UM Bi" = "#FD6F30", "UM Mono" = "#F9A729",
     "UM/BM" = "#F9D23C", "BM" = "#5FBB68", "BM/TP" = "#64CDCC", "BM_shed" = "#BBC9E5", "TP" = "#A4A4D5", "TP_shed" = "#BBC9E5"))
}




bootFn = function(x, trim = 0.2, nboot = 2000) {
	out = sapply(1:nboot, function(i) { mean(sample(x, replace = T), trim = trim, na.rm=T) })
	return(c(mean(out), mean(out) - sd(out), mean(out) + sd(out)))
}


LAT52ByStage = data.frame(t(matrix(unlist(by(Sl_UMIs["Solyc10g007270",], Sl_stages, bootFn)), nrow = 3)))
colnames(LAT52ByStage) = c('val', 'lower', 'upper')
LAT52ByStage$Stage = unique(Sl_stages[order(Sl_stages)])

svg('LAT52_MeanExp_by_stage.svg', width=4, height=2.5)
ggplot(LAT52ByStage) + geom_bar(aes(x=Stage, y=val), stat='identity', fill='gray') + geom_errorbar(aes(x=Stage, ymin=lower, ymax=upper), width=0.3, size=.9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90)) + scale_y_continuous(expand=c(0,0)) + labs(x= 'Stage', y = 'Transcripts per\npollen precursor', title ='Tomato LAT52 Expression') + theme(plot.title = element_text(hjust = 0.5))
dev.off()



Sl_ratio=Sl_g1/(Sl_g1+Sl_g2)
Sl_geneFrac = rowMeans(abs(Sl_ratio - 0.5)+0.5, na.rm=T)
Sl_geneFrac[rowSums(Sl_g1+Sl_g2 >10)<3] = NA
length(which(Sl_geneFrac >= 0.99))                    
[1] 593
length(which(Sl_geneFrac >= 0.95))                      
[1] 1290
length(which(rowSums(Sl_g1+Sl_g2 >10)>=3))
[1] 1477    #40.1% >=0.99, 87.3% >=0.95
hist(Sl_geneFrac)
SlHapGs = rownames(Sl_UMIs[which(Sl_geneFrac >= 0.99),])
SlHapGsMeans = by(t(Sl_TPM[SlHapGs,]), Sl_stages, colMeans)
calcMeanExp_HapGs = function(x){
    rowMeans(Sl_TPM[SlHapGs,which(Sl_stages==x)], na.rm=T) #mean TPM of each gene for each group
}
SlHapGsByStage = sapply(levels(Sl_stages), calcMeanExp_HapGs)
SlHapGs2 = names(which(apply(SlHapGsByStage, 1, max) >= 100)) #377 genes
SlHapGs3 = names(which(rowMeans(Sl_ratio[SlHapGs2,], na.rm=T) >0.2 & rowMeans(Sl_ratio[SlHapGs2,], na.rm=T) <0.8)) #306 genes

HapGplot = function(Exp, UMIs, gene, stages, keep){
Exp = data.frame(t(matrix(unlist(by(UMIs[gene,], stages, bootFn)), nrow = 3)))
colnames(Exp) = c('val', 'lower', 'upper')
Exp$Stage = unique(stages[order(stages)])
Exp = Exp[Exp$Stage %in% keep,]
ggplot(Exp) + geom_bar(aes(x=Stage, y=val), stat='identity', fill='gray') + geom_errorbar(aes(x=Stage, ymin=lower, ymax=upper), width=0.3, size=.9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90)) + scale_y_continuous(expand=c(0,0)) + labs(x= 'Stage', y = 'Transcripts per\npollen precursor', title = gene) + theme(plot.title = element_text(hjust = 0.5))
}

Sl_keep = levels(Sl_stages)
Sl_HapGsPlot = function(x){
    HapGplot(SlHapGsByStage, Sl_UMIs, SlHapGs3[x], Sl_stages, Sl_keep)
}

ggarrange(Sl_HapGsPlot(1), Sl_HapGsPlot(2), Sl_HapGsPlot(3), Sl_HapGsPlot(4), Sl_HapGsPlot(5), Sl_HapGsPlot(6), 
    Sl_HapGsPlot(7), Sl_HapGsPlot(8), Sl_HapGsPlot(9), Sl_HapGsPlot(10), Sl_HapGsPlot(11), Sl_HapGsPlot(12),
    Sl_HapGsPlot(13), Sl_HapGsPlot(14), Sl_HapGsPlot(15), Sl_HapGsPlot(16), Sl_HapGsPlot(17), Sl_HapGsPlot(18), 
    Sl_HapGsPlot(19), Sl_HapGsPlot(20), nrow=4, ncol=5)



Os_ratio=Os_g1/(Os_g1+Os_g2)
Os_geneFrac = rowMeans(abs(Os_ratio - 0.5)+0.5, na.rm=T)
Os_geneFrac[rowSums(Os_g1+Os_g2 >10)<3] = NA
length(which(Os_geneFrac >= 0.99))                    
[1] 56
length(which(Os_geneFrac >= 0.95))                      
[1] 281
length(which(rowSums(Os_g1+Os_g2 >10)>=3))
[1] 529     #10.6% >=0.99, 53.1% >=0.95
hist(Os_geneFrac)
OsHapGs = rownames(Os_UMIs[which(Os_geneFrac >= 0.99),])
OsHapGsMeans = by(t(Os_TPM[OsHapGs,]), Os_stages, colMeans)
calcMeanExp_HapGs = function(x){
    rowMeans(Os_TPM[OsHapGs,which(Os_stages==x)], na.rm=T) #mean TPM of each gene for each group
}
OsHapGsByStage = sapply(levels(Os_stages)[c(1:5,7)], calcMeanExp_HapGs)
OsHapGs2 = names(which(apply(OsHapGsByStage, 1, max) >= 100))
OsHapGs3 = names(which(rowMeans(Os_ratio[OsHapGs2,], na.rm=T) >0.2 & rowMeans(Os_ratio[OsHapGs2,], na.rm=T) <0.8))

HapGplot = function(Exp, UMIs, gene, stages, keep){
Exp = data.frame(t(matrix(unlist(by(UMIs[gene,], stages, bootFn)), nrow = 3)))
colnames(Exp) = c('val', 'lower', 'upper')
Exp$Stage = unique(stages[order(stages)])
Exp = Exp[Exp$Stage %in% keep,]
ggplot(Exp) + geom_bar(aes(x=Stage, y=val), stat='identity', fill='gray') + geom_errorbar(aes(x=Stage, ymin=lower, ymax=upper), width=0.3, size=.9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90)) + scale_y_continuous(expand=c(0,0)) + labs(x= 'Stage', y = 'Transcripts per\npollen precursor', title = gene) + theme(plot.title = element_text(hjust = 0.5))
}

Os_keep = c("Tetrad", "UM Bi", "UM Mono", "UM/BM", "BM", "TP")
Os_HapGsPlot = function(x){
    HapGplot(OsHapGsByStage, Os_UMIs, OsHapGs3[x], Os_stages, Os_keep)
}
#ggarrange(sapply(1:19, Os_HapGsPlot), nrow=4, ncol=5)

ggarrange(Os_HapGsPlot(1), Os_HapGsPlot(2), Os_HapGsPlot(3), Os_HapGsPlot(4), Os_HapGsPlot(5), Os_HapGsPlot(6), 
    Os_HapGsPlot(7), Os_HapGsPlot(8), Os_HapGsPlot(9), Os_HapGsPlot(10), Os_HapGsPlot(11), Os_HapGsPlot(12),
    Os_HapGsPlot(13), Os_HapGsPlot(14), Os_HapGsPlot(15), Os_HapGsPlot(16), Os_HapGsPlot(17), Os_HapGsPlot(18), 
    Os_HapGsPlot(19), nrow=4, ncol=5)



Zm_merge_ratio=Zm_merge_g1/(Zm_merge_g1+Zm_merge_g2)
Zm_merge_geneFrac = rowMeans(abs(Zm_merge_ratio - 0.5)+0.5, na.rm=T)
Zm_merge_geneFrac[rowSums(Zm_merge_g1+Zm_merge_g2 >10)<3] = NA
length(which(Zm_merge_geneFrac >= 0.99))  
[1] 238
length(which(Zm_merge_geneFrac >= 0.95))   
[1] 1185
length(which(rowSums(Zm_merge_g1+Zm_merge_g2 >10)>=3))
[1] 3194    #7.45% >=0.99, 37.1% >=0.95
hist(Zm_merge_geneFrac)

At_ratio=At_g1/(At_g1+At_g2)
At_geneFrac = rowMeans(abs(At_ratio - 0.5)+0.5, na.rm=T)
At_geneFrac[rowSums(At_g1+At_g2 >10)<3] = NA
length(which(At_geneFrac >= 0.99))                    
[1] 31
length(which(At_geneFrac >= 0.95))                    
[1] 210
length(which(rowSums(At_g1+At_g2 >10)>=3))
[1] 1429    #2.17% >=0.99, 14.7% >=0.95
hist(At_geneFrac)

df <- data.frame(
  geneFrac = c(Zm_merge_geneFrac, Os_geneFrac, 
            Sl_geneFrac, At_geneFrac),
  species = c(rep("Maize", times=length(Zm_merge_geneFrac)), rep("Rice", times=length(Os_geneFrac)), 
  rep("Tomato", times=length(Sl_geneFrac)), rep("Arabidopsis", times=length(At_geneFrac)))
)
df$species = factor(df$species, levels=c("Maize", "Tomato", "Arabidopsis", "Rice"))

ggplot(df, aes(x = geneFrac, fill = species)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, color = "white") +
  scale_fill_manual(values = c("#FFD93D", "#FF6B6B","#6BCB77", "#4D96FF")) + 
  labs(title = "", x = "Fraction of Monoallelic Transcripts", y = "# of Genes") +
  theme_classic()


At_stages2 = as.character(At_stages)
names(At_stages2) = At_stages
At_stages2[which(At_HapFrac >=40)] = paste(At_stages2[which(At_HapFrac >=40)], "mono", sep="_")
At_stages2[which(At_HapFrac <40)] = paste(At_stages2[which(At_HapFrac <40)], "bi", sep="_")
At_stages2 = factor(At_stages2, levels=c("Tetrad_bi", "UM_bi", "UM_mono", "UM/BM_bi", "UM/BM_mono", "BM_bi", "BM_mono", "BM/TP_bi", "BM/TP_mono", "TP_bi", "TP_mono"))

AtHapGs = rownames(At_UMIs[which(At_geneFrac >= 0.99),])
AtHapGsMeans = by(t(At_TPM[AtHapGs,]), At_stages2, colMeans)
calcMeanExp_HapGs = function(x){
    rowMeans(At_TPM[AtHapGs,which(At_stages2==x)], na.rm=T) #mean TPM of each gene for each group
}
AtHapGsByStage = sapply(c("Tetrad_bi", "UM_bi", "UM_mono", "BM_bi", "BM_mono", "TP_bi", "TP_mono"), calcMeanExp_HapGs)
AtHapGs2 = names(which(apply(AtHapGsByStage, 1, max) >= 100))
AtHapGs3 = names(which(rowMeans(At_ratio[AtHapGs2,], na.rm=T) >0.2 & rowMeans(At_ratio[AtHapGs2,], na.rm=T) <0.8))


HapGplot = function(x){
AtHapGsByStage = data.frame(t(matrix(unlist(by(At_UMIs[x,], At_stages2, bootFn)), nrow = 3)))
colnames(AtHapGsByStage) = c('val', 'lower', 'upper')
AtHapGsByStage$Stage = unique(At_stages2[order(At_stages2)])
AtHapGsByStage = AtHapGsByStage[AtHapGsByStage$Stage %in% c("Tetrad_bi", "UM_bi", "UM_mono", "BM_bi", "BM_mono", "TP_bi", "TP_mono"),]
ggplot(AtHapGsByStage) + geom_bar(aes(x=Stage, y=val), stat='identity', fill='gray') + geom_errorbar(aes(x=Stage, ymin=lower, ymax=upper), width=0.3, size=.9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90)) + scale_y_continuous(expand=c(0,0)) + labs(x= 'Stage', y = 'Transcripts per\npollen precursor', title =x) + theme(plot.title = element_text(hjust = 0.5))
}

ggarrange(HapGplot(AtHapGs3[1]), HapGplot(AtHapGs3[2]), HapGplot(AtHapGs3[3]), HapGplot(AtHapGs3[4]), HapGplot(AtHapGs3[5]), HapGplot(AtHapGs3[6]), HapGplot(AtHapGs3[7]), nrow=2, ncol=4)



AtUMIsByStage = data.frame(t(matrix(unlist(by(colSums(At_UMIs), At_stages2, bootFn)), nrow = 3)))
colnames(AtUMIsByStage) = c('val', 'lower', 'upper')
AtUMIsByStage$Stage = unique(At_stages2[order(At_stages2)])
AtUMIsByStage = AtUMIsByStage[AtUMIsByStage$Stage %in% c("Tetrad_bi", "UM_bi", "UM_mono", "BM_bi", "BM_mono", "TP_bi", "TP_mono"),]
ggplot(AtUMIsByStage) + geom_bar(aes(x=Stage, y=val), stat='identity', fill='gray') + geom_errorbar(aes(x=Stage, ymin=lower, ymax=upper), width=0.3, size=.9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90)) + scale_y_continuous(expand=c(0,0)) + labs(x= 'Stage', y = 'Transcripts per\npollen precursor', title ='Arabidopsis') + theme(plot.title = element_text(hjust = 0.5))



SlUMIsByStage = data.frame(t(matrix(unlist(by(colSums(Sl_UMIs), Sl_stages, bootFn)), nrow = 3)))
colnames(SlUMIsByStage) = c('val', 'lower', 'upper')
SlUMIsByStage$Stage = unique(Sl_stages[order(Sl_stages)])
#AtUMIsByStage = AtUMIsByStage[AtUMIsByStage$Stage %in% c("Tetrad_bi", "UM_bi", "UM_mono", "BM_bi", "BM_mono", "TP_bi", "TP_mono"),]
ggplot(SlUMIsByStage) + geom_bar(aes(x=Stage, y=val), stat='identity', fill='gray') + geom_errorbar(aes(x=Stage, ymin=lower, ymax=upper), width=0.3, size=.9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90)) + scale_y_continuous(expand=c(0,0)) + labs(x= 'Stage', y = 'Transcripts per\npollen precursor', title ='Tomato') + theme(plot.title = element_text(hjust = 0.5))
