#Zm_Os_ortholog_comparison3
#Maize vs. Rice ortholog comparison
setwd('C:/Users/julia/OneDrive/Desktop/Grad School/Nelms lab/Bioinformatics/R')

HomTab4 <- read.delim("Zm_Os_At_Sl_mart_export.txt")
colnames(HomTab4) = c("Zm_geneID", "Zm_gene_name", "Os_geneID", "Os_gene_name", "Os_homology_type", "Os_percent_ID", "Os_confidence", 
"At_geneID", "At_gene_name", "At_homology_type", "At_percent_ID", "At_confidence",
"Sl_geneID", "Sl_gene_name", "Sl_homology_type", "Sl_percent_ID", "Sl_confidence")
ConfHomZmOs = HomTab4[which(HomTab4$Os_confidence == 1),]
length(unique(ConfHomZmOs$Os_geneID))
#[1] 16925

One2oneZmOs = ConfHomZmOs[which(ConfHomZmOs$Os_homology_type == "ortholog_one2one"),]
One2oneZmOs = One2oneZmOs[which(!duplicated(One2oneZmOs$Zm_geneID)),1:7]
rownames(One2oneZmOs) = One2oneZmOs[,1]
dim(One2oneZmOs)
#[1] 9541    7 #9541 1-1 orthologs rated as confident by Ensembl

ZmOs = read.table('Zm_Os_dNdS.tab', skip = 2, sep = '\t', as.is = T)
ZmOs = cbind(ZmOs, sub("[.].*", "", sub(".*LOC", "LOC", ZmOs[,4])), sub("_.*", "", sub("[|].*", '', sub(".*CDS:", "", ZmOs[,8]))), as.numeric(sub('.*[|]','',ZmOs[,4])))
colnames(ZmOs)[13:15] = c('Os_G','Zm_G','Identity')

ZmOs_NS = ZmOs[ZmOs$Zm_G %in% rownames(One2oneZmOs),]
ZmOs_NS = ZmOs_NS[which(!duplicated(ZmOs_NS[,14])),]
ZmOs_NS = ZmOs_NS[which(!duplicated(ZmOs_NS[,13])),]

rownames(ZmOs_NS) = ZmOs_NS[,14]
ZmOs_NS$dNdS = as.numeric(ZmOs_NS[,2]) / as.numeric(ZmOs_NS[,1])
length(ZmOs_NS$dNdS)
#[1] 6889 #6889 maize genes that are confident 1-1 orthologs w/ rice and have dN/dS data from CoGe


load("Rice.RData")
load("Maize_merge.RData")
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
#Zm_merge_stages[which(Zm_merge_stages == "TP_shed" | Zm_merge_stages == "TP")] = "Pollen"
Zm_merge_stages = factor(Zm_merge_stages, levels = c("Tetrad", "UM", "UM/BM", "BM", "BM/TP", "TP", "TP_shed"))

ZmOs_NS = ZmOs_NS[ZmOs_NS$Zm_G %in% rownames(Zm_merge_TPM),]

Zm_StageMean = function(stage){
    rowMeans(Zm_merge_TPM[ZmOs_NS$Zm_G,names(Zm_merge_stages[which(Zm_merge_stages == stage)])])
}
Os_StageMean = function(stage){
    rowMeans(Os_TPM[ZmOs_NS$Os_G,names(Os_stages[which(Os_stages == stage)])])
}

ZmOs_comp = as.data.frame(cbind(sapply(levels(Zm_merge_stages), Zm_StageMean), sapply(levels(Os_stages), Os_StageMean)))
colnames(ZmOs_comp) = c(paste0("Zm_", levels(Zm_merge_stages)), paste0("Os_", levels(Os_stages)))
#ZmOs_comp$dNdS = ZmOs_NS$dNdS

ZmOs_comp$maxBi_Zm = apply(ZmOs_comp[,1:2], 1, max)
ZmOs_comp$maxMono_Zm = apply(ZmOs_comp[,4:7], 1, max) #skipping UM/BM in maize for now as it is a mix of bi and mono samps
ZmOs_comp$maxBi_Os = apply(ZmOs_comp[,8:9], 1, max)
ZmOs_comp$maxMono_Os = apply(ZmOs_comp[,10:15], 1, max)



ZmSb = read.table('Zm_v4_sorghum_dNdS.tab', skip = 2, sep = '\t', as.is = T)
ZmSb = cbind(ZmSb, sub("_.*", "", sub(".*CDS:", "", ZmSb[,4])), sub("_.*", "", sub("[|].*", '', sub(".*CDS:", "", ZmSb[,8]))), as.numeric(sub('.*[|]','',ZmSb[,4])))
colnames(ZmSb)[13:15] = c('Zm_G','Sb_G','Identity')

## Remove putative orthologs that are not reciprocal best hits ##
removes = rep(FALSE, nrow(ZmSb))
for (sgene in as.character(unique(ZmSb[duplicated(ZmSb[,13]),13]))) {
	tr = which(ZmSb[,13] == sgene)
	removes[tr[rank(-ZmSb[tr,15]) != 1]] = TRUE
}
for (mgene in as.character(unique(ZmSb[duplicated(ZmSb[,14]),14]))) {
	tr = which(ZmSb[,14] == mgene)
	removes[tr[rank(-ZmSb[tr,15]) != 1]] = TRUE
}
ZmSb = ZmSb[!removes,]
##

ZmSb = ZmSb[ZmSb[,15] >= 80,]  # Require at least 80% identity between orthologs when calculating dNdS
ZmSb = ZmSb[-which(is.na(as.numeric(ZmSb[,1]))),] # Remove genes that lack calculated synonymous and non-synonymous rates or are undefined
ZmSb = ZmSb[-which(as.numeric(ZmSb[,1]) == 0),] # Remove genes with 0 synonymous subs 
ZmNS = data.frame(ZmSb_dN = as.numeric(ZmSb[,2]), ZmSb_dS = as.numeric(ZmSb[,1]))
rownames(ZmNS) = ZmSb[,13]


#Downloaded table with all gene model versions: https://download.maizegdb.org/Pan-genes/Pan-Zea/
PanGeneTab = read.delim("pan-zea.v4.pan-genes_table.tsv")
V5toV4tab = PanGeneTab[,6:7]
V5toV4tab = V5toV4tab[which(rowSums(V5toV4tab == "NONE") == 0),]
V5toV4tab$Zm00001d = sub("_.*", "", V5toV4tab$Zm00001d)
V5toV4tab$Zm00001eb = sub("_.*", "", V5toV4tab$Zm00001eb)

vConvert = V5toV4tab[which(V5toV4tab$Zm00001d %in% rownames(ZmNS)),]
rownames(vConvert) = vConvert[,1]

dim(ZmNS)
#[1] 14892      2
ZmNSv5 = ZmNS[rownames(ZmNS) %in% vConvert[,1],]
rownames(ZmNSv5) = vConvert[rownames(ZmNSv5),2]
dim(ZmNSv5)
#[1] 13194     2
ZmNSv5 = ZmNSv5[-which(is.na(ZmNSv5[,1])),]

#### Making 30TPMcat and 100TPMcat
makeCat = function(TPM, cat, ZmOs_comp){
  ZmOs_comp[[cat]] = ifelse(
    ZmOs_comp$maxMono_Os > TPM & ZmOs_comp$maxMono_Zm <= TPM, "OsOnly",
    ifelse(
      ZmOs_comp$maxMono_Zm > TPM & ZmOs_comp$maxMono_Os <= TPM, "ZmOnly",
      ifelse(
        ZmOs_comp$maxMono_Os > TPM & ZmOs_comp$maxMono_Zm > TPM, "BothOsZm",
        "NeitherOsZm"
      )
    )
  )
  return(ZmOs_comp)
}
ZmOs_comp = makeCat(TPM = 30, cat = "30TPMcat", ZmOs_comp = ZmOs_comp)
ZmOs_comp = makeCat(TPM = 100, cat = "100TPMcat", ZmOs_comp = ZmOs_comp)
####


ZmOs_comp$Os_G = ZmOs_NS$Os_G
ZmOs_comp$Zm_G = ZmOs_NS$Zm_G



ZmNSv5 = ZmNSv5[rownames(ZmNSv5) %in% rownames(ZmOs_comp),]

ZmOs_comp$ZmSb_dN = NA
ZmOs_comp[rownames(ZmNSv5),"ZmSb_dN"] = ZmNSv5$ZmSb_dN
ZmOs_comp$ZmSb_dS = NA
ZmOs_comp[rownames(ZmNSv5),"ZmSb_dS"] = ZmNSv5$ZmSb_dS
ZmOs_comp$"30TPMcat" = factor(ZmOs_comp$"30TPMcat", levels=c("NeitherOsZm","OsOnly","ZmOnly","BothOsZm"))
ZmOs_comp$"100TPMcat" = factor(ZmOs_comp$"100TPMcat", levels=c("NeitherOsZm","OsOnly","ZmOnly","BothOsZm"))






library(dplyr)
library(ggplot2)

ggplot(
  ZmOs_comp,
  aes(x = `100TPMcat`, y = ZmSb_dN / ZmSb_dS)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.25) +
  scale_y_log10() +
  theme_classic() +
  labs(
    x = "Haploid expression category",
    y = "dN/dS (log scale)"
  )

pairwise.wilcox.test(ZmOs_comp$ZmSb_dN / ZmOs_comp$ZmSb_dS, ZmOs_comp$"30TPMcat", p.adjust = 'holm')  # p-values
pairwise.wilcox.test(ZmOs_comp$ZmSb_dN / ZmOs_comp$ZmSb_dS, ZmOs_comp$"100TPMcat", p.adjust = 'holm')  # p-values




ZmOsNS = data.frame(ZmOs_dN = as.numeric(ZmOs_NS[,2]), ZmOs_dS = as.numeric(ZmOs_NS[,1]))
rownames(ZmOsNS) = ZmOs_NS[,14]

ZmOs_comp$ZmOs_dN = NA
ZmOs_comp[rownames(ZmOsNS),"ZmOs_dN"] = ZmOsNS$ZmOs_dN
ZmOs_comp$ZmOs_dS = NA
ZmOs_comp[rownames(ZmOsNS),"ZmOs_dS"] = ZmOsNS$ZmOs_dS



ggplot(
  ZmOs_comp,
  aes(x = `100TPMcat`, y = ZmOs_dN / ZmOs_dS)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.25) +
  scale_y_log10() +
  theme_classic() +
  labs(
    x = "Haploid expression category",
    y = "dN/dS (log scale)"
  )

pairwise.wilcox.test(ZmOs_comp$ZmOs_dN / ZmOs_comp$ZmOs_dS, ZmOs_comp$"30TPMcat", p.adjust = 'holm')  # p-values
pairwise.wilcox.test(ZmOs_comp$ZmOs_dN / ZmOs_comp$ZmOs_dS, ZmOs_comp$"100TPMcat", p.adjust = 'holm')  # p-values





#Rice CoGe: https://genomevolution.org/coge/SynMap.pl 
#Oryza sativa japonica (Rice) (v7, id16890): masked

    RiceNSdata = read.table('Rice_Oryza_longistaminata_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    RiceNSdata = cbind(RiceNSdata, sub("[.].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[|].*", '', sub(".*OL", "OL", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
    colnames(RiceNSdata)[13:15] = c('RiceGene','longiGene','Identity')

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
RiceNSdata = RiceNSdata[-which(as.numeric(RiceNSdata[,1]) == 0),] # Remove genes with 0 synonymous subs #don't use for Lp or Bd
RiceNS = data.frame(OsOl_dN = as.numeric(RiceNSdata[,2]), OsOl_dS = as.numeric(RiceNSdata[,1]))
rownames(RiceNS) = RiceNSdata[,13]

rownames(ZmOs_comp) = ZmOs_comp$Os_G

RiceNS = RiceNS[rownames(RiceNS) %in% rownames(ZmOs_comp),]

ZmOs_comp$OsOl_dN = NA
ZmOs_comp[rownames(RiceNS),"OsOl_dN"] = RiceNS$OsOl_dN
ZmOs_comp$OsOl_dS = NA
ZmOs_comp[rownames(RiceNS),"OsOl_dS"] = RiceNS$OsOl_dS




    RiceNSdata = read.table('Rice_Leersia_perrieri_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    RiceNSdata = cbind(RiceNSdata, sub("[.].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[|].*", '', sub(".*CDS:", "", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
    colnames(RiceNSdata)[13:15] = c('RiceGene','LeersiaGene','Identity')
   
   
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
RiceNS = data.frame(OsLp_dN = as.numeric(RiceNSdata[,2]), OsLp_dS = as.numeric(RiceNSdata[,1]))
rownames(RiceNS) = RiceNSdata[,13]

RiceNS = RiceNS[rownames(RiceNS) %in% rownames(ZmOs_comp),]

ZmOs_comp$OsLp_dN = NA
ZmOs_comp[rownames(RiceNS),"OsLp_dN"] = RiceNS$OsLp_dN
ZmOs_comp$OsLp_dS = NA
ZmOs_comp[rownames(RiceNS),"OsLp_dS"] = RiceNS$OsLp_dS


rownames(ZmOs_comp) = ZmOs_comp$Zm_G
save(ZmOs_comp, file = "ZmOs_comp.RData")



library(dplyr)
library(ggplot2)
library(tidyr)

# Identify maize and rice columns
zm_cols <- grep("^Zm_", colnames(ZmOs_comp[,1:15]), value = TRUE)
os_cols <- grep("^Os_", colnames(ZmOs_comp[,1:15]), value = TRUE)

# Check the columns being used
zm_cols
os_cols

# Transform normalized mean TPM expression
expr_log <- log10(
  as.matrix(ZmOs_comp[, c(zm_cols, os_cols)]) + 100
)

# Calculate Pearson correlations across genes
cor_all <- cor(
  expr_log,
  use = "pairwise.complete.obs",
  method = "pearson"
)

# Extract only cross-species correlations:
# Rice stages on rows, maize stages on columns
cor_Zm_Os <- cor_all[os_cols, zm_cols]

# Inspect the correlation matrix
round(cor_Zm_Os, 3)


# Convert correlation matrix to long format
cor_long <- as.data.frame(as.table(cor_Zm_Os)) %>%
  rename(
    Rice_stage = Var1,
    Maize_stage = Var2,
    Correlation = Freq
  ) %>%
  mutate(
    Rice_stage = factor(Rice_stage, levels = os_cols),
    Maize_stage = factor(Maize_stage, levels = zm_cols)
  )

# Create cleaner stage labels
rice_labels <- sub("^Os_", "", os_cols)
maize_labels <- sub("^Zm_", "", zm_cols)

# Your preferred palette
heatmap_colors <- colorRampPalette(
  c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020")
)(100)


plot_Zm_Os_correlation <- function(
    cor_long,
    color_scale = "fixed"
) {

  if (!color_scale %in% c("fixed", "minmax")) {
    stop("color_scale must be 'fixed' or 'minmax'")
  }

  # Color limits
  if (color_scale == "fixed") {

    color_limits <- c(-1, 1)

  } else {

    color_limits <- range(
      cor_long$Correlation,
      na.rm = TRUE
    )

    if (diff(color_limits) == 0) {
      color_limits <- color_limits + c(-0.001, 0.001)
    }
  }

  ggplot(
    cor_long,
    aes(
      x = Maize_stage,
      y = Rice_stage,
      fill = Correlation
    )
  ) +

    geom_tile(
      color = "white",
      linewidth = 0.8
    ) +

    geom_text(
      aes(label = sprintf("%.2f", Correlation)),
      size = 4,
      color = "black"
    ) +

    scale_fill_gradientn(
      colors = heatmap_colors,
      limits = color_limits,
      na.value = "grey90",
      name = "Pearson\ncorrelation"
    ) +

    # Maize: original order, left to right
    scale_x_discrete(
      limits = zm_cols,
      labels = maize_labels,
      drop = FALSE
    ) +

    # Rice: first stage at bottom, latest stage at top
    scale_y_discrete(
      limits = os_cols,
      labels = rice_labels,
      drop = FALSE
    ) +

    coord_equal() +

    labs(
      title = "Maize vs. Rice expression correlation",
      subtitle = paste(
        "Pearson correlation of log10(mean TPM + 100)",
        "| Color scale:", color_scale
      ),
      x = "Maize developmental stage",
      y = "Rice developmental stage"
    ) +

    theme_classic() +

    theme(
      axis.text = element_text(
        size = 11,
        color = "black"
      ),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        size = 14,
        face = "bold"
      ),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
}


# Fixed correlation scale (-1 to 1)
p_cor_fixed <- plot_Zm_Os_correlation(
  cor_long,
  color_scale = "fixed"
)

p_cor_fixed

# Min-max scale
p_cor_minmax <- plot_Zm_Os_correlation(
  cor_long,
  color_scale = "minmax"
)

p_cor_minmax