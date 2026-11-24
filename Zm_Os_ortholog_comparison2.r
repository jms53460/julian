#Zm_Os_ortholog_comparison2
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
ZmNS = as.numeric(ZmSb[,2]) / as.numeric(ZmSb[,1])
names(ZmNS) = ZmSb[,13]


#Downloaded table with all gene model versions: https://download.maizegdb.org/Pan-genes/Pan-Zea/
PanGeneTab = read.delim("pan-zea.v4.pan-genes_table.tsv")
V5toV4tab = PanGeneTab[,6:7]
V5toV4tab = V5toV4tab[which(rowSums(V5toV4tab == "NONE") == 0),]
V5toV4tab$Zm00001d = sub("_.*", "", V5toV4tab$Zm00001d)
V5toV4tab$Zm00001eb = sub("_.*", "", V5toV4tab$Zm00001eb)

vConvert = V5toV4tab[which(V5toV4tab$Zm00001d %in% names(ZmNS)),]
rownames(vConvert) = vConvert[,1]

length(ZmNS)
#[1] 14892
ZmNSv5 = ZmNS[names(ZmNS) %in% vConvert[,1]]
names(ZmNSv5) = vConvert[names(ZmNSv5),2]
length(ZmNSv5)
#[1] 13194
ZmNSv5 = ZmNSv5[-which(is.na(ZmNSv5))]

TPM = 30
OsOnly = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Os > TPM & ZmOs_comp$maxMono_Zm <=TPM),])
ZmOnly = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm > TPM & ZmOs_comp$maxMono_Os <=TPM),])
BothOsZm = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm > TPM & ZmOs_comp$maxMono_Os >TPM),])
NeitherOsZm = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm <= TPM & ZmOs_comp$maxMono_Os <=TPM),])

ZmOs_comp$cat = rownames(ZmOs_comp)
ZmOs_comp$cat[ZmOs_comp$cat %in% OsOnly] = "OsOnly"
ZmOs_comp$cat[ZmOs_comp$cat %in% ZmOnly] = "ZmOnly"
ZmOs_comp$cat[ZmOs_comp$cat %in% BothOsZm] = "BothOsZm"
ZmOs_comp$cat[ZmOs_comp$cat %in% NeitherOsZm] = "NeitherOsZm"
table(ZmOs_comp$cat)

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

ZmOs_comp2 = ZmOs_comp[,20:22]

ZmOs_comp2 = ZmOs_comp2[ZmOs_comp2$Zm_G %in% names(ZmNSv5),]
ZmOs_comp2$ZmSb_dNdS = ZmNSv5[ZmOs_comp2$Zm_G]
ZmOs_comp2$cat = factor(ZmOs_comp2$cat, levels=c("NeitherOsZm","OsOnly","ZmOnly","BothOsZm"))

bxp_info <- boxplot(ZmOs_comp2$ZmSb_dNdS ~ ZmOs_comp2$cat, plot = FALSE)
boxplot(ZmOs_comp2$ZmSb_dNdS ~ ZmOs_comp2$cat, ylim = c(0,1), pch=19, 
main=paste("Maize vs. Sorghum dN/dS ", TPM, " TPM", sep=""), ylab="dN/dS", 
xlab=paste("Expressed >", TPM, " TPM post SGT in:", sep=""), 
names = paste0(bxp_info$names, "\n(n = ", bxp_info$n, ")"))

pairwise.wilcox.test(ZmOs_comp2$ZmSb_dNdS, ZmOs_comp2$cat, p.adjust = 'holm')  # p-values
table(ZmOs_comp2$cat)
         NeitherOsZm OsOnly ZmOnly
OsOnly   0.002       -      -     
ZmOnly   0.023       0.846  -     
BothOsZm 8.8e-09     0.064  0.121
NeitherOsZm      OsOnly      ZmOnly    BothOsZm 
       4328         275         134         212 #100 TPM


         NeitherOsZm OsOnly  ZmOnly 
OsOnly   0.78152     -       -      
ZmOnly   0.34054     0.78152 -      
BothOsZm 2.9e-09     0.00019 0.09808
NeitherOsZm      OsOnly      ZmOnly    BothOsZm 
       3739         512         222         476 #50 TPM




#	Zea mays B73 Ab10 (Zea mays B73 Reference NAM 5.0): Genome fasta file from Ensembl plants (v5.0, id64521): unmasked
#	Zea mays subsp. mexicana (vZx-PI566673-REFERENCE-YAN-1.0, id36526): unmasked
ZmTeo = read.table('Zm_v5_teosinte_dNdS.tab', skip = 2, sep = '\t', as.is = T)
ZmTeo = cbind(ZmTeo, sub("_.*", "", sub(".*CDS:", "", ZmTeo[,8])), sub("_.*", "", sub(".*Zx", "Zx", ZmTeo[,4])), as.numeric(sub('.*[|]','',ZmTeo[,4])))
colnames(ZmTeo)[13:15] = c('Zm_G','Teo_G','Identity')

## Remove putative orthologs that are not reciprocal best hits ##
removes = rep(FALSE, nrow(ZmTeo))
for (sgene in as.character(unique(ZmTeo[duplicated(ZmTeo[,13]),13]))) {
	tr = which(ZmTeo[,13] == sgene)
	removes[tr[rank(-ZmTeo[tr,15]) != 1]] = TRUE
}
for (mgene in as.character(unique(ZmTeo[duplicated(ZmTeo[,14]),14]))) {
	tr = which(ZmTeo[,14] == mgene)
	removes[tr[rank(-ZmTeo[tr,15]) != 1]] = TRUE
}
ZmTeo = ZmTeo[!removes,]
##

ZmTeo = ZmTeo[ZmTeo[,15] >= 80,]  # Require at least 80% identity between orthologs when calculating dNdS
ZmTeo = ZmTeo[-which(is.na(as.numeric(ZmTeo[,1]))),] # Remove genes that lack calculated synonymous and non-synonymous rates or are undefined
ZmTeo = ZmTeo[-which(as.numeric(ZmTeo[,1]) == 0),] # Remove genes with 0 synonymous subs 
ZmNSteo = as.numeric(ZmTeo[,2]) / as.numeric(ZmTeo[,1])
names(ZmNSteo) = ZmTeo[,13]
length(ZmNSteo)
#[1] 16654

TPM = 50
OsOnly = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Os > TPM & ZmOs_comp$maxMono_Zm <=TPM),])
ZmOnly = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm > TPM & ZmOs_comp$maxMono_Os <=TPM),])
BothOsZm = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm > TPM & ZmOs_comp$maxMono_Os >TPM),])
NeitherOsZm = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm <= TPM & ZmOs_comp$maxMono_Os <=TPM),])

ZmOs_comp$cat = rownames(ZmOs_comp)
ZmOs_comp$cat[ZmOs_comp$cat %in% OsOnly] = "OsOnly"
ZmOs_comp$cat[ZmOs_comp$cat %in% ZmOnly] = "ZmOnly"
ZmOs_comp$cat[ZmOs_comp$cat %in% BothOsZm] = "BothOsZm"
ZmOs_comp$cat[ZmOs_comp$cat %in% NeitherOsZm] = "NeitherOsZm"
table(ZmOs_comp$cat)

ZmOs_comp$Os_G = ZmOs_NS$Os_G
ZmOs_comp$Zm_G = ZmOs_NS$Zm_G

ZmOs_comp2 = ZmOs_comp[,20:22]

ZmOs_comp2 = ZmOs_comp2[ZmOs_comp2$Zm_G %in% names(ZmNSteo),]
ZmOs_comp2$ZmTeo_dNdS = ZmNSteo[ZmOs_comp2$Zm_G]
ZmOs_comp2$cat = factor(ZmOs_comp2$cat, levels=c("NeitherOsZm","OsOnly","ZmOnly","BothOsZm"))

boxplot(ZmOs_comp2$ZmTeo_dNdS ~ ZmOs_comp2$cat, ylim = c(0,1), pch=19, main=paste("Maize vs. Teosinte dN/dS ", TPM, " TPM", sep=""), ylab="dN/dS", xlab=paste("Expressed >", TPM, " TPM post SGT in:", sep=""))

bxp_info <- boxplot(ZmOs_comp2$ZmTeo_dNdS ~ ZmOs_comp2$cat, plot = FALSE)
boxplot(ZmOs_comp2$ZmTeo_dNdS ~ ZmOs_comp2$cat, ylim = c(0,1), pch=19, 
main=paste("Maize vs. Teosinte dN/dS ", TPM, " TPM", sep=""), ylab="dN/dS", 
xlab=paste("Expressed >", TPM, " TPM post SGT in:", sep=""), 
names = paste0(bxp_info$names, "\n(n = ", bxp_info$n, ")"))

pairwise.wilcox.test(ZmOs_comp2$ZmTeo_dNdS, ZmOs_comp2$cat, p.adjust = 'holm')  # p-values
table(ZmOs_comp2$cat)

         NeitherOsZm OsOnly  ZmOnly
OsOnly   1.0000      -       -     
ZmOnly   2.2e-05     0.0011  -     
BothOsZm 3.1e-08     8.5e-05 1.0000
NeitherOsZm      OsOnly      ZmOnly    BothOsZm 
       3614         223         116         159     #100 TPM

         NeitherOsZm OsOnly  ZmOnly 
OsOnly   0.66106     -       -      
ZmOnly   0.02053     0.13760 -      
BothOsZm 4e-08       0.00058 0.64194
NeitherOsZm      OsOnly      ZmOnly    BothOsZm 
       3161         405         177         369     #50 TPM



#Rice CoGe: https://genomevolution.org/coge/SynMap.pl 
#Oryza sativa japonica (Rice) (v7, id16890): masked

    RiceNSdata = read.table('Rice_Leersia_perrieri_CoGe.tab', skip = 2, sep = '\t', as.is = T)
    RiceNSdata = cbind(RiceNSdata, sub("[.].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[|].*", '', sub(".*CDS:", "", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
    colnames(RiceNSdata)[13:15] = c('RiceGene','LeersiaGene','Identity')
   
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
RiceNS = as.numeric(RiceNSdata[,2]) / as.numeric(RiceNSdata[,1])
names(RiceNS) = RiceNSdata[,13]



TPM = 100
OsOnly = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Os > TPM & ZmOs_comp$maxMono_Zm <=TPM),])
ZmOnly = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm > TPM & ZmOs_comp$maxMono_Os <=TPM),])
BothOsZm = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm > TPM & ZmOs_comp$maxMono_Os >TPM),])
NeitherOsZm = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm <= TPM & ZmOs_comp$maxMono_Os <=TPM),])

ZmOs_comp$cat = rownames(ZmOs_comp)
ZmOs_comp$cat[ZmOs_comp$cat %in% OsOnly] = "OsOnly"
ZmOs_comp$cat[ZmOs_comp$cat %in% ZmOnly] = "ZmOnly"
ZmOs_comp$cat[ZmOs_comp$cat %in% BothOsZm] = "BothOsZm"
ZmOs_comp$cat[ZmOs_comp$cat %in% NeitherOsZm] = "NeitherOsZm"
table(ZmOs_comp$cat)

ZmOs_comp2 = ZmOs_comp[,20:22]
rownames(ZmOs_comp2) = ZmOs_comp2$Os_G
ZmOs_comp2 = ZmOs_comp2[ZmOs_comp2$Os_G %in% names(RiceNS),]
ZmOs_comp2$OsLp_dNdS = RiceNS[ZmOs_comp2$Os_G]
ZmOs_comp2$cat = factor(ZmOs_comp2$cat, levels=c("NeitherOsZm","ZmOnly","OsOnly","BothOsZm"))

bxp_info <- boxplot(ZmOs_comp2$OsLp_dNdS ~ ZmOs_comp2$cat, plot = FALSE)
boxplot(ZmOs_comp2$OsLp_dNdS ~ ZmOs_comp2$cat, ylim = c(0,1), pch=19, 
main=paste("Rice vs. O. longistaminata dN/dS ", TPM, " TPM", sep=""), ylab="dN/dS", 
xlab=paste("Expressed >", TPM, " TPM post SGT in:", sep=""), 
names = paste0(bxp_info$names, "\n(n = ", bxp_info$n, ")"))

pairwise.wilcox.test(ZmOs_comp2$OsLp_dNdS, ZmOs_comp2$cat, p.adjust = 'holm')  # p-values
table(ZmOs_comp2$cat)

L. perrieri
         NeitherOsZm ZmOnly OsOnly
ZmOnly   0.067       -      -     
OsOnly   0.023       0.793  -     
BothOsZm 7.3e-07     0.142  0.067 
NeitherOsZm      ZmOnly      OsOnly    BothOsZm 
       4609         143         294         216     #100 TPM

         NeitherOsZm ZmOnly OsOnly 
ZmOnly   0.16        -      -      
OsOnly   0.82        0.17   -      
BothOsZm 2.1e-09     0.17   2.7e-06
NeitherOsZm      ZmOnly      OsOnly    BothOsZm 
       3993         229         548         492     #50 TPM


#Rice CoGe: https://genomevolution.org/coge/SynMap.pl 
#Oryza sativa japonica (Rice) (v7, id16890): masked
#Probably used IR64

    RiceNSdata = read.table('Rice_japonica_indica_dNdS.tab', skip = 2, sep = '\t', as.is = T)
    RiceNSdata = cbind(RiceNSdata, sub("[.].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[.].*", '', sub(".*Osir", "Osir", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
    colnames(RiceNSdata)[13:15] = c('JaponicaGene','IndicaGene','Identity')

    RiceNSdata = read.table('Rice_japonica_indicaR498_dNdS.tab', skip = 2, sep = '\t', as.is = T)
    RiceNSdata = cbind(RiceNSdata, sub("[.].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[.].*", '', sub(".*Osir", "Osir", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
    colnames(RiceNSdata)[13:15] = c('JaponicaGene','IndicaGene','Identity')

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
RiceNS = as.numeric(RiceNSdata[,2]) / as.numeric(RiceNSdata[,1])
names(RiceNS) = RiceNSdata[,13]
dim(RiceNSdata)                                                             
#[1] 18223    15

TPM = 100
OsOnly = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Os > TPM & ZmOs_comp$maxMono_Zm <=TPM),])
ZmOnly = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm > TPM & ZmOs_comp$maxMono_Os <=TPM),])
BothOsZm = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm > TPM & ZmOs_comp$maxMono_Os >TPM),])
NeitherOsZm = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm <= TPM & ZmOs_comp$maxMono_Os <=TPM),])

ZmOs_comp$cat = rownames(ZmOs_comp)
ZmOs_comp$cat[ZmOs_comp$cat %in% OsOnly] = "OsOnly"
ZmOs_comp$cat[ZmOs_comp$cat %in% ZmOnly] = "ZmOnly"
ZmOs_comp$cat[ZmOs_comp$cat %in% BothOsZm] = "BothOsZm"
ZmOs_comp$cat[ZmOs_comp$cat %in% NeitherOsZm] = "NeitherOsZm"
table(ZmOs_comp$cat)

ZmOs_comp2 = ZmOs_comp[,20:22]
rownames(ZmOs_comp2) = ZmOs_comp2$Os_G
ZmOs_comp2 = ZmOs_comp2[ZmOs_comp2$Os_G %in% names(RiceNS),]
ZmOs_comp2$JpIn_dNdS = RiceNS[ZmOs_comp2$Os_G]
ZmOs_comp2$cat = factor(ZmOs_comp2$cat, levels=c("NeitherOsZm","ZmOnly","OsOnly","BothOsZm"))

bxp_info <- boxplot(ZmOs_comp2$JpIn_dNdS ~ ZmOs_comp2$cat, plot = FALSE)
boxplot(ZmOs_comp2$JpIn_dNdS ~ ZmOs_comp2$cat, ylim = c(0,1), pch=19, 
main=paste("Rice japonica vs. indica dN/dS ", TPM, " TPM", sep=""), ylab="dN/dS", 
xlab=paste("Expressed >", TPM, " TPM post SGT in:", sep=""), 
names = paste0(bxp_info$names, "\n(n = ", bxp_info$n, ")"))

pairwise.wilcox.test(ZmOs_comp2$JpIn_dNdS, ZmOs_comp2$cat, p.adjust = 'holm')  # p-values
table(ZmOs_comp2$cat)





    RiceNSdata = read.table('Rice_NP_KT_dNdS.tab', skip = 2, sep = '\t', as.is = T)
    RiceNSdata = cbind(RiceNSdata, sub("[.].*", "", sub(".*LOC", "LOC", RiceNSdata[,4])), sub("[.].*", '', sub(".*OsK", "OsK", RiceNSdata[,8])), as.numeric(sub('.*[|]','',RiceNSdata[,4])))
    colnames(RiceNSdata)[13:15] = c('NipponbareGene','KitaakeGene','Identity')


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
RiceNS = as.numeric(RiceNSdata[,2]) / as.numeric(RiceNSdata[,1])
names(RiceNS) = RiceNSdata[,13]
dim(RiceNSdata) 
#[1] 6102   15

TPM = 100
OsOnly = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Os > TPM & ZmOs_comp$maxMono_Zm <=TPM),])
ZmOnly = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm > TPM & ZmOs_comp$maxMono_Os <=TPM),])
BothOsZm = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm > TPM & ZmOs_comp$maxMono_Os >TPM),])
NeitherOsZm = rownames(ZmOs_comp[which(ZmOs_comp$maxMono_Zm <= TPM & ZmOs_comp$maxMono_Os <=TPM),])

ZmOs_comp$cat = rownames(ZmOs_comp)
ZmOs_comp$cat[ZmOs_comp$cat %in% OsOnly] = "OsOnly"
ZmOs_comp$cat[ZmOs_comp$cat %in% ZmOnly] = "ZmOnly"
ZmOs_comp$cat[ZmOs_comp$cat %in% BothOsZm] = "BothOsZm"
ZmOs_comp$cat[ZmOs_comp$cat %in% NeitherOsZm] = "NeitherOsZm"
table(ZmOs_comp$cat)

ZmOs_comp$Os_G = ZmOs_NS$Os_G
ZmOs_comp$Zm_G = ZmOs_NS$Zm_G

ZmOs_comp2 = ZmOs_comp[,20:22]
rownames(ZmOs_comp2) = ZmOs_comp2$Os_G
ZmOs_comp2 = ZmOs_comp2[ZmOs_comp2$Os_G %in% names(RiceNS),]
ZmOs_comp2$NpKt_dNdS = RiceNS[ZmOs_comp2$Os_G]
ZmOs_comp2$cat = factor(ZmOs_comp2$cat, levels=c("NeitherOsZm","ZmOnly","OsOnly","BothOsZm"))
table(ZmOs_comp2$cat)
#NeitherOsZm      ZmOnly      OsOnly    BothOsZm 
#        901          15          45          38

bxp_info <- boxplot(ZmOs_comp2$NpKt_dNdS ~ ZmOs_comp2$cat, plot = FALSE)
boxplot(ZmOs_comp2$NpKt_dNdS ~ ZmOs_comp2$cat, ylim = c(0,1), pch=19, 
main=paste("Rice Nipponbare vs. Kitaake dN/dS ", TPM, " TPM", sep=""), ylab="dN/dS", 
xlab=paste("Expressed >", TPM, " TPM post SGT in:", sep=""), 
names = paste0(bxp_info$names, "\n(n = ", bxp_info$n, ")"))

pairwise.wilcox.test(ZmOs_comp2$NpKt_dNdS, ZmOs_comp2$cat, p.adjust = 'holm')  # p-values


