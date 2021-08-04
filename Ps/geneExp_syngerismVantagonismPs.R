load("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/lfc/output/PsyringaeLFC_kallisto.RData")

lib <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/initial_files/P-syringae_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/initial_files/Psraw_counts_kallisoNCBI.txt",sep="\t",header=TRUE)
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/initial_files/Ps_diffExp_all.csv",header=TRUE,sep=",")

library(DESeq2)
library(ImpulseDE2)

#Obtain only significantly differentially expressed genes
impulse_resultsFull <- readRDS("diffExp/output/impulseDiffGenesFullvsMono_kallisto.rds")
impulse_resultsBt <- readRDS("diffExp/output/impulseDiffGenesPsBtvsMono_kallisto.rds")
impulse_resultsCv <- readRDS("diffExp/output/impulseDiffGenesPsCvvsMono_kallisto.rds")
#Note: you need to load library(ImpulseDE2) before loading impulseRDS object otherwise there's issues with S4 objects.

Fullgenes <- impulse_resultsFull$vecDEGenes
PsBtgenes <- impulse_resultsBt$vecDEGenes
PsCvgenes <- impulse_resultsCv$vecDEGenes

allSigGenes <- c(Fullgenes,PsBtgenes,PsCvgenes)
allSigGenes.final <- unique(allSigGenes)

#If you want to check if these are unique to only the Full community differentially expressed genes
sigGenesFullOnly <- readRDS("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/diffExp/output/sigGenesFullOnly_Ps_kallisto.rds")

#######Synergistic gene expression#######

lfc_pairwiseCombined_Int <- 2^(lfc_PsBtInt_comparedtoMonoInt.fc) + 2^(lfc_PsCvInt_comparedtoMonoInt.fc)
lfc_pairwiseCombined_25 <- 2^(lfc_PsBt25_comparedtoMono25.fc) + 2^(lfc_PsCv25_comparedtoMono25.fc)
lfc_pairwiseCombined_30 <- 2^(lfc_PsBt30_comparedtoMono30.fc) + 2^(lfc_PsCv30_comparedtoMono30.fc)
lfc_pairwiseCombined_35 <- 2^(lfc_PsBt35_comparedtoMono35.fc) + 2^(lfc_PsCv35_comparedtoMono35.fc)
lfc_pairwiseCombined_40 <- 2^(lfc_PsBt40_comparedtoMono40.fc) + 2^(lfc_PsCv40_comparedtoMono40.fc)
lfc_pairwiseCombined_45 <- 2^(lfc_PsBt45_comparedtoMono45.fc) + 2^(lfc_PsCv45_comparedtoMono45.fc)

lfc_Full_vs_pwCombined_Int <- (2^(lfc_3memInt_comparedtoMonoInt.fc))/lfc_pairwiseCombined_Int
lfc_Full_vs_pwCombined_25 <- (2^(lfc_3mem25_comparedtoMono25.fc))/lfc_pairwiseCombined_25
lfc_Full_vs_pwCombined_30 <- (2^(lfc_3mem30_comparedtoMono30.fc))/lfc_pairwiseCombined_30
lfc_Full_vs_pwCombined_35 <- (2^(lfc_3mem35_comparedtoMono35.fc))/lfc_pairwiseCombined_35
lfc_Full_vs_pwCombined_40 <- (2^(lfc_3mem40_comparedtoMono40.fc))/lfc_pairwiseCombined_40
lfc_Full_vs_pwCombined_45 <- (2^(lfc_3mem45_comparedtoMono45.fc))/lfc_pairwiseCombined_45

#Create side by side for filtering low Full community lfc values.
lfc_Full_synergyInt <- as.data.frame(cbind(lfc_3memInt_comparedtoMonoInt.fc,lfc_Full_vs_pwCombined_Int))
row.names(lfc_Full_synergyInt) <- row.names(lfc_3memInt_comparedtoMonoInt)
#Filter only significant genes
lfc_Full_synergyInt <- lfc_Full_synergyInt[which(row.names(lfc_Full_synergyInt) %in% allSigGenes.final),]

lfc_Full_synergy25 <- as.data.frame(cbind(lfc_3mem25_comparedtoMono25.fc,lfc_Full_vs_pwCombined_25))
row.names(lfc_Full_synergy25) <- row.names(lfc_3mem25_comparedtoMono25)
#Filter only significant genes
lfc_Full_synergy25 <- lfc_Full_synergy25[which(row.names(lfc_Full_synergy25) %in% allSigGenes.final),]

lfc_Full_synergy30 <- as.data.frame(cbind(lfc_3mem30_comparedtoMono30.fc,lfc_Full_vs_pwCombined_30))
row.names(lfc_Full_synergy30) <- row.names(lfc_3mem30_comparedtoMono30)
#Filter only significant genes
lfc_Full_synergy30 <- lfc_Full_synergy30[which(row.names(lfc_Full_synergy30) %in% allSigGenes.final),]

lfc_Full_synergy35 <- as.data.frame(cbind(lfc_3mem35_comparedtoMono35.fc,lfc_Full_vs_pwCombined_35))
row.names(lfc_Full_synergy35) <- row.names(lfc_3mem35_comparedtoMono35)
#Filter only significant genes
lfc_Full_synergy35 <- lfc_Full_synergy35[which(row.names(lfc_Full_synergy35) %in% allSigGenes.final),]

lfc_Full_synergy40 <- as.data.frame(cbind(lfc_3mem40_comparedtoMono40.fc,lfc_Full_vs_pwCombined_40))
row.names(lfc_Full_synergy40) <- row.names(lfc_3mem40_comparedtoMono40)
#Filter only significant genes
lfc_Full_synergy40 <- lfc_Full_synergy40[which(row.names(lfc_Full_synergy40) %in% allSigGenes.final),]

lfc_Full_synergy45 <- as.data.frame(cbind(lfc_3mem45_comparedtoMono45.fc,lfc_Full_vs_pwCombined_45))
row.names(lfc_Full_synergy45) <- row.names(lfc_3mem45_comparedtoMono45)
#Filter only significant genes
lfc_Full_synergy45 <- lfc_Full_synergy45[which(row.names(lfc_Full_synergy45) %in% allSigGenes.final),]

#Filter out low lfc values from dataframe
lfc_Full_synergyInt_filt <- lfc_Full_synergyInt[lfc_Full_synergyInt$lfc_3memInt_comparedtoMonoInt.fc>0,]
lfc_Full_synergy25_filt <- lfc_Full_synergy25[lfc_Full_synergy25$lfc_3mem25_comparedtoMono25.fc>0,]
lfc_Full_synergy30_filt <- lfc_Full_synergy30[lfc_Full_synergy30$lfc_3mem30_comparedtoMono30.fc>0,]
lfc_Full_synergy35_filt <- lfc_Full_synergy35[lfc_Full_synergy35$lfc_3mem35_comparedtoMono35.fc>0,]
lfc_Full_synergy40_filt <- lfc_Full_synergy40[lfc_Full_synergy40$lfc_3mem40_comparedtoMono40.fc>0,]
lfc_Full_synergy45_filt <- lfc_Full_synergy45[lfc_Full_synergy45$lfc_3mem45_comparedtoMono45.fc>0,]

#Filter lfc >1.5 comparing Full community to additive pairwise communities
lfc_Full_synergyInt_final <- lfc_Full_synergyInt_filt[lfc_Full_synergyInt_filt$lfc_Full_vs_pwCombined_Int>1.2,]
lfc_Full_synergy25_final <- lfc_Full_synergy25_filt[lfc_Full_synergy25_filt$lfc_Full_vs_pwCombined_25>1.2,]
lfc_Full_synergy30_final <- lfc_Full_synergy30_filt[lfc_Full_synergy30_filt$lfc_Full_vs_pwCombined_30>1.2,]
lfc_Full_synergy35_final <- lfc_Full_synergy35_filt[lfc_Full_synergy35_filt$lfc_Full_vs_pwCombined_35>1.2,]
lfc_Full_synergy40_final <- lfc_Full_synergy40_filt[lfc_Full_synergy40_filt$lfc_Full_vs_pwCombined_40>1.2,]
lfc_Full_synergy45_final <- lfc_Full_synergy45_filt[lfc_Full_synergy45_filt$lfc_Full_vs_pwCombined_45>1.2,]

#Obtain locus tags
#lfc_Full_synergyInt_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_synergyInt_final)),]
#lfc_Full_synergy25_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_synergy25_final)),]
#lfc_Full_synergy30_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_synergy30_final)),]
#lfc_Full_synergy35_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_synergy35_final)),]
#lfc_Full_synergy40_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_synergy40_final)),]
#lfc_Full_synergy45_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_synergy45_final)),]

#Put all locus tags together
#genes_syn <- c(as.character(lfc_Full_synergyInt_finalGenes$Locus_Tag),as.character(lfc_Full_synergy25_finalGenes$Locus_Tag),
#           as.character(lfc_Full_synergy30_finalGenes$Locus_Tag),as.character(lfc_Full_synergy35_finalGenes$Locus_Tag),
#           as.character(lfc_Full_synergy40_finalGenes$Locus_Tag),as.character(lfc_Full_synergy45_finalGenes$Locus_Tag))

genes_syn <- c(as.character(row.names(lfc_Full_synergyInt_final)),as.character(row.names(lfc_Full_synergy25_final)),
           as.character(row.names(lfc_Full_synergy30_final)),as.character(row.names(lfc_Full_synergy35_final)),
           as.character(row.names(lfc_Full_synergy40_final)),as.character(row.names(lfc_Full_synergy45_final)))

#Check for the amount of unique elements
synergyGenes <- as.data.frame(unique(genes_syn))

write.csv(synergyGenes,"synergyGenes_Ps_kallisto.csv")

#Prepare list for figure.
Full_synergistic <- c(nrow(lfc_Full_synergyInt_final),nrow(lfc_Full_synergy25_final),nrow(lfc_Full_synergy30_final),
                      nrow(lfc_Full_synergy35_final),nrow(lfc_Full_synergy40_final),nrow(lfc_Full_synergy45_final))

#######Antagonistic gene expression#######

#Create side by side for filtering low Full community lfc values.
lfc_Full_antagInt <- as.data.frame(cbind(lfc_3memInt_comparedtoMonoInt.fc,lfc_PsBtInt_comparedtoMonoInt.fc,lfc_PsCvInt_comparedtoMonoInt.fc))
row.names(lfc_Full_antagInt) <- row.names(lfc_3memInt_comparedtoMonoInt)
colnames(lfc_Full_antagInt) <- c("Full","A","B")
#Filter only significant genes
lfc_Full_antagInt <- lfc_Full_antagInt[which(row.names(lfc_Full_antagInt) %in% allSigGenes.final),]

lfc_Full_antag25 <- as.data.frame(cbind(lfc_3mem25_comparedtoMono25.fc,lfc_PsBt25_comparedtoMono25.fc,lfc_PsCv25_comparedtoMono25.fc))
row.names(lfc_Full_antag25) <- row.names(lfc_3mem25_comparedtoMono25)
colnames(lfc_Full_antag25) <- c("Full","A","B")
#Filter only significant genes
lfc_Full_antag25 <- lfc_Full_antag25[which(row.names(lfc_Full_antag25) %in% allSigGenes.final),]

lfc_Full_antag30 <- as.data.frame(cbind(lfc_3mem30_comparedtoMono30.fc,lfc_PsBt30_comparedtoMono30.fc,lfc_PsCv30_comparedtoMono30.fc))
row.names(lfc_Full_antag30) <- row.names(lfc_3mem30_comparedtoMono30)
colnames(lfc_Full_antag30) <- c("Full","A","B")
#Filter only significant genes
lfc_Full_antag30 <- lfc_Full_antag30[which(row.names(lfc_Full_antag30) %in% allSigGenes.final),]

lfc_Full_antag35 <- as.data.frame(cbind(lfc_3mem35_comparedtoMono35.fc,lfc_PsBt35_comparedtoMono35.fc,lfc_PsCv35_comparedtoMono35.fc))
row.names(lfc_Full_antag35) <- row.names(lfc_3mem35_comparedtoMono35)
colnames(lfc_Full_antag35) <- c("Full","A","B")
#Filter only significant genes
lfc_Full_antag35 <- lfc_Full_antag35[which(row.names(lfc_Full_antag35) %in% allSigGenes.final),]

lfc_Full_antag40 <- as.data.frame(cbind(lfc_3mem40_comparedtoMono40.fc,lfc_PsBt40_comparedtoMono40.fc,lfc_PsCv40_comparedtoMono40.fc))
row.names(lfc_Full_antag40) <- row.names(lfc_3mem40_comparedtoMono40)
colnames(lfc_Full_antag40) <- c("Full","A","B")
#Filter only significant genes
lfc_Full_antag40 <- lfc_Full_antag40[which(row.names(lfc_Full_antag40) %in% allSigGenes.final),]

lfc_Full_antag45 <- as.data.frame(cbind(lfc_3mem45_comparedtoMono45.fc,lfc_PsBt45_comparedtoMono45.fc,lfc_PsCv45_comparedtoMono45.fc))
row.names(lfc_Full_antag45) <- row.names(lfc_3mem45_comparedtoMono45)
colnames(lfc_Full_antag45) <- c("Full","A","B")
#Filter only significant genes
lfc_Full_antag45 <- lfc_Full_antag45[which(row.names(lfc_Full_antag45) %in% allSigGenes.final),]

#Filter out genes where at least 1.5 fc is observed in either pairwise condition
lfc_Full_antagInt_filt <- lfc_Full_antagInt[(lfc_Full_antagInt$A>0.585 | lfc_Full_antagInt$B>0.585),]
lfc_Full_antag25_filt <- lfc_Full_antag25[(lfc_Full_antag25$A>0.585 | lfc_Full_antag25$B>0.585),]
lfc_Full_antag30_filt <- lfc_Full_antag30[(lfc_Full_antag30$A>0.585 | lfc_Full_antag30$B>0.585),]
lfc_Full_antag35_filt <- lfc_Full_antag35[(lfc_Full_antag35$A>0.585 | lfc_Full_antag35$B>0.585),]
lfc_Full_antag40_filt <- lfc_Full_antag40[(lfc_Full_antag40$A>0.585 | lfc_Full_antag40$B>0.585),]
lfc_Full_antag45_filt <- lfc_Full_antag45[(lfc_Full_antag45$A>0.585 | lfc_Full_antag45$B>0.585),]

#Finalize genes where lfc in Full is <1
lfc_Full_antagInt_final <- lfc_Full_antagInt_filt[lfc_Full_antagInt_filt$Full<0,]
lfc_Full_antag25_final <- lfc_Full_antag25_filt[lfc_Full_antag25_filt$Full<0,]
lfc_Full_antag30_final <- lfc_Full_antag30_filt[lfc_Full_antag30_filt$Full<0,]
lfc_Full_antag35_final <- lfc_Full_antag35_filt[lfc_Full_antag35_filt$Full<0,]
lfc_Full_antag40_final <- lfc_Full_antag40_filt[lfc_Full_antag40_filt$Full<0,]
lfc_Full_antag45_final <- lfc_Full_antag45_filt[lfc_Full_antag45_filt$Full<0,]

#Obtain locus tags
#lfc_Full_antagInt_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_antagInt_final)),]
#lfc_Full_antag25_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_antag25_final)),]
#lfc_Full_antag30_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_antag30_final)),]
#lfc_Full_antag35_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_antag35_final)),]
#lfc_Full_antag40_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_antag40_final)),]
#lfc_Full_antag45_finalGenes <- counts[which(counts$geneID %in% row.names(lfc_Full_antag45_final)),]

#Put all locus tags together
#genes_antag <- c(as.character(lfc_Full_antagInt_finalGenes$Locus_Tag),as.character(lfc_Full_antag25_finalGenes$Locus_Tag),
#           as.character(lfc_Full_antag30_finalGenes$Locus_Tag),as.character(lfc_Full_antag35_finalGenes$Locus_Tag),
#           as.character(lfc_Full_antag40_finalGenes$Locus_Tag),as.character(lfc_Full_antag45_finalGenes$Locus_Tag))

genes_antag <- c(as.character(row.names(lfc_Full_antagInt_final)),as.character(row.names(lfc_Full_antag25_final)),
           as.character(row.names(lfc_Full_antag30_final)),as.character(row.names(lfc_Full_antag35_final)),
           as.character(row.names(lfc_Full_antag40_final)),as.character(row.names(lfc_Full_antag45_final)))

#Check for the amount of unique elements
antagGenes <- as.data.frame(unique(genes_antag))

write.csv(antagGenes,"antagonismGenes_Ps_kallisto.csv")

#Prepare lists for figure

Full_antagonism <- c(nrow(lfc_Full_antagInt_final),nrow(lfc_Full_antag25_final),nrow(lfc_Full_antag30_final),
                      nrow(lfc_Full_antag35_final),nrow(lfc_Full_antag40_final),nrow(lfc_Full_antag45_final))

Full_antagonism_CondA <- c(sum(lfc_Full_antagInt_final$A<0),sum(lfc_Full_antag25_final$A<0),sum(lfc_Full_antag30_final$A<0),
                           sum(lfc_Full_antag35_final$A<0),sum(lfc_Full_antag40_final$A<0),sum(lfc_Full_antag45_final$A<0))

Full_antagonism_CondB <- c(sum(lfc_Full_antagInt_final$B<0),sum(lfc_Full_antag25_final$B<0),sum(lfc_Full_antag30_final$B<0),
                           sum(lfc_Full_antag35_final$B<0),sum(lfc_Full_antag40_final$B<0),sum(lfc_Full_antag45_final$B<0))

Full_antagonism_Full <- Full_antagonism - (Full_antagonism_CondA + Full_antagonism_CondB)

#List for syn and antag
FinalGeneTrends <- data.frame(condition=rep(c("Non-additive Up", "Non-additive Down"), each=6),
                time=rep(c("12.5", "25", "30", "35", "40", "45"),2),
                Genes=c(Full_synergistic,Full_antagonism))

#List for syn and antag categories
FinalGeneTrends_f <- data.frame(condition=rep(c("Non-additive Up", "Non-additive Down (PsBtCv)","Non-additive Down (PsBt)","Non-additive Down (PsCv)"), each=6),
                time=rep(c("12.5", "25", "30", "35", "40", "45"),4),
                Genes=c(Full_synergistic, Full_antagonism_Full, Full_antagonism_CondA,Full_antagonism_CondB))
library(ggplot2)
FinalGeneTrends$trend <- factor(FinalGeneTrends$time,levels = c("12.5", "25", "30", "35", "40", "45"),ordered = TRUE)
FinalGeneTrends_f$trend <- factor(FinalGeneTrends_f$time,levels = c("12.5", "25", "30", "35", "40", "45"),ordered = TRUE)

saveRDS(FinalGeneTrends,"FinalGeneTrends_kallisto.rds")
saveRDS(FinalGeneTrends_f,"FinalGeneTrends_antagCategories_kallisto.rds")

library (ggplot2)

FinalGeneTrends <- readRDS("FinalGeneTrends_kallisto.rds")
FinalGeneTrends_f <- readRDS("FinalGeneTrends_antagCategories_kallisto.rds")

p <- ggplot(data=FinalGeneTrends, aes(x=time,.desc =TRUE, y=Genes, fill=condition)) +
geom_bar(stat="identity", color="black", position=position_stack()) + ylim(0,300) +
  theme_minimal()

FinalGeneTrends_f$condition <- factor(FinalGeneTrends_f$condition,levels = c('Non-additive Up','Non-additive Down (PsBtCv)','Non-additive Down (PsBt)','Non-additive Down (PsCv)'),ordered = TRUE)
p2 <- ggplot(data=FinalGeneTrends_f, aes(x=time,.desc =TRUE, y=Genes, fill=condition)) +
geom_bar(stat="identity", color="black", position='dodge') + ylim(0,70) +
        labs(y="# of genes with non-additive expression", x = "Time (h)") +
        scale_fill_manual(values=c("black", "gray","darkgoldenrod2","burlywood4"))+
        theme(legend.position = "bottom", legend.key.size = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size=3),
        axis.text = element_text(size=6),
        axis.title = element_text(size=8))

ggsave("Ps_synAntag_trends.eps",plot=p2,device="eps",width=8,height=8, units="cm",dpi=300)


ggsave("PsgenesSynvsAntag_kallisto.tif",plot=p,device="tiff",width=15, units="cm",dpi=600)
ggsave("PsgenesSynvsAntagCat_kallisto.tif",plot=p2,device="tiff",width=15, units="cm",dpi=600)

#Obtain most consistent non-additive up

lfc_Full_synergyInt_final[order(lfc_Full_synergyInt_final$lfc_Full_vs_pwCombined_Int),]
lfc_Full_synergy25_final[order(lfc_Full_synergy25_final$lfc_Full_vs_pwCombined_25),]
lfc_Full_synergy30_final[order(lfc_Full_synergy30_final$lfc_Full_vs_pwCombined_30),]
lfc_Full_synergy35_final[order(lfc_Full_synergy35_final$lfc_Full_vs_pwCombined_35),]
lfc_Full_synergy40_final[order(lfc_Full_synergy40_final$lfc_Full_vs_pwCombined_40),]
lfc_Full_synergy45_final[order(lfc_Full_synergy45_final$lfc_Full_vs_pwCombined_45),]
