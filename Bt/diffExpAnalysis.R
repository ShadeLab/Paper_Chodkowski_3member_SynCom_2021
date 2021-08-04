```

###Package Info###

R version 3.5.0 (2018-04-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS: /opt/software/R/3.5.0-iomkl-2018a-X11-20180131/lib64/R/lib/libR.so
LAPACK: /opt/software/R/3.5.0-iomkl-2018a-X11-20180131/lib64/R/modules/lapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils
 [8] datasets  methods   base

other attached packages:
 [1] VennDiagram_1.6.20          futile.logger_1.4.3
 [3] stringr_1.3.1               ImpulseDE2_1.6.0
 [5] DESeq2_1.22.1               SummarizedExperiment_1.12.0
 [7] DelayedArray_0.8.0          BiocParallel_1.16.1
 [9] matrixStats_0.54.0          Biobase_2.42.0
[11] GenomicRanges_1.34.0        GenomeInfoDb_1.18.1
[13] IRanges_2.16.0              S4Vectors_0.20.1
[15] BiocGenerics_0.28.0

loaded via a namespace (and not attached):
 [1] bit64_0.9-7            splines_3.5.0          Formula_1.2-3
 [4] assertthat_0.2.0       latticeExtra_0.6-28    blob_1.1.1
 [7] GenomeInfoDbData_1.2.0 pillar_1.3.0           RSQLite_2.1.1
[10] backports_1.1.2        lattice_0.20-35        glue_1.3.0
[13] digest_0.6.18          RColorBrewer_1.1-2     XVector_0.22.0
[16] checkmate_1.8.5        colorspace_1.3-2       cowplot_0.9.3
[19] htmltools_0.3.6        Matrix_1.2-14          plyr_1.8.4
[22] XML_3.98-1.16          pkgconfig_2.0.2        GetoptLong_0.1.7
[25] genefilter_1.64.0      zlibbioc_1.28.0        purrr_0.2.5
[28] xtable_1.8-3           scales_1.0.0           htmlTable_1.12
[31] tibble_1.4.2           annotate_1.60.0        ggplot2_3.1.0
[34] nnet_7.3-12            lazyeval_0.2.1         survival_2.42-3
[37] magrittr_1.5           crayon_1.3.4           memoise_1.1.0
[40] foreign_0.8-70         tools_3.5.0            data.table_1.11.8
[43] GlobalOptions_0.1.0    formatR_1.5            ComplexHeatmap_1.20.0
[46] locfit_1.5-9.1         munsell_0.5.0          cluster_2.0.7-1
[49] lambda.r_1.2.3         AnnotationDbi_1.44.0   bindrcpp_0.2.2
[52] compiler_3.5.0         rlang_0.3.0.1          RCurl_1.95-4.11
[55] rstudioapi_0.8         circlize_0.4.4         rjson_0.2.20
[58] htmlwidgets_1.3        bitops_1.0-6           base64enc_0.1-3
[61] gtable_0.2.0           DBI_1.0.0              R6_2.3.0
[64] gridExtra_2.3          knitr_1.20             dplyr_0.7.8
[67] bit_1.1-14             bindr_0.1.1            Hmisc_4.1-1
[70] futile.options_1.0.1   shape_1.4.4            stringi_1.2.4
[73] Rcpp_1.0.0             geneplotter_1.60.0     rpart_4.1-13
[76] acepack_1.4.1          tidyselect_0.2.5

library(DESeq2)
library(ImpulseDE2)

#Now we can start Impulse analysis

#Load objects containing size factor and dispersion estimates

dispersions_deseq <- readRDS("diffExp/initial_files/dispersions_deseq_kallistoNCBI.rds")
sizeEst_deseq <- readRDS("diffExp/initial_files/sizeEst_deseq_kallistoNCBI.rds")
diffGenesF <- readRDS("diffExp/initial_files/diffGenesF_kallistoNCBI.rds")

###Now we can prep for Impulse
library(ImpulseDE2)
lib <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/initial_files/B-thailandensis_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/initial_files/Btraw_counts_kallisoNCBI.txt",sep="\t",header=TRUE)

#Remove rRNA, tRNA, and miscRNA
library(rtracklayer)
gff <- readGFF("initial_files/GCA_000012365.1_ASM1236v1_genomic.gff")
CDS <- gff[gff$type=="CDS",]
CDS_parent <- as.data.frame(CDS$Parent)

countsF <- counts[which(counts$target_id %in% CDS_parent$value),]

#Split first column by "-" delimiter so we're left with only locus tags
library(tidyr)
library(dplyr)

countsFf <- countsF %>% separate(target_id, c("gene","Locus"),1, sep = "-",remove=TRUE)
countsFf <- select(countsFf,-c(gene))

#remove non-numerical columns from count matrix
genecounts <- countsFf[2:ncol(countsFf)]
row.names(genecounts) <- countsFf[,1]

#Convert numeric to integers
genecounts[1:ncol(genecounts)] <- lapply(genecounts[1:ncol(genecounts)],as.integer)


#Replace NA values with 0s
genecounts[is.na(genecounts)] <- 0

#Remove 0 count genes
genesZeroC <- rownames(genecounts)[which(rowSums(genecounts) == 0)]
genesRemove <- c(genesZeroC)

#Extract unique identifiers from each condition
lib$libraryName = as.character(lib$libraryName)
genecounts.sort=genecounts[,match(lib$libraryName, names(genecounts))]
library(stringr)
Extract <- c("3mem","Cv","Ps","mono")
keywords <- str_extract(lib$sampleName, paste(Extract,collapse="|"))

################# Monoculture- case only  ########################
mapping <- read.csv(file="initial_files/Bt_diffExp_3mem.csv",header=TRUE,sep=",")

#Get mono only
mapping <- mapping[25:46,]

keys <- which(keywords %in% c("mono"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesFull=diffGenesF[,match(align, colnames(diffGenesF))]
colnames(diffGenesFull) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="case", design$Condition)
diffGenesFull <- as.matrix(diffGenesFull)

names(dispersions_deseq) <- row.names(diffGenesFull)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify transients
#NOTE: matCountData should be raw counts, not normalized orelse you'll be normalizing normalized data.
impulse_results_trends <- runImpulseDE2(matCountData = diffGenesFull,dfAnnotation =design,
boolCaseCtrl = FALSE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_results_trends,"diffExp/output/impulseDiffGenesMonoCaseOnly_kallisto.rds")

################# Only Full community vs control  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/initial_files/Bt_diffExp_3mem.csv",header=TRUE,sep=",")

keys <- which(keywords %in% c("3mem","mono"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesFull=diffGenesF[,match(align, colnames(diffGenesF))]
colnames(diffGenesFull) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="control", design$Condition)
design$Condition <- sub(pattern="exp",replacement="case", design$Condition)
diffGenesFull <- as.matrix(diffGenesFull)

names(dispersions_deseq) <- row.names(diffGenesFull)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify transients
#NOTE: matCountData should be raw counts, not normalized orelse you'll be normalizing normalized data.
impulse_results_trends <- runImpulseDE2(matCountData = diffGenesFull,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_results_trends,"diffExp/output/impulseDiffGenesFullvsMono_kallisto.rds")

################# Only BtCv coculture vs control  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/initial_files/Bt_diffExp_2mem_Chromo.csv",header=TRUE,sep=",")

keys <- which(keywords %in% c("Cv","mono"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesCv=diffGenesF[,match(align, colnames(diffGenesF))]
colnames(diffGenesCv) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="control", design$Condition)
design$Condition <- sub(pattern="exp",replacement="case", design$Condition)
diffGenesCv <- as.matrix(diffGenesCv)

names(dispersions_deseq) <- row.names(diffGenesCv)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify Transients
impulse_resultsCv_trends <- runImpulseDE2(matCountData = diffGenesCv,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsCv_trends,"diffExp/output/impulseDiffGenesBtCvvsMono_kallisto.rds")

################# Only Bt-Ps coculture vs control  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/initial_files/Bt_diffExp_2mem_Pseudo.csv",header=TRUE,sep=",")

keys <- which(keywords %in% c("Ps","mono"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesPs=diffGenesF[,match(align, colnames(diffGenesF))]
colnames(diffGenesPs) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="control", design$Condition)
design$Condition <- sub(pattern="exp",replacement="case", design$Condition)
diffGenesPs <- as.matrix(diffGenesPs)

names(dispersions_deseq) <- row.names(diffGenesPs)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify transients
impulse_resultsPs_trends <- runImpulseDE2(matCountData = diffGenesPs,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsPs_trends,"diffExp/output/impulseDiffGenesBtPsvsMono_kallisto.rds")

########################################

#Create venn diagram of differential gene expression with all coculture conditions compared to monoculture control

diffGenesF <- readRDS("diffExp/initial_files/diffGenesF_kallistoNCBI.rds")

library(ImpulseDE2)

impulse_resultsFull <- readRDS("diffExp/output/impulseDiffGenesFullvsMono_kallisto.rds")
impulse_resultsCv <- readRDS("diffExp/output/impulseDiffGenesBtCvvsMono_kallisto.rds")
impulse_resultsPs <- readRDS("diffExp/output/impulseDiffGenesBtPsvsMono_kallisto.rds")

Fullgenes <- impulse_resultsFull$vecDEGenes
BtCvgenes <- impulse_resultsCv$vecDEGenes
BtPsgenes <- impulse_resultsPs$vecDEGenes


library(VennDiagram)

x <- list(
  BtCvPs = Fullgenes,
  BtCv = BtCvgenes,
  BtPs = BtPsgenes
  )

venn.plot <- venn.diagram(x = list(BtCvPs =Fullgenes,
BtCv=BtCvgenes, BtPs=BtPsgenes),
filename = "ImpulseDE2Summary_comparedToMonoControl.tiff",
    #fill = c("cornflowerblue", "green", "yellow"),
    fill = NULL,
    col = rep("black",3),
    cex = 1.0,
    fontfamily = "serif",
    fontface = "bold",
    cat.cex = 1.5,
    cat.fontfamily = "serif")

#There appear to be 862 differentially expressed genes unique to the full community.
#Are these genes differentially regulated in the Full community compared to the BtCv and BtPs cocultures?

#Functions to pull of genes of interest from venn diagram.
    Intersect <- function (x) {
      # Multiple set version of intersect
      # x is a list
      if (length(x) == 1) {
        unlist(x)
      } else if (length(x) == 2) {
        intersect(x[[1]], x[[2]])
      } else if (length(x) > 2){
        intersect(x[[1]], Intersect(x[-1]))
      }
    }

    Union <- function (x) {
      # Multiple set version of union
      # x is a list
      if (length(x) == 1) {
        unlist(x)
      } else if (length(x) == 2) {
        union(x[[1]], x[[2]])
      } else if (length(x) > 2) {
        union(x[[1]], Union(x[-1]))
      }
    }

    Setdiff <- function (x, y) {
      # Remove the union of the y's from the common x's.
      # x and y are lists of characters.
      xx <- Intersect(x)
      yy <- Union(y)
      setdiff(xx, yy)
    }


genesF <- Setdiff(x = list(Full =Fullgenes),y=list(BtCv=BtCvgenes, BtPs=BtPsgenes))

#Extract these genes by sample matrix

diffGenesFOnly <- diffGenesF[which(row.names(diffGenesF) %in% genesF),]

################# Only BtFull vs BtCv  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/initial_files/Bt_diffExp_3mem_v_BtCv.csv",header=TRUE,sep=",")

keys <- which(keywords %in% c("3mem","Cv"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesFO_BtCv=diffGenesFOnly[,match(align, colnames(diffGenesFOnly))]
colnames(diffGenesFO_BtCv) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="control", design$Condition)
design$Condition <- sub(pattern="exp",replacement="case", design$Condition)
diffGenesFO_BtCv <- as.matrix(diffGenesFO_BtCv)

dispersions_deseq <- readRDS("diffExp/initial_files/dispersions_deseq_kallistoNCBI.rds")
dispersions_deseq <- dispersions_deseq[which(row.names(diffGenesF) %in% genesF)]
names(dispersions_deseq) <- row.names(diffGenesFO_BtCv)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify transients
impulse_resultsFO_BtCv_trends <- runImpulseDE2(matCountData = diffGenesFO_BtCv,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsFO_BtCv_trends,"diffExp/output/impulseUniqueDiffGenesFullvsBtCv_kallisto.rds")

################# Only Bt-Full vs Bt-Ps coculture ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/initial_files/Bt_diffExp_3mem_v_BtPs.csv",header=TRUE,sep=",")

keys <- which(keywords %in% c("3mem","Ps"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesFO_BtPs=diffGenesFOnly[,match(align, colnames(diffGenesFOnly))]
colnames(diffGenesFO_BtPs) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="control", design$Condition)
design$Condition <- sub(pattern="exp",replacement="case", design$Condition)
diffGenesFO_BtPs <- as.matrix(diffGenesFO_BtPs)

dispersions_deseq <- readRDS("diffExp/initial_files/dispersions_deseq_kallistoNCBI.rds")
dispersions_deseq <- dispersions_deseq[which(row.names(diffGenesF) %in% genesF)]
names(dispersions_deseq) <- row.names(diffGenesFO_BtPs)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify transients
impulse_resultsFO_BtPs_trends <- runImpulseDE2(matCountData = diffGenesFO_BtPs,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsFO_BtPs_trends,"diffExp/output/impulseUniqueDiffGenesFullvsBtPs_kallisto.rds")

#################################################################
library(VennDiagram)

impulse_resultsFO_BtCv <- readRDS("diffExp/output/impulseUniqueDiffGenesFullvsBtCv_kallisto.rds")
impulse_resultsFO_BtPs <- readRDS("diffExp/output/impulseUniqueDiffGenesFullvsBtPs_kallisto.rds")

Full_monoCon <- genesF
Full_BtCvCon <- impulse_resultsFO_BtCv$vecDEGenes
Full_BtPsCon <- impulse_resultsFO_BtPs$vecDEGenes

venn.plot <- venn.diagram(x = list(FullvMono =Full_monoCon,
FullvBtCv=Full_BtCvCon, FullvBtPs=Full_BtPsCon),
filename = "ImpulseDE2FulluniqueDiffGenes_kallistoBt.tiff",
col = "transparent",
    fill = c("cornflowerblue", "green", "yellow"),
    #cex = 1.0,
    fontfamily = "serif",
    fontface = "bold",
    #cat.cex = 1.5,
    cat.fontfamily = "serif")

#This venn represents genes that are uniquely differentially regulated in the Full community compared to monoculture. These genes were then tested for significance using BtCv as the control and BtPs as the control.
#So genes in venn diagram show genes that were differentially regulated in the Full community compared to mono and BtCv,mono and BtPs, mono only, and all other conditions.

final_sigGeneFullOnly <- intersect(Full_monoCon,intersect(Full_BtCvCon, Full_BtPsCon))
saveRDS(final_sigGeneFullOnly,"diffExp/output/sigGenesFullOnly_Bt_kallisto.rds")

#Parse out other genes of interest
#overlap <- calculate.overlap(list("FullvMono" =Full_monoCon,
#"FullvBtCv"=Full_BtCvCon, "FullvBtPs"=Full_BtPsCon))

#final_sigGenesBtCv_Mono <- overlap$a2
#final_sigGenesBtPs_Mono <- overlap$a4

#Write these genes to csv
write.csv(final_sigGeneFullOnly, "SigGenesUniqueToBt3mem.csv")

######################Significant gene trends##########################

library(ImpulseDE2)

impulse_resultsFull <- readRDS("diffExp/output/impulseDiffGenesFullvsMono_kallisto.rds")
impulse_resultsCv <- readRDS("diffExp/output/impulseDiffGenesBtCvvsMono_kallisto.rds")
impulse_resultsPs <- readRDS("diffExp/output/impulseDiffGenesBtPsvsMono_kallisto.rds")
impulse_resultsMono <- readRDS("diffExp/output/impulseDiffGenesMonoCaseOnly_kallisto.rds")


load("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/lfc/output/BthailandensisLFC_kallisto.RData")

lfc_analysisFull <- cbind(lfc_3memInt_comparedtoMonoInt.fc,lfc_3mem25_comparedtoMono25.fc,lfc_3mem30_comparedtoMono30.fc,
                          lfc_3mem35_comparedtoMono35.fc,lfc_3mem40_comparedtoMono40.fc,lfc_3mem45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisFull) <- row.names(lfc_3memInt_comparedtoMonoInt)

lfc_analysisBtCv <- cbind(lfc_BtCvInt_comparedtoMonoInt.fc,lfc_BtCv25_comparedtoMono25.fc,lfc_BtCv30_comparedtoMono30.fc,
                          lfc_BtCv35_comparedtoMono35.fc,lfc_BtCv40_comparedtoMono40.fc,lfc_BtCv45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisBtCv) <- row.names(lfc_BtCvInt_comparedtoMonoInt)

lfc_analysisBtPs <- cbind(lfc_BtPsInt_comparedtoMonoInt.fc,lfc_BtPs25_comparedtoMono25.fc,lfc_BtPs30_comparedtoMono30.fc,
                          lfc_BtPs35_comparedtoMono35.fc,lfc_BtPs40_comparedtoMono40.fc,lfc_BtPs45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisBtPs) <- row.names(lfc_BtPsInt_comparedtoMonoInt)

lfc_analysisMono <- cbind(lfc_mono_25withinComparisionToTime0.fc,lfc_mono_30withinComparisionToTime0.fc,lfc_mono_35withinComparisionToTime0.fc,
                          lfc_mono_40withinComparisionToTime0.fc,lfc_mono_45withinComparisionToTime0.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisMono) <- row.names(lfc_mono25_withinComparision_ToTime0)

FullSG <- lfc_analysisFull[match(impulse_resultsFull$vecDEGenes,row.names(lfc_analysisFull)),]
BtCvSG <- lfc_analysisBtCv[match(impulse_resultsCv$vecDEGenes,row.names(lfc_analysisBtCv)),]
BtPsSG <- lfc_analysisBtPs[match(impulse_resultsPs$vecDEGenes,row.names(lfc_analysisBtPs)),]
MonoSG <- lfc_analysisMono[match(impulse_resultsMono$vecDEGenes,row.names(lfc_analysisMono)),]

FullSG_lfc_analysis <- as.data.frame(FullSG)
BtCvSG_lfc_analysis <- as.data.frame(BtCvSG)
BtPsSG_lfc_analysis <- as.data.frame(BtPsSG)
MonoSG_lfc_analysis <- as.data.frame(MonoSG)

library(dplyr)
library(tibble)

#####Full community significant gene trends
#up-regulated
Full_up <- FullSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(Full_up)

#downregulated
Full_down <-  FullSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(Full_down)

#Exponential Phase TP removal. It happens to be the first column in FullSG_lfc_analysis
FullSG_lfc_analysisSP <- FullSG_lfc_analysis[,-1]

#up-regulated in stationary phase
Full_upSP <- FullSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(Full_upSP)

#downregulated in stationary phase
Full_downSP <-  FullSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(Full_downSP)

saveRDS(Full_up,"diffExp/output/Full_up_kallisto.rds")
saveRDS(Full_upSP,"diffExp/output/Full_upSP_kallisto.rds")
saveRDS(Full_down,"diffExp/output/Full_down_kallisto.rds")
saveRDS(Full_downSP,"diffExp/output/Full_downSP_kallisto.rds")

FullgeneTrends <- data.frame(trends=c("Up","Down","UpSp","DownSP","VSP"),
  sigGenes=c(nrow(Full_up),nrow(Full_down),(nrow(Full_upSP)-nrow(Full_up)),(nrow(Full_downSP)-nrow(Full_down)),(nrow(FullSG)-sum(nrow(Full_upSP),nrow(Full_downSP)))))

FullSigGenes <- c(nrow(Full_up),nrow(Full_down),(nrow(Full_upSP)-nrow(Full_up)),(nrow(Full_downSP)-nrow(Full_down)),(nrow(FullSG)-sum(nrow(Full_upSP),nrow(Full_downSP))))

#####BtCv significant gene trends
#up-regulated
BtCv_up <- BtCvSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(BtCv_up)

#downregulated
BtCv_down <-  BtCvSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(BtCv_down)

#Exponential Phase TP removal. It happens to be the first column in FullSG_lfc_analysis
BtCvSG_lfc_analysisSP <- BtCvSG_lfc_analysis[,-1]

#up-regulated in stationary phase
BtCv_upSP <- BtCvSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(BtCv_upSP)

#downregulated in stationary phase
BtCv_downSP <-  BtCvSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(BtCv_downSP)

saveRDS(BtCv_up,"diffExp/output/BtCv_up_kallisto.rds")
saveRDS(BtCv_upSP,"diffExp/output/BtCv_upSP_kallisto.rds")
saveRDS(BtCv_down,"diffExp/output/BtCv_down_kallisto.rds")
saveRDS(BtCv_downSP,"diffExp/output/BtCv_downSP_kallisto.rds")

BtCvgeneTrends <- data.frame(trends=c("Up","Down","UpSp","DownSP","VSP"),
  sigGenes=c(nrow(BtCv_up),nrow(BtCv_down),(nrow(BtCv_upSP)-nrow(BtCv_up)),(nrow(BtCv_downSP)-nrow(BtCv_down)),(nrow(BtCvSG)-sum(nrow(BtCv_upSP),nrow(BtCv_downSP)))))

BtCvSigGenes <- c(nrow(BtCv_up),nrow(BtCv_down),(nrow(BtCv_upSP)-nrow(BtCv_up)),(nrow(BtCv_downSP)-nrow(BtCv_down)),(nrow(BtCvSG)-sum(nrow(BtCv_upSP),nrow(BtCv_downSP))))

#####BtPs significant gene trends
#up-regulated
BtPs_up <- BtPsSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(BtPs_up)

#downregulated
BtPs_down <-  BtPsSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(BtPs_down)

#Exponential Phase TP removal. It happens to be the first column in FullSG_lfc_analysis
BtPsSG_lfc_analysisSP <- BtPsSG_lfc_analysis[,-1]

#up-regulated in stationary phase
BtPs_upSP <- BtPsSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(BtPs_upSP)

#downregulated in stationary phase
BtPs_downSP <-  BtPsSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(BtPs_downSP)

saveRDS(BtPs_up,"diffExp/output/BtPs_up_kallisto.rds")
saveRDS(BtPs_upSP,"diffExp/output/BtPs_upSP_kallisto.rds")
saveRDS(BtPs_down,"diffExp/output/BtPs_down_kallisto.rds")
saveRDS(BtPs_downSP,"diffExp/output/BtPs_downSP_kallisto.rds")

BtPsgeneTrends <- data.frame(trends=c("Up","Down","UpSp","DownSP","VSP"),
  sigGenes=c(nrow(BtPs_up),nrow(BtPs_down),(nrow(BtPs_upSP)-nrow(BtPs_up)),(nrow(BtPs_downSP)-nrow(BtPs_down)),(nrow(BtPsSG)-sum(nrow(BtPs_upSP),nrow(BtPs_downSP)))))

BtPsSigGenes <- c(nrow(BtPs_up),nrow(BtPs_down),(nrow(BtPs_upSP)-nrow(BtPs_up)),(nrow(BtPs_downSP)-nrow(BtPs_down)),(nrow(BtPsSG)-sum(nrow(BtPs_upSP),nrow(BtPs_downSP))))

FinalGeneTrends <- data.frame(condition=rep(c("Full", "BtCv","BtPs"), each=5),
                trend=rep(c("Up", "Down", "UpSP","DownSP", "VSP"),3),
                sigGenes=c(FullSigGenes,BtCvSigGenes,BtPsSigGenes))

write.csv(FinalGeneTrends,"lfc/output/FinalGeneTrends_kallisto.csv")

#####################

#####Mono significant gene trends
#up-regulated
Mono_up <- MonoSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(Mono_up)

Mono_upF <- Mono_up %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,any_vars(. > 1)) %>%
    column_to_rownames('gene')

#downregulated
Mono_down <-  MonoSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(Mono_down)


Mono_downF <- Mono_down %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,any_vars(. < -1)) %>%
    column_to_rownames('gene')


#Variable genes

#Combine up and down only
Mono_oneD <- c(row.names(Mono_up),row.names(Mono_down))

genesToRemove <- which(row.names(MonoSG_lfc_analysis) %in% Mono_oneD)

Mono_var <- MonoSG_lfc_analysis[-c(genesToRemove),]

Mono_varF <- Mono_var %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,any_vars(. > 1 | . < -1)) %>%
    column_to_rownames('gene')

#Observe these variable genes. Check which genes are only >1, <1 or both.
Mono_varFG <- Mono_varF > 1
Mono_varFL <- Mono_varF < -1

apply(Mono_varFG,1,any)
apply(Mono_varFL,1,any)

#In this case, they are only eventually > 1.
#Add these genes to the Mono_up category

Mono_upF <- rbind(Mono_upF,Mono_varF)

saveRDS(Mono_upF,"diffExp/output/Mono_upF_kallisto.rds")
saveRDS(Mono_downF,"diffExp/output/Mono_downF_kallisto.rds")
saveRDS(Mono_varF,"diffExp/output/Mono_varF_kallisto.rds")

library(ggplot2)
FinalGeneTrends <- read.csv(file="lfc/output/FinalGeneTrends_kallisto.csv", sep=",",header=TRUE)

FinalGeneTrends$trend <- factor(FinalGeneTrends$trend,levels = c('Up','Down','UpSP','DownSP','VSP'),ordered = TRUE)
FinalGeneTrends$condition <- factor(FinalGeneTrends$condition,levels = c('Full','BtCv','BtPs'),ordered = TRUE)

p <- ggplot(data=FinalGeneTrends, aes(x=trend,.desc =TRUE, y=sigGenes, fill=condition)) +
geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values=c("black", "brown2","darkgoldenrod2")) +
  theme_minimal()
ggsave("BtSigGeneTrends.tiff",plot=p,device="tiff",width=15, units="cm",dpi=600)


#MA plots
library(ashr)
load("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses/lfc/output/BthailandensisLFC.RData")

#For DESEQ2, plotting MA and getting summarys
#summary(lfc_3mem_int)
#sum(lfc_3mem_int$padj < 0.1, na.rm=TRUE)
#plotMA(lfc_3mem_int,ylim=c(-5,5),alpha = 0.01)

###############COG analysis##########################

Bt_Cog_pre <- read.csv("diffExp/initial_files/Bt_COG_EggNogMapper.csv",header=TRUE,sep=",")

#Remove Categories with less that 10 genes
Bt_Cog <- Bt_Cog_pre[!Bt_Cog_pre$COG %in% names(which(table(Bt_Cog_pre$COG) < 10)), ]

#Analyzing coculture conditions COG categories to monoculture- edit to only include signicant genes
Full_upSP <- readRDS("diffExp/output/Full_upSP_kallisto.rds")
Full_downSP <- readRDS("diffExp/output/Full_downSP_kallisto.rds")
BtCv_upSP <- readRDS("diffExp/output/BtCv_upSP_kallisto.rds")
BtCv_downSP <- readRDS("diffExp/output/BtCv_downSP_kallisto.rds")
BtPs_upSP <- readRDS("diffExp/output/BtPs_upSP_kallisto.rds")
BtPs_downSP <- readRDS("diffExp/output/BtPs_downSP_kallisto.rds")

library(ImpulseDE2)
#Obtain significant genes of interest
impulse_resultsFull <- readRDS("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/diffExp/output/impulseDiffGenesFullvsMono_kallisto.rds")
impulse_resultsCv <- readRDS("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/diffExp/output/impulseDiffGenesBtCvvsMono_kallisto.rds")
impulse_resultsPs <- readRDS("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/diffExp/output/impulseDiffGenesBtPsvsMono_kallisto.rds")

Fullgenes <- impulse_resultsFull@vecDEGenes
BtCvgenes <- impulse_resultsCv@vecDEGenes
BtPsgenes <- impulse_resultsPs@vecDEGenes

#Confirm selected genes are differentially expressed.
table(is.na(match(rownames(Full_upSP),Fullgenes)))
table(is.na(match(rownames(Full_downSP),Fullgenes)))
table(is.na(match(rownames(BtCv_upSP),BtCvgenes)))
table(is.na(match(rownames(BtCv_downSP),BtCvgenes)))
table(is.na(match(rownames(BtPs_upSP),BtPsgenes)))
table(is.na(match(rownames(BtPs_downSP),BtPsgenes)))

#All match differentially expressed genes

#Plot COG categories by up and down regulated genes compared to monoculture

Bt_FullU <- Bt_Cog[which(Bt_Cog$Locus %in% row.names(Full_upSP)),]
Bt_FullD <- Bt_Cog[which(Bt_Cog$Locus %in% row.names(Full_downSP)),]
Bt_CvU <- Bt_Cog[which(Bt_Cog$Locus %in% row.names(BtCv_upSP)),]
Bt_CvD <- Bt_Cog[which(Bt_Cog$Locus %in% row.names(BtCv_downSP)),]
Bt_PsU <- Bt_Cog[which(Bt_Cog$Locus %in% row.names(BtPs_upSP)),]
Bt_PsD <- Bt_Cog[which(Bt_Cog$Locus %in% row.names(BtPs_downSP)),]

Bt_FullU_sum <- as.data.frame(table(Bt_FullU$COG)/table(Bt_Cog$COG)*100)
Bt_FullD_sum <- as.data.frame(table(Bt_FullD$COG)/table(Bt_Cog$COG)*100)
Bt_CvU_sum <- as.data.frame(table(Bt_CvU$COG)/table(Bt_Cog$COG)*100)
Bt_CvD_sum <- as.data.frame(table(Bt_CvD$COG)/table(Bt_Cog$COG)*100)
Bt_PsU_sum <- as.data.frame(table(Bt_PsU$COG)/table(Bt_Cog$COG)*100)
Bt_PsD_sum <- as.data.frame(table(Bt_PsD$COG)/table(Bt_Cog$COG)*100)

#Merge dataframes
library(tidyverse)
Bt_Results <- list(Bt_FullU_sum, Bt_FullD_sum, Bt_CvU_sum,Bt_CvD_sum,Bt_PsU_sum,Bt_PsD_sum) %>% reduce(left_join, by = "Var1")

#Change names
names(Bt_Results) <- c("Category","Full_Upregulation","Full_Downregulation","Cv_Upregulation","Cv_Downregulation","Ps_Upregulation","Ps_Downregulation")

#Reshape
library(reshape2)
library(dplyr)
Bt_ResultsF <- melt(Bt_Results)

#Add condition for facet
Bt_ResultsF$Condition <- c(rep("BtCvPs",nrow(Bt_Results)*2),rep("BtCv",nrow(Bt_Results)*2),rep("BtPs",nrow(Bt_Results)*2))
#Add regulation for facet
Bt_ResultsF$Regulation <- rep(c(rep("Upregulation",nrow(Bt_Results)),rep("Downregulation",nrow(Bt_Results))),3)

#Remove NaN rows
Bt_ResultsF <- Bt_ResultsF[complete.cases(Bt_ResultsF), ]

#Drop levels with no data
Bt_ResultsF$Category <- factor(Bt_ResultsF$Category)

#Let's plot stacked bar
library(ggplot2)

## set the order categories
temp_df_1 <-
  Bt_ResultsF %>%
  arrange(Category)
the_order <- temp_df_1$Category

temp_df_2 <-
  Bt_ResultsF %>%
  arrange(desc(Regulation))
the_order2 <- temp_df_2$Regulation

#Set order of facets
#temp$size_f = factor(temp$size, levels=c('50%','100%','150%','200%'))

#Set order of regulation
#Bt_ResultsF$new = factor(Bt_ResultsF$Regulation, levels=c("Upregulation","Downregulation"), labels=c("Upregulation","Downregulation"))

#Plot
gg <- Bt_ResultsF %>% ggplot(aes(x=Category,weight=value,fill=Regulation)) + geom_bar(position="dodge") +  scale_x_discrete(limits = levels(the_order)) +
  labs(x = "", y = "% genes related to category") + scale_fill_brewer(breaks=c("local", "imported"), palette = "Set1") +
  facet_grid(rows = vars(Condition)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        plot.title = element_text(size=14),
        legend.text = element_text(size=9),
        panel.background = element_rect(fill =  "grey90"))

#Order data
Bt_ResultsF_new <- Bt_ResultsF

Bt_ResultsF_new$Regulation <- factor(Bt_ResultsF_new$Regulation, levels = c("Upregulation","Downregulation"))
Bt_ResultsF_new$Condition = factor(Bt_ResultsF$Condition, levels=c("BtPs","BtCv","BtCvPs"))

#Re-plot
gg <- Bt_ResultsF_new %>% ggplot(aes(x=Category,weight=value,fill=Regulation)) + geom_bar(position="dodge") +  scale_x_discrete(limits = levels(the_order)) +
  labs(x = "", y = "% genes related to category") + scale_fill_brewer(breaks=c("local", "imported"), palette = "Set1") +
  facet_grid(rows = vars(Condition)) +
  scale_fill_manual("legend", values = c("Upregulation" = "blue", "Downregulation" = "red"))
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        plot.title = element_text(size=14),
        legend.text = element_text(size=9),
        panel.background = element_rect(fill =  "grey90"))


gg <- Bt_ResultsF_new %>% ggplot(aes(x=Category,weight=value,fill=Regulation)) + geom_bar(position="dodge") +  scale_x_discrete(limits = levels(the_order)) +
  labs(x = "", y = "% genes related to category") + scale_fill_brewer(breaks=c("local", "imported"), palette = "Set1") +
  facet_grid(rows = vars(Condition)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        plot.title = element_text(size=14),
        legend.text = element_text(size=12),
        strip.text.y = element_text(size = 12),
        panel.background = element_rect(fill =  "grey90"))

gg_final <- gg + scale_fill_manual(values=c("black","grey"))

ggsave("Bt_COG_categories.eps",plot=gg_final,device="eps",width=30, units="cm",dpi=300)


```
