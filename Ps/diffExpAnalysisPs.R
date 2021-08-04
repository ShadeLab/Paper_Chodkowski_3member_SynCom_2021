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
lib <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/initial_files/P-syringae_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/initial_files/Psraw_counts_kallisoNCBI.txt",sep="\t",header=TRUE)

#Remove rRNA, tRNA, and miscRNA
library(rtracklayer)
gff <- readGFF("initial_files/Psyringae_genomic.gff") #if NCBI
CDS <- gff[gff$type=="CDS",]

#For NCBI
library(tidyr)
library(dplyr)
countsS <- counts %>% separate(target_id, c("gene","Locus"),1, sep = "-",remove=TRUE)
countsS <- select(countsS,-c(gene))
countsF <- countsS[which(countsS$Locus %in% CDS$locus_tag),] #For NCBI

#remove non-numerical columns from count matrix
genecounts <- countsF[2:ncol(countsF)]
row.names(genecounts) <- countsF[,1]

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
Extract <- c("3mem","Bt","Cv","mono")
keywords <- str_extract(lib$sampleName, paste(Extract,collapse="|"))

################# Monoculture- case only  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/initial_files/Ps_diffExp_3mem.csv",header=TRUE,sep=",")

#Get mono only
mapping <- mapping[25:48,]

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
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/initial_files/Ps_diffExp_3mem.csv",header=TRUE,sep=",")

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

################# Only PsCv coculture vs control  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/initial_files/Ps_diffExp_2mem_Cv.csv",header=TRUE,sep=",")

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

saveRDS(impulse_resultsCv_trends,"diffExp/output/impulseDiffGenesPsCvvsMono_kallisto.rds")

################# Only PsBt coculture vs control  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/initial_files/Ps_diffExp_2mem_Bt.csv",header=TRUE,sep=",")

keys <- which(keywords %in% c("Bt","mono"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesBt=diffGenesF[,match(align, colnames(diffGenesF))]
colnames(diffGenesBt) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="control", design$Condition)
design$Condition <- sub(pattern="exp",replacement="case", design$Condition)
diffGenesBt <- as.matrix(diffGenesBt)

names(dispersions_deseq) <- row.names(diffGenesBt)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify transients
impulse_resultsBt_trends <- runImpulseDE2(matCountData = diffGenesBt,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsBt_trends,"diffExp/output/impulseDiffGenesPsBtvsMono_kallisto.rds")

########################################

#Create venn diagram of differential gene expression with all coculture conditions compared to monoculture control

diffGenesF <- readRDS("diffExp/initial_files/diffGenesF_kallistoNCBI.rds")

library(ImpulseDE2)

impulse_resultsFull <- readRDS("diffExp/output/impulseDiffGenesFullvsMono_kallisto.rds")
impulse_resultsBt <- readRDS("diffExp/output/impulseDiffGenesPsBtvsMono_kallisto.rds")
impulse_resultsCv <- readRDS("diffExp/output/impulseDiffGenesPsCvvsMono_kallisto.rds")

Fullgenes <- impulse_resultsFull$vecDEGenes
PsBtgenes <- impulse_resultsBt$vecDEGenes
PsCvgenes <- impulse_resultsCv$vecDEGenes

library(VennDiagram)

venn.plot <- venn.diagram(x = list(PsBtCv =Fullgenes,
PsBt=PsBtgenes, PsCv=PsCvgenes),
filename = "ImpulseDE2Summary_comparedToMonoControl.tiff",
    #fill = c("cornflowerblue", "green", "yellow"),
    fill = NULL,
    col = rep("black",3),
    cex = 1.0,
    fontfamily = "serif",
    fontface = "bold",
    cat.cex = 1.5,
    cat.fontfamily = "serif")

#There appear to be 210 differentially expressed genes unique to the full community.
#Are these genes differentially regulated in the Full community compared to the PsBt and PsCv cocultures?

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


genesF <- Setdiff(x = list(Full =Fullgenes),y=list(PsBt=PsBtgenes, PsCv=PsCvgenes))

#Extract these 859 genes out from the gene by sample matrix- 219 from kallisto

diffGenesFOnly <- diffGenesF[which(row.names(diffGenesF) %in% genesF),]

################# Only PsFull vs PsBt  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/initial_files/Ps_diffExp_3mem_v_PsBt.csv",header=TRUE,sep=",")

keys <- which(keywords %in% c("3mem","Bt"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesFO_PsBt=diffGenesFOnly[,match(align, colnames(diffGenesFOnly))]
colnames(diffGenesFO_PsBt) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="control", design$Condition)
design$Condition <- sub(pattern="exp",replacement="case", design$Condition)
diffGenesFO_PsBt <- as.matrix(diffGenesFO_PsBt)

dispersions_deseq <- readRDS("diffExp/initial_files/dispersions_deseq_kallistoNCBI.rds")
dispersions_deseq <- dispersions_deseq[which(row.names(diffGenesF) %in% genesF)]
names(dispersions_deseq) <- row.names(diffGenesFO_PsBt)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify transients
impulse_resultsFO_PsBt_trends <- runImpulseDE2(matCountData = diffGenesFO_PsBt,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsFO_PsBt_trends,"diffExp/output/impulseUniqueDiffGenesFullvsPsBt_kallisto.rds")

################# Only PsFull vs PsCv coculture ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/initial_files/Ps_diffExp_3mem_v_PsCv.csv",header=TRUE,sep=",")

keys <- which(keywords %in% c("3mem","Cv"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesFO_PsCv=diffGenesFOnly[,match(align, colnames(diffGenesFOnly))]
colnames(diffGenesFO_PsCv) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="control", design$Condition)
design$Condition <- sub(pattern="exp",replacement="case", design$Condition)
diffGenesFO_PsCv <- as.matrix(diffGenesFO_PsCv)

dispersions_deseq <- readRDS("diffExp/initial_files/dispersions_deseq_kallistoNCBI.rds")
dispersions_deseq <- dispersions_deseq[which(row.names(diffGenesF) %in% genesF)]
names(dispersions_deseq) <- row.names(diffGenesFO_PsCv)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify transients
impulse_resultsFO_PsCv_trends <- runImpulseDE2(matCountData = diffGenesFO_PsCv,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsFO_PsCv_trends,"diffExp/output/impulseUniqueDiffGenesFullvsPsCv_kallisto.rds")

#################################################################

library(VennDiagram)

impulse_resultsFO_PsBt <- readRDS("diffExp/output/impulseUniqueDiffGenesFullvsPsBt_kallisto.rds")
impulse_resultsFO_PsCv <- readRDS("diffExp/output/impulseUniqueDiffGenesFullvsPsCv_kallisto.rds")

Full_monoCon <- genesF
Full_PsBtCon <- impulse_resultsFO_PsBt$vecDEGenes
Full_PsCvCon <- impulse_resultsFO_PsCv$vecDEGenes

venn.plot <- venn.diagram(x = list(FullvMono =Full_monoCon,
FullvPsBt=Full_PsBtCon, FullvPsCv=Full_PsCvCon),
filename = "ImpulseDE2FulluniqueDiffGenes_kallistoPs.tiff",
col = "transparent",
    fill = c("cornflowerblue", "green", "yellow"),
    #cex = 1.0,
    fontfamily = "serif",
    fontface = "bold",
    #cat.cex = 1.5,
    cat.fontfamily = "serif")

#This venn represents genes that are uniquely differentially regulated in the Full community compared to monoculture. These genes were then tested for significance using BtCv as the control and BtPs as the control.
#So genes in venn diagram show genes that were differentially regulated in the Full community compared to mono and BtCv,mono and BtPs, mono only, and all other conditions.

final_sigGeneFullOnly <- intersect(Full_monoCon,intersect(Full_PsBtCon, Full_PsCvCon))
saveRDS(final_sigGeneFullOnly,"diffExp/output/sigGenesFullOnly_Ps_kallisto.rds")

#Parse out other genes of interest
#overlap <- calculate.overlap(list("FullvMono" =Full_monoCon,
#"FullvPsBt"=Full_PsBtCon, "FullvPsCv"=Full_PsCvCon))

#final_sigGenesBtCv_Mono <- overlap$a2
#final_sigGenesBtPs_Mono <- overlap$a4

#Write these genes to csv
write.csv(final_sigGeneFullOnly, "SigGenesUniqueToPs3mem.csv")

######################Significant gene trends##########################

library(ImpulseDE2)

impulse_resultsFull <- readRDS("diffExp/output/impulseDiffGenesFullvsMono_kallisto.rds")
impulse_resultsBt <- readRDS("diffExp/output/impulseDiffGenesPsBtvsMono_kallisto.rds")
impulse_resultsCv <- readRDS("diffExp/output/impulseDiffGenesPsCvvsMono_kallisto.rds")
impulse_resultsMono <- readRDS("diffExp/output/impulseDiffGenesMonoCaseOnly_kallisto.rds")

load("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Pseudomonas/final_analyses_Rev/lfc/output/PsyringaeLFC_kallisto.RData")

lfc_analysisFull <- cbind(lfc_3memInt_comparedtoMonoInt.fc,lfc_3mem25_comparedtoMono25.fc,lfc_3mem30_comparedtoMono30.fc,
                          lfc_3mem35_comparedtoMono35.fc,lfc_3mem40_comparedtoMono40.fc,lfc_3mem45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisFull) <- row.names(lfc_3memInt_comparedtoMonoInt)

lfc_analysisPsBt <- cbind(lfc_PsBtInt_comparedtoMonoInt.fc,lfc_PsBt25_comparedtoMono25.fc,lfc_PsBt30_comparedtoMono30.fc,
                          lfc_PsBt35_comparedtoMono35.fc,lfc_PsBt40_comparedtoMono40.fc,lfc_PsBt45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisPsBt) <- row.names(lfc_PsBtInt_comparedtoMonoInt)

lfc_analysisPsCv <- cbind(lfc_PsCvInt_comparedtoMonoInt.fc,lfc_PsCv25_comparedtoMono25.fc,lfc_PsCv30_comparedtoMono30.fc,
                          lfc_PsCv35_comparedtoMono35.fc,lfc_PsCv40_comparedtoMono40.fc,lfc_PsCv45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisPsCv) <- row.names(lfc_PsCvInt_comparedtoMonoInt)

lfc_analysisMono <- cbind(lfc_mono_25withinComparisionToTime0.fc,lfc_mono_30withinComparisionToTime0.fc,lfc_mono_35withinComparisionToTime0.fc,
                          lfc_mono_40withinComparisionToTime0.fc,lfc_mono_45withinComparisionToTime0.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisMono) <- row.names(lfc_mono25_withinComparision_ToTime0)


FullSG <- lfc_analysisFull[match(impulse_resultsFull$vecDEGenes,row.names(lfc_analysisFull)),]
PsBtSG <- lfc_analysisPsBt[match(impulse_resultsBt$vecDEGenes,row.names(lfc_analysisPsBt)),]
PsCvSG <- lfc_analysisPsCv[match(impulse_resultsCv$vecDEGenes,row.names(lfc_analysisPsCv)),]
MonoSG <- lfc_analysisMono[match(impulse_resultsMono$vecDEGenes,row.names(lfc_analysisMono)),]

FullSG_lfc_analysis <- as.data.frame(FullSG)
PsBtSG_lfc_analysis <- as.data.frame(PsBtSG)
PsCvSG_lfc_analysis <- as.data.frame(PsCvSG)
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

#####PsBt significant gene trends
#up-regulated
PsBt_up <- PsBtSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(PsBt_up)

#downregulated
PsBt_down <-  PsBtSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(PsBt_down)

#Exponential Phase TP removal. It happens to be the first column in FullSG_lfc_analysis
PsBtSG_lfc_analysisSP <- PsBtSG_lfc_analysis[,-1]

#up-regulated in stationary phase
PsBt_upSP <- PsBtSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(PsBt_upSP)

#downregulated in stationary phase
PsBt_downSP <-  PsBtSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(PsBt_downSP)

saveRDS(PsBt_up,"diffExp/output/PsBt_up_kallisto.rds")
saveRDS(PsBt_upSP,"diffExp/output/PsBt_upSP_kallisto.rds")
saveRDS(PsBt_down,"diffExp/output/PsBt_down_kallisto.rds")
saveRDS(PsBt_downSP,"diffExp/output/PsBt_downSP_kallisto.rds")

PsBtgeneTrends <- data.frame(trends=c("Up","Down","UpSp","DownSP","VSP"),
  sigGenes=c(nrow(PsBt_up),nrow(PsBt_down),(nrow(PsBt_upSP)-nrow(PsBt_up)),(nrow(PsBt_downSP)-nrow(PsBt_down)),(nrow(PsBtSG)-sum(nrow(PsBt_upSP),nrow(PsBt_downSP)))))

PsBtSigGenes <- c(nrow(PsBt_up),nrow(PsBt_down),(nrow(PsBt_upSP)-nrow(PsBt_up)),(nrow(PsBt_downSP)-nrow(PsBt_down)),(nrow(PsBtSG)-sum(nrow(PsBt_upSP),nrow(PsBt_downSP))))

#####PsCv significant gene trends
#up-regulated
PsCv_up <- PsCvSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(PsCv_up)

#downregulated
PsCv_down <-  PsCvSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(PsCv_down)

#Exponential Phase TP removal. It happens to be the first column in FullSG_lfc_analysis
PsCvSG_lfc_analysisSP <- PsCvSG_lfc_analysis[,-1]

#up-regulated in stationary phase
PsCv_upSP <- PsCvSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(PsCv_upSP)

#downregulated in stationary phase
PsCv_downSP <-  PsCvSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(PsCv_downSP)

saveRDS(PsCv_up,"diffExp/output/PsCv_up_kallisto.rds")
saveRDS(PsCv_upSP,"diffExp/output/PsCv_upSP_kallisto.rds")
saveRDS(PsCv_down,"diffExp/output/PsCv_down_kallisto.rds")
saveRDS(PsCv_downSP,"diffExp/output/PsCv_downSP_kallisto.rds")

PsCvgeneTrends <- data.frame(trends=c("Up","Down","UpSp","DownSP","VSP"),
  sigGenes=c(nrow(PsCv_up),nrow(PsCv_down),(nrow(PsCv_upSP)-nrow(PsCv_up)),(nrow(PsCv_downSP)-nrow(PsCv_down)),(nrow(PsCvSG)-sum(nrow(PsCv_upSP),nrow(PsCv_downSP)))))

PsCvSigGenes <- c(nrow(PsCv_up),nrow(PsCv_down),(nrow(PsCv_upSP)-nrow(PsCv_up)),(nrow(PsCv_downSP)-nrow(PsCv_down)),(nrow(PsCvSG)-sum(nrow(PsCv_upSP),nrow(PsCv_downSP))))

FinalGeneTrends <- data.frame(condition=rep(c("Full", "PsBt","PsCv"), each=5),
                trend=rep(c("Up", "Down", "UpSP","DownSP", "VSP"),3),
                sigGenes=c(FullSigGenes,PsBtSigGenes,PsCvSigGenes))

write.csv(FinalGeneTrends,"lfc/output/FinalGeneTrends_kallisto.csv")

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

muo <- apply(Mono_varFG,1,any)
mdo <- apply(Mono_varFL,1,any)

#Which genes are >1 and <-1 during SP?
muo == mdo

#No variable genes

saveRDS(Mono_upF,"diffExp/output/Mono_upF_kallisto.rds")
saveRDS(Mono_downF,"diffExp/output/Mono_downF_kallisto.rds")
#saveRDS(Mono_varF,"diffExp/output/Mono_varF_kallisto.rds")

library(ggplot2)
FinalGeneTrends <- read.csv(file="lfc/output/FinalGeneTrends_kallisto.csv", sep=",",header=TRUE)

FinalGeneTrends$trend <- factor(FinalGeneTrends$trend,levels = c('Up','Down','UpSP','DownSP','VSP'),ordered = TRUE)
FinalGeneTrends$condition <- factor(FinalGeneTrends$condition,levels = c('Full','PsBt','PsCv'),ordered = TRUE)

p <- ggplot(data=FinalGeneTrends, aes(x=trend,.desc =TRUE, y=sigGenes, fill=condition)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_manual(values=c("black", "darkgoldenrod2","burlywood4")) +
  theme_minimal()
ggsave("PsSigGeneTrends.tiff",plot=p,device="tiff",width=15, units="cm",dpi=600)


#For DESEQ2, plotting MA and getting summarys
#summary(Full_int)
#sum(Full_int$padj < 0.1, na.rm=TRUE)
#plotMA(Full_int,ylim=c(-5,5))

###############COG analysis##################

Ps_Cog_pre <- read.csv("diffExp/initial_files/Ps_COG_EggNogMapper.csv",header=TRUE,sep=",")

#Remove Categories with less that 10 genes
Ps_Cog <- Ps_Cog_pre[!Ps_Cog_pre$COG %in% names(which(table(Ps_Cog_pre$COG) < 10)), ]

#Analyzing coculture conditions COG categories to monoculture- edit to only include significant genes
Full_upSP <- readRDS("diffExp/output/Full_upSP_kallisto.rds")
Full_downSP <- readRDS("diffExp/output/Full_downSP_kallisto.rds")
PsBt_upSP <- readRDS("diffExp/output/PsBt_upSP_kallisto.rds")
PsBt_downSP <- readRDS("diffExp/output/PsBt_downSP_kallisto.rds")
PsCv_upSP <- readRDS("diffExp/output/PsCv_upSP_kallisto.rds")
PsCv_downSP <- readRDS("diffExp/output/PsCv_downSP_kallisto.rds")

library(ImpulseDE2)
#Obtain significant genes of interest
impulse_resultsFull <- readRDS("diffExp/output/impulseDiffGenesFullvsMono_kallisto.rds")
impulse_resultsBt <- readRDS("diffExp/output/impulseDiffGenesPsBtvsMono_kallisto.rds")
impulse_resultsCv <- readRDS("diffExp/output/impulseDiffGenesPsCvvsMono_kallisto.rds")

Fullgenes <- impulse_resultsFull@vecDEGenes
PsBtgenes <- impulse_resultsBt@vecDEGenes
PsCvgenes <- impulse_resultsCv@vecDEGenes

#Confirm selected genes are differentially expressed.
table(is.na(match(rownames(Full_upSP),Fullgenes)))
table(is.na(match(rownames(Full_downSP),Fullgenes)))
table(is.na(match(rownames(PsBt_upSP),PsBtgenes)))
table(is.na(match(rownames(PsBt_downSP),PsBtgenes)))
table(is.na(match(rownames(PsCv_upSP),PsCvgenes)))
table(is.na(match(rownames(PsCv_downSP),PsCvgenes)))

#All match differentially expressed genes

#Plot COG categories by up and down regulated genes compared to monoculture

Ps_FullU <- Ps_Cog[which(Ps_Cog$Locus %in% row.names(Full_upSP)),]
Ps_FullD <- Ps_Cog[which(Ps_Cog$Locus %in% row.names(Full_downSP)),]
Ps_BtU <- Ps_Cog[which(Ps_Cog$Locus %in% row.names(PsBt_upSP)),]
Ps_BtD <- Ps_Cog[which(Ps_Cog$Locus %in% row.names(PsBt_downSP)),]
Ps_CvU <- Ps_Cog[which(Ps_Cog$Locus %in% row.names(PsCv_upSP)),]
Ps_CvD <- Ps_Cog[which(Ps_Cog$Locus %in% row.names(PsCv_downSP)),]

Ps_FullU_sum <- as.data.frame(table(Ps_FullU$COG)/table(Ps_Cog$COG)*100)
Ps_FullD_sum <- as.data.frame(table(Ps_FullD$COG)/table(Ps_Cog$COG)*100)
Ps_BtU_sum <- as.data.frame(table(Ps_BtU$COG)/table(Ps_Cog$COG)*100)
Ps_BtD_sum <- as.data.frame(table(Ps_BtD$COG)/table(Ps_Cog$COG)*100)
Ps_CvU_sum <- as.data.frame(table(Ps_CvU$COG)/table(Ps_Cog$COG)*100)
Ps_CvD_sum <- as.data.frame(table(Ps_CvD$COG)/table(Ps_Cog$COG)*100)

#Merge dataframes
library(tidyverse)
Ps_Results <- list(Ps_FullU_sum, Ps_FullD_sum, Ps_BtU_sum,Ps_BtD_sum,Ps_CvU_sum,Ps_CvD_sum) %>% reduce(left_join, by = "Var1")

#Change names
names(Ps_Results) <- c("Category","Full_Upregulation","Full_Downregulation","Bt_Upregulation","Bt_Downregulation","Cv_Upregulation","Cv_Downregulation")

#Reshape
library(reshape2)
library(dplyr)
Ps_ResultsF <- melt(Ps_Results)

#Add condition for facet
Ps_ResultsF$Condition <- c(rep("PsBtCv",nrow(Ps_Results)*2),rep("PsBt",nrow(Ps_Results)*2),rep("PsCv",nrow(Ps_Results)*2))
#Add regulation for facet
Ps_ResultsF$Regulation <- rep(c(rep("Upregulation",nrow(Ps_Results)),rep("Downregulation",nrow(Ps_Results))),3)

#Remove NaN rows
Ps_ResultsF <- Ps_ResultsF[complete.cases(Ps_ResultsF), ]

#Drop levels with no data
Ps_ResultsF$Category <- factor(Ps_ResultsF$Category)

#Let's plot stacked bar
library(ggplot2)

## set the order categories
temp_df_1 <-
  Ps_ResultsF %>%
  arrange(Category)
the_order <- temp_df_1$Category

#Plot
gg <- Ps_ResultsF %>% ggplot(aes(x=Category,weight=value,fill=Regulation)) + geom_bar(position="dodge") +  scale_x_discrete(limits = levels(the_order)) +
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
Ps_ResultsF_new <- Ps_ResultsF

Ps_ResultsF_new$Regulation <- factor(Ps_ResultsF_new$Regulation, levels = c("Upregulation","Downregulation"))
Ps_ResultsF_new$Condition = factor(Ps_ResultsF$Condition, levels=c("PsCv","PsBt","PsBtCv"))

#Re-plot
gg <- Ps_ResultsF_new %>% ggplot(aes(x=Category,weight=value,fill=Regulation)) + geom_bar(position="dodge") +  scale_x_discrete(limits = levels(the_order)) +
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

ggsave("Ps_COG_categories.eps",plot=gg_final,device="eps",width=30, units="cm",dpi=300)





```
