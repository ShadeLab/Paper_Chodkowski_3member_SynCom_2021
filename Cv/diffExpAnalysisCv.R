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
lib <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/initial_files/Cviolaceum_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/initial_files/Cvraw_counts_kallisoNCBI.txt",sep="\t",header=TRUE)

#Remove rRNA, tRNA, and miscRNA
library(rtracklayer)
gff <- readGFF("initial_files/Cviolaceum_genomic.gff") #for NCBI
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
Extract <- c("3mem","Bt","Ps","mono")
keywords <- str_extract(lib$sampleName, paste(Extract,collapse="|"))

################# Monoculture- case only  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/initial_files/Cv_diffExp_3mem.csv",header=TRUE,sep=",")

#Get mono only
mapping <- mapping[24:46,]

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
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/initial_files/Cv_diffExp_3mem.csv",header=TRUE,sep=",")

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

#Identify DEGs
#NOTE: matCountData should be raw counts, not normalized orelse you'll be normalizing normalized data.
impulse_results_trends <- runImpulseDE2(matCountData = diffGenesFull,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_results_trends,"diffExp/output/impulseDiffGenesFullvsMono_kallisto.rds")

################# Only CvBt coculture vs control  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/initial_files/Cv_diffExp_2mem_Bt.csv",header=TRUE,sep=",")

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

#Identify DEGs
impulse_resultsCv_trends <- runImpulseDE2(matCountData = diffGenesBt,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsCv_trends,"diffExp/output/impulseDiffGenesCvBtvsMono_kallisto.rds")

################# Only CvPs coculture vs control  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/initial_files/Cv_diffExp_2mem_Ps.csv",header=TRUE,sep=",")

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

#Identify DEGs
impulse_resultsPs_trends <- runImpulseDE2(matCountData = diffGenesPs,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsPs_trends,"diffExp/output/impulseDiffGenesCvPsvsMono_kallisto.rds")

########################################
#Now that we have significant genes determined by ImpulseDE2 model, let's filter DEGs by DEGs that have at least one TP with LFC >= 1 or <= -1
library(ImpulseDE2)

impulse_resultsFull <- readRDS("diffExp/output/impulseDiffGenesFullvsMono_kallisto.rds")
impulse_resultsBt <- readRDS("diffExp/output/impulseDiffGenesCvBtvsMono_kallisto.rds")
impulse_resultsPs <- readRDS("diffExp/output/impulseDiffGenesCvPsvsMono_kallisto.rds")

load("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/lfc/output/CviolaceumLFC_kallisto.RData")

lfc_analysisFull <- cbind(lfc_3memInt_comparedtoMonoInt.fc,lfc_3mem25_comparedtoMono25.fc,lfc_3mem30_comparedtoMono30.fc,
                          lfc_3mem35_comparedtoMono35.fc,lfc_3mem40_comparedtoMono40.fc,lfc_3mem45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisFull) <- row.names(lfc_3memInt_comparedtoMonoInt)

lfc_analysisCvBt <- cbind(lfc_CvBtInt_comparedtoMonoInt.fc,lfc_CvBt25_comparedtoMono25.fc,lfc_CvBt30_comparedtoMono30.fc,
                          lfc_CvBt35_comparedtoMono35.fc,lfc_CvBt40_comparedtoMono40.fc,lfc_CvBt45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisCvBt) <- row.names(lfc_CvBtInt_comparedtoMonoInt)

lfc_analysisCvPs <- cbind(lfc_CvPsInt_comparedtoMonoInt.fc,lfc_CvPs25_comparedtoMono25.fc,lfc_CvPs30_comparedtoMono30.fc,
                          lfc_CvPs35_comparedtoMono35.fc,lfc_CvPs40_comparedtoMono40.fc,lfc_CvPs45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisCvPs) <- row.names(lfc_CvPsInt_comparedtoMonoInt)

lfc_analysisMono <- cbind(lfc_mono_25withinComparisionToTime0.fc,lfc_mono_30withinComparisionToTime0.fc,lfc_mono_35withinComparisionToTime0.fc,
                          lfc_mono_40withinComparisionToTime0.fc,lfc_mono_45withinComparisionToTime0.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisMono) <- row.names(lfc_mono25_withinComparision_ToTime0)

FullSG <- lfc_analysisFull[match(impulse_resultsFull$vecDEGenes,row.names(lfc_analysisFull)),]
CvBtSG <- lfc_analysisCvBt[match(impulse_resultsBt$vecDEGenes,row.names(lfc_analysisCvBt)),]
CvPsSG <- lfc_analysisCvPs[match(impulse_resultsPs$vecDEGenes,row.names(lfc_analysisCvPs)),]

FullSG_lfc_analysis <- as.data.frame(FullSG)
CvBtSG_lfc_analysis <- as.data.frame(CvBtSG)
CvPsSG_lfc_analysis <- as.data.frame(CvPsSG)

library(dplyr)
library(tibble)

#####Full community LFC thresholds
#up-regulated
Full_up <- FullSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,any_vars(. >= 1 )) %>%
    column_to_rownames('gene')

nrow(Full_up)

Full_down <- FullSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,any_vars(. <= -1 )) %>%
    column_to_rownames('gene')

nrow(Full_down)

#Find duplicate genes from both up/downregulated
duprows <- rownames(Full_down) %in% rownames(Full_up)

#Combine lists
Full_LFC_thresh <- rbind(Full_up,Full_down[!duprows,])

#Save these as that final DEGs with a LFC threshold
saveRDS(rownames(Full_LFC_thresh),"diffExp/output/impulseResultsFull_FDR-LFC_kallisto.rds")

#####CvBt LFC threshold
#up-regulated
CvBt_up <- CvBtSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,any_vars(. >= 1 )) %>%
    column_to_rownames('gene')

nrow(CvBt_up)

#downregulated
CvBt_down <-  CvBtSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,any_vars(. <= -1 )) %>%
    column_to_rownames('gene')

nrow(CvBt_down)

#Find duplicate genes from both up/downregulated
duprows <- rownames(CvBt_down) %in% rownames(CvBt_up)

#Combine lists
CvBt_LFC_thresh <- rbind(CvBt_up,CvBt_down[!duprows,])

#Save these as that final DEGs with a LFC threshold
saveRDS(rownames(CvBt_LFC_thresh),"diffExp/output/impulseResultsCvBt_FDR-LFC_kallisto.rds")

#####CvPs LFC threshold
#up-regulated
CvPs_up <- CvPsSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,any_vars(. >= 1 )) %>%
    column_to_rownames('gene')

nrow(CvPs_up)

#downregulated
CvPs_down <-  CvPsSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,any_vars(. <= -1 )) %>%
    column_to_rownames('gene')

nrow(CvPs_down)

#Find duplicate genes from both up/downregulated
duprows <- rownames(CvPs_down) %in% rownames(CvPs_up)

#Combine lists
CvPs_LFC_thresh <- rbind(CvPs_up,CvPs_down[!duprows,])

#Save these as that final DEGs with a LFC threshold
saveRDS(rownames(CvPs_LFC_thresh),"diffExp/output/impulseResultsCvPs_FDR-LFC_kallisto.rds")

#############################
#Create venn diagram of differential gene expression with all coculture conditions compared to monoculture control

#These will be DEGS with a LFC threshold of at least 1 or -1 at 1 TP

diffGenesF <- readRDS("diffExp/initial_files/diffGenesF_kallistoNCBI.rds")

library(ImpulseDE2)

Fullgenes <- readRDS("diffExp/output/impulseResultsFull_FDR-LFC_kallisto.rds")
CvBtgenes <- readRDS("diffExp/output/impulseResultsCvBt_FDR-LFC_kallisto.rds")
CvPsgenes <- readRDS("diffExp/output/impulseResultsCvPs_FDR-LFC_kallisto.rds")

#Fullgenes <- impulse_resultsFull$vecDEGenes
#CvBtgenes <- impulse_resultsBt$vecDEGenes
#CvPsgenes <- impulse_resultsPs$vecDEGenes

library(VennDiagram)

venn.plot <- venn.diagram(x = list(CvPsBt =Fullgenes,
CvBt=CvBtgenes, CvPs=CvPsgenes),
filename = "ImpulseDE2Summary_comparedToMonoControl.tiff",
    #fill = c("cornflowerblue", "green", "yellow"),
    fill = NULL,
    col = rep("black",3),
    cex = 1.0,
    fontfamily = "serif",
    fontface = "bold",
    cat.cex = 1.5,
    cat.fontfamily = "serif")

#There appear to be 189 differentially expressed genes unique to the full community.
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


genesF <- Setdiff(x = list(Full =Fullgenes),y=list(CvBt=CvBtgenes, CvPs=CvPsgenes))

#Extract these 139 genes out from the gene by sample matrix- 167 from kallisto

diffGenesFOnly <- diffGenesF[which(row.names(diffGenesF) %in% genesF),]

################# Only CvFull vs CvBt  ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/initial_files/Cv_diffExp_3mem_v_CvBt.csv",header=TRUE,sep=",")

keys <- which(keywords %in% c("3mem","Bt"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesFO_CvBt=diffGenesFOnly[,match(align, colnames(diffGenesFOnly))]
colnames(diffGenesFO_CvBt) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="control", design$Condition)
design$Condition <- sub(pattern="exp",replacement="case", design$Condition)
diffGenesFO_CvBt <- as.matrix(diffGenesFO_CvBt)

dispersions_deseq <- readRDS("diffExp/initial_files/dispersions_deseq_kallistoNCBI.rds")
dispersions_deseq <- dispersions_deseq[which(row.names(diffGenesF) %in% genesF)]
names(dispersions_deseq) <- row.names(diffGenesFO_CvBt)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify transients
impulse_resultsFO_CvBt_trends <- runImpulseDE2(matCountData = diffGenesFO_CvBt,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsFO_CvBt_trends,"diffExp/output/impulseUniqueDiffGenesFullvsCvBt_kallisto.rds")

################# Only CvFull vs CvPs coculture ########################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/initial_files/Cv_diffExp_3mem_v_CvPs.csv",header=TRUE,sep=",")

keys <- which(keywords %in% c("3mem","Ps"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesFO_CvPs=diffGenesFOnly[,match(align, colnames(diffGenesFOnly))]
colnames(diffGenesFO_CvPs) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="control", design$Condition)
design$Condition <- sub(pattern="exp",replacement="case", design$Condition)
diffGenesFO_CvPs <- as.matrix(diffGenesFO_CvPs)

dispersions_deseq <- readRDS("diffExp/initial_files/dispersions_deseq_kallistoNCBI.rds")
dispersions_deseq <- dispersions_deseq[which(row.names(diffGenesF) %in% genesF)]
names(dispersions_deseq) <- row.names(diffGenesFO_CvPs)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify transients
impulse_resultsFO_CvPs_trends <- runImpulseDE2(matCountData = diffGenesFO_CvPs,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsFO_CvPs_trends,"diffExp/output/impulseUniqueDiffGenesFullvsCvPs_kallisto.rds")

#################################################################
library(VennDiagram)

impulse_resultsFO_CvBt <- readRDS("diffExp/output/impulseUniqueDiffGenesFullvsCvBt_kallisto.rds")
impulse_resultsFO_CvPs <- readRDS("diffExp/output/impulseUniqueDiffGenesFullvsCvPs_kallisto.rds")

Full_monoCon <- genesF
Full_CvBtCon <- impulse_resultsFO_CvBt$vecDEGenes
Full_CvPsCon <- impulse_resultsFO_CvPs$vecDEGenes

venn.plot <- venn.diagram(x = list(FullvMono =Full_monoCon,
FullvCvBt=Full_CvBtCon, FullvCvPs=Full_CvPsCon),
filename = "ImpulseDE2FulluniqueDiffGenes_kallistoCv.tiff",
col = "transparent",
    fill = c("cornflowerblue", "green", "yellow"),
    #cex = 1.0,
    fontfamily = "serif",
    fontface = "bold",
    #cat.cex = 1.5,
    cat.fontfamily = "serif")

#This venn represents genes that are uniquely differentially regulated in the Full community compared to monoculture. These genes were then tested for significance using CvBt as the control and CvPs as the control.
#So genes in venn diagram show genes that were differentially regulated in the Full community compared to mono and CvBt,mono and CvPs, mono only, and all other conditions.

final_sigGeneFullOnly <- intersect(Full_monoCon,intersect(Full_CvBtCon, Full_CvPsCon))
saveRDS(final_sigGeneFullOnly,"diffExp/output/sigGenesFullOnly_Cv_kallisto.rds")

#Parse out other genes of interest
#overlap <- calculate.overlap(list("FullvMono" =Full_monoCon,
#"FullvCvBt"=Full_CvBtCon, "FullvCvPs"=Full_CvPsCon))

#final_sigGenesBtCv_Mono <- overlap$a2
#final_sigGenesBtPs_Mono <- overlap$a4

#Write these genes to csv
write.csv(final_sigGeneFullOnly, "SigGenesUniqueToCv3mem.csv")

######################Significant gene trends##########################

library(ImpulseDE2)

impulse_resultsFull <- readRDS("diffExp/output/impulseResultsFull_FDR-LFC_kallisto.rds")
impulse_resultsBt <- readRDS("diffExp/output/impulseResultsCvBt_FDR-LFC_kallisto.rds")
impulse_resultsPs <- readRDS("diffExp/output/impulseResultsCvPs_FDR-LFC_kallisto.rds")

load("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/lfc/output/CviolaceumLFC_kallisto.RData")

lfc_analysisFull <- cbind(lfc_3memInt_comparedtoMonoInt.fc,lfc_3mem25_comparedtoMono25.fc,lfc_3mem30_comparedtoMono30.fc,
                          lfc_3mem35_comparedtoMono35.fc,lfc_3mem40_comparedtoMono40.fc,lfc_3mem45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisFull) <- row.names(lfc_3memInt_comparedtoMonoInt)

lfc_analysisCvBt <- cbind(lfc_CvBtInt_comparedtoMonoInt.fc,lfc_CvBt25_comparedtoMono25.fc,lfc_CvBt30_comparedtoMono30.fc,
                          lfc_CvBt35_comparedtoMono35.fc,lfc_CvBt40_comparedtoMono40.fc,lfc_CvBt45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisCvBt) <- row.names(lfc_CvBtInt_comparedtoMonoInt)

lfc_analysisCvPs <- cbind(lfc_CvPsInt_comparedtoMonoInt.fc,lfc_CvPs25_comparedtoMono25.fc,lfc_CvPs30_comparedtoMono30.fc,
                          lfc_CvPs35_comparedtoMono35.fc,lfc_CvPs40_comparedtoMono40.fc,lfc_CvPs45_comparedtoMono45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisCvPs) <- row.names(lfc_CvPsInt_comparedtoMonoInt)

lfc_analysisMono <- cbind(lfc_mono_25withinComparisionToTime0.fc,lfc_mono_30withinComparisionToTime0.fc,lfc_mono_35withinComparisionToTime0.fc,
                          lfc_mono_40withinComparisionToTime0.fc,lfc_mono_45withinComparisionToTime0.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisMono) <- row.names(lfc_mono25_withinComparision_ToTime0)

FullSG <- lfc_analysisFull[match(impulse_resultsFull,row.names(lfc_analysisFull)),]
CvBtSG <- lfc_analysisCvBt[match(impulse_resultsBt,row.names(lfc_analysisCvBt)),]
CvPsSG <- lfc_analysisCvPs[match(impulse_resultsPs,row.names(lfc_analysisCvPs)),]

FullSG_lfc_analysis <- as.data.frame(FullSG)
CvBtSG_lfc_analysis <- as.data.frame(CvBtSG)
CvPsSG_lfc_analysis <- as.data.frame(CvPsSG)

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

#####CvBt significant gene trends
#up-regulated
CvBt_up <- CvBtSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(CvBt_up)

#downregulated
CvBt_down <-  CvBtSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(CvBt_down)

#Exponential Phase TP removal. It happens to be the first column in FullSG_lfc_analysis
CvBtSG_lfc_analysisSP <- CvBtSG_lfc_analysis[,-1]

#up-regulated in stationary phase
CvBt_upSP <- CvBtSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(CvBt_upSP)

#downregulated in stationary phase
CvBt_downSP <-  CvBtSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(CvBt_downSP)

saveRDS(CvBt_up,"diffExp/output/CvBt_up_kallisto.rds")
saveRDS(CvBt_upSP,"diffExp/output/CvBt_upSP_kallisto.rds")
saveRDS(CvBt_down,"diffExp/output/CvBt_down_kallisto.rds")
saveRDS(CvBt_downSP,"diffExp/output/CvBt_downSP_kallisto.rds")

CvBtgeneTrends <- data.frame(trends=c("Up","Down","UpSp","DownSP","VSP"),
  sigGenes=c(nrow(CvBt_up),nrow(CvBt_down),(nrow(CvBt_upSP)-nrow(CvBt_up)),(nrow(CvBt_downSP)-nrow(CvBt_down)),(nrow(CvBtSG)-sum(nrow(CvBt_upSP),nrow(CvBt_downSP)))))

CvBtSigGenes <- c(nrow(CvBt_up),nrow(CvBt_down),(nrow(CvBt_upSP)-nrow(CvBt_up)),(nrow(CvBt_downSP)-nrow(CvBt_down)),(nrow(CvBtSG)-sum(nrow(CvBt_upSP),nrow(CvBt_downSP))))

#####CvPs significant gene trends
#up-regulated
CvPs_up <- CvPsSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(CvPs_up)

#downregulated
CvPs_down <-  CvPsSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(CvPs_down)

#Exponential Phase TP removal. It happens to be the first column in FullSG_lfc_analysis
CvPsSG_lfc_analysisSP <- CvPsSG_lfc_analysis[,-1]

#up-regulated in stationary phase
CvPs_upSP <- CvPsSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(CvPs_upSP)

#downregulated in stationary phase
CvPs_downSP <-  CvPsSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(CvPs_downSP)

saveRDS(CvPs_up,"diffExp/output/CvPs_up_kallisto.rds")
saveRDS(CvPs_upSP,"diffExp/output/CvPs_upSP_kallisto.rds")
saveRDS(CvPs_down,"diffExp/output/CvPs_down_kallisto.rds")
saveRDS(CvPs_downSP,"diffExp/output/CvPs_downSP_kallisto.rds")

CvPsgeneTrends <- data.frame(trends=c("Up","Down","UpSp","DownSP","VSP"),
  sigGenes=c(nrow(CvPs_up),nrow(CvPs_down),(nrow(CvPs_upSP)-nrow(CvPs_up)),(nrow(CvPs_downSP)-nrow(CvPs_down)),(nrow(CvPsSG)-sum(nrow(CvPs_upSP),nrow(CvPs_downSP)))))

CvPsSigGenes <- c(nrow(CvPs_up),nrow(CvPs_down),(nrow(CvPs_upSP)-nrow(CvPs_up)),(nrow(CvPs_downSP)-nrow(CvPs_down)),(nrow(CvPsSG)-sum(nrow(CvPs_upSP),nrow(CvPs_downSP))))

FinalGeneTrends <- data.frame(condition=rep(c("Full", "CvBt","CvPs"), each=5),
                trend=rep(c("Up", "Down", "UpSP","DownSP", "VSP"),3),
                sigGenes=c(FullSigGenes,CvBtSigGenes,CvPsSigGenes))

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

#Two for now, ignore but fix later!
Mono_varFGF <- Mono_varF[muo,]
Mono_varFLF <- Mono_varF[mdo,]

Mono_upF <- rbind(Mono_upF,Mono_varFGF)
Mono_downF <- rbind(Mono_downF,Mono_varFLF)

saveRDS(Mono_upF,"diffExp/output/Mono_upF_kallisto.rds")
saveRDS(Mono_downF,"diffExp/output/Mono_downF_kallisto.rds")
saveRDS(Mono_varF,"diffExp/output/Mono_varF_kallisto.rds")

library(ggplot2)
FinalGeneTrends <- read.csv(file="lfc/output/FinalGeneTrends_kallisto.csv", sep=",",header=TRUE)

FinalGeneTrends$trend <- factor(FinalGeneTrends$trend,levels = c('Up','Down','UpSP','DownSP','VSP'),ordered = TRUE)
FinalGeneTrends$condition <- factor(FinalGeneTrends$condition,levels = c('Full','CvBt','CvPs'),ordered = TRUE)

p <- ggplot(data=FinalGeneTrends, aes(x=trend,.desc =TRUE, y=sigGenes, fill=condition)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_manual(values=c("black", "brown2","burlywood4")) +
  theme_minimal()
ggsave("CvSigGeneTrends.tiff",plot=p,device="tiff",width=15, units="cm",dpi=600)


#For DESEQ2, plotting MA and getting summarys
#summary(Full_int)
#sum(Full_int$padj < 0.1, na.rm=TRUE)
#plotMA(Full_int,ylim=c(-5,5))

##############COG analysis#########################

Cv_Cog_pre <- read.csv("diffExp/initial_files/Cv_COG_EggNogMapper.csv",header=TRUE,sep=",")

#Remove Categories with less that 10 genes
Cv_Cog <- Cv_Cog_pre[!Cv_Cog_pre$COG %in% names(which(table(Cv_Cog_pre$COG) < 10)), ]

#Analyzing coculture conditions COG categories to monoculture- edit to only include signicant genes
Full_upSP <- readRDS("diffExp/output/Full_upSP_kallisto.rds")
Full_downSP <- readRDS("diffExp/output/Full_downSP_kallisto.rds")
CvBt_upSP <- readRDS("diffExp/output/CvBt_upSP_kallisto.rds")
CvBt_downSP <- readRDS("diffExp/output/CvBt_downSP_kallisto.rds")
CvPs_upSP <- readRDS("diffExp/output/CvPs_upSP_kallisto.rds")
CvPs_downSP <- readRDS("diffExp/output/CvPs_downSP_kallisto.rds")

library(ImpulseDE2)
#Obtain significant genes of interest
impulse_resultsFull <- readRDS("diffExp/output/impulseDiffGenesFullvsMono_kallisto.rds")
impulse_resultsBt <- readRDS("diffExp/output/impulseDiffGenesCvBtvsMono_kallisto.rds")
impulse_resultsPs <- readRDS("diffExp/output/impulseDiffGenesCvPsvsMono_kallisto.rds")

Fullgenes <- readRDS("diffExp/output/impulseResultsFull_FDR-LFC_kallisto.rds")
CvBtgenes <- readRDS("diffExp/output/impulseResultsCvBt_FDR-LFC_kallisto.rds")
CvPsgenes <- readRDS("diffExp/output/impulseResultsCvPs_FDR-LFC_kallisto.rds")

#Confirm selected genes are differentially expressed.
table(is.na(match(rownames(Full_upSP),Fullgenes)))
table(is.na(match(rownames(Full_downSP),Fullgenes)))
table(is.na(match(rownames(CvBt_upSP),CvBtgenes)))
table(is.na(match(rownames(CvBt_downSP),CvBtgenes)))
table(is.na(match(rownames(CvPs_upSP),CvPsgenes)))
table(is.na(match(rownames(CvPs_downSP),CvPsgenes)))

#All match differentially expressed genes

#Plot COG categories by up and down regulated genes compared to monoculture

Cv_FullU <- Cv_Cog[which(Cv_Cog$Locus %in% row.names(Full_upSP)),]
Cv_FullD <- Cv_Cog[which(Cv_Cog$Locus %in% row.names(Full_downSP)),]
Cv_BtU <- Cv_Cog[which(Cv_Cog$Locus %in% row.names(CvBt_upSP)),]
Cv_BtD <- Cv_Cog[which(Cv_Cog$Locus %in% row.names(CvBt_downSP)),]
Cv_PsU <- Cv_Cog[which(Cv_Cog$Locus %in% row.names(CvPs_upSP)),]
Cv_PsD <- Cv_Cog[which(Cv_Cog$Locus %in% row.names(CvPs_downSP)),]

Cv_FullU_sum <- as.data.frame(table(Cv_FullU$COG)/table(Cv_Cog$COG)*100)
Cv_FullD_sum <- as.data.frame(table(Cv_FullD$COG)/table(Cv_Cog$COG)*100)
Cv_BtU_sum <- as.data.frame(table(Cv_BtU$COG)/table(Cv_Cog$COG)*100)
Cv_BtD_sum <- as.data.frame(table(Cv_BtD$COG)/table(Cv_Cog$COG)*100)
Cv_PsU_sum <- as.data.frame(table(Cv_PsU$COG)/table(Cv_Cog$COG)*100)
Cv_PsD_sum <- as.data.frame(table(Cv_PsD$COG)/table(Cv_Cog$COG)*100)

#Merge dataframes
library(tidyverse)
Cv_Results <- list(Cv_FullU_sum, Cv_FullD_sum, Cv_BtU_sum,Cv_BtD_sum,Cv_PsU_sum,Cv_PsD_sum) %>% reduce(left_join, by = "Var1")

#Change names
names(Cv_Results) <- c("Category","Full_Upregulation","Full_Downregulation","Bt_Upregulation","Bt_Downregulation","Ps_Upregulation","Ps_Downregulation")

#Reshape
library(reshape2)
library(dplyr)
Cv_ResultsF <- melt(Cv_Results)

#Add condition for facet
Cv_ResultsF$Condition <- c(rep("CvPsBt",nrow(Cv_Results)*2),rep("CvBt",nrow(Cv_Results)*2),rep("CvPs",nrow(Cv_Results)*2))
#Add regulation for facet
Cv_ResultsF$Regulation <- rep(c(rep("Upregulation",nrow(Cv_Results)),rep("Downregulation",nrow(Cv_Results))),3)

#Remove NaN rows
Cv_ResultsF <- Cv_ResultsF[complete.cases(Cv_ResultsF), ]

#Drop levels with no data
Cv_ResultsF$Category <- factor(Cv_ResultsF$Category)

#Let's plot stacked bar
library(ggplot2)

## set the order categories
temp_df_1 <-
  Cv_ResultsF %>%
  arrange(Category)
the_order <- temp_df_1$Category

#Plot
gg <- Cv_ResultsF %>% ggplot(aes(x=Category,weight=value,fill=Regulation)) + geom_bar(position="dodge") +  scale_x_discrete(limits = levels(the_order)) +
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
Cv_ResultsF_new <- Cv_ResultsF

Cv_ResultsF_new$Regulation <- factor(Cv_ResultsF_new$Regulation, levels = c("Upregulation","Downregulation"))
Cv_ResultsF_new$Condition = factor(Cv_ResultsF$Condition, levels=c("CvPs","CvBt","CvPsBt"))

#Re-plot
gg <- Cv_ResultsF_new %>% ggplot(aes(x=Category,weight=value,fill=Regulation)) + geom_bar(position="dodge") +  scale_x_discrete(limits = levels(the_order)) +
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

ggsave("Cv_COG_categories.eps",plot=gg_final,device="eps",width=30, units="cm",dpi=300)

######### Coculture comparison

library(DESeq2)
library(ImpulseDE2)

#Now we can start Impulse analysis

#Load objects containing size factor and dispersion estimates

dispersions_deseq <- readRDS("diffExp/initial_files/dispersions_deseq_kallistoNCBI.rds")
sizeEst_deseq <- readRDS("diffExp/initial_files/sizeEst_deseq_kallistoNCBI.rds")
diffGenesF <- readRDS("diffExp/initial_files/diffGenesF_kallistoNCBI.rds")

###Now we can prep for Impulse
library(ImpulseDE2)
lib <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/initial_files/Cviolaceum_LIBRARIES.txt", sep="\t",header=TRUE)
counts <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/initial_files/Cvraw_counts_kallisoNCBI.txt",sep="\t",header=TRUE)

#Remove rRNA, tRNA, and miscRNA
library(rtracklayer)
gff <- readGFF("initial_files/Cviolaceum_genomic.gff") #for NCBI
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
Extract <- c("3mem","Bt","Ps","mono")
keywords <- str_extract(lib$sampleName, paste(Extract,collapse="|"))

#######CvBt vs CvPs "control"####################
mapping <- read.csv(file="/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Chromobacterium/final_analyses_Rev/initial_files/Cv_diffExp_CvBt_v_CvPs.csv",header=TRUE,sep=",")

keys <- which(keywords %in% c("Bt","Ps"))
samps <- as.vector(lib$sampleName[keys])
samples <- match(samps,lib$sampleName)
geneMat <- genecounts.sort[,samples]
genemat <- as.data.frame(geneMat)
align <- as.vector(mapping[,1])

diffGenesCo=diffGenesF[,match(align, colnames(diffGenesF))]
colnames(diffGenesCo) <- mapping$Sample
#rownames(diffGenes) <- rownames(genecounts.sort)
mapping$time <- as.factor(mapping$time)

design <- data.frame("Sample"=mapping$Sample,"Condition"=mapping$class,
"Time"=mapping$time, "Batch"=rep("B_NULL",nrow(mapping),row.names=mapping$Sample))
design$Time <- as.numeric(as.character(design$Time))*60
design$Condition <- sub(pattern="con",replacement="control", design$Condition)
design$Condition <- sub(pattern="exp",replacement="case", design$Condition)
diffGenesCo <- as.matrix(diffGenesCo)

names(dispersions_deseq) <- row.names(diffGenesCo)
#For size estimates, we first need to parse out our samples of interest
sizeFactors <- sizeEst_deseq[which(names(sizeEst_deseq) %in% align)]
names(sizeFactors) <- mapping$Sample

#Identify DEGs
impulse_resultsCo_trends <- runImpulseDE2(matCountData = diffGenesCo,dfAnnotation =design,
boolCaseCtrl = TRUE,vecDispersionsExternal=dispersions_deseq,
vecSizeFactorsExternal=sizeFactors,scaQThres=0.01,boolIdentifyTransients=TRUE)

saveRDS(impulse_resultsCo_trends,"diffExp/output/impulseDiffGenesCvBtvsCvPs_kallisto.rds")

###########Filter DEGs by LFC threshold

#Now that we have significant genes determined by ImpulseDE2 model, let's filter DEGs by DEGs that have at least one TP with LFC >= 1 or <= -1
library(ImpulseDE2)

impulse_resultsCo <- readRDS("diffExp/output/impulseDiffGenesCvBtvsCvPs_kallisto.rds")

load("/mnt/research/ShadeLab/Chodkowski/JGI_SynCom/RNAseq/Burkholderia/final_analyses_Rev/lfc/output/BthailandensisLFC_kallisto.RData")

lfc_analysisCo <- cbind(lfc_CvBtInt_comparedtoCvPsInt.fc,lfc_CvBt25_comparedtoCvPs25.fc,lfc_CvBt30_comparedtoCvPs30.fc,
                          lfc_CvBt35_comparedtoCvPs35.fc,lfc_CvBt40_comparedtoCvPs40.fc,lfc_CvBt45_comparedtoCvPs45.fc)
#Just pick one sample to place the row names back in.
row.names(lfc_analysisCo) <- row.names(lfc_CvBtInt_comparedtoCvPsInt)

CoSG <- lfc_analysisCo[match(impulse_resultsCo$vecDEGenes,row.names(lfc_analysisCo)),]

CoSG_lfc_analysis <- as.data.frame(CoSG)

library(dplyr)
library(tibble)

#####Coculture comparisons LFC thresholds
#up-regulated
Co_up <- CoSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,any_vars(. >= 1 )) %>%
    column_to_rownames('gene')

nrow(Co_up)

Co_down <- CoSG_lfc_analysis %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,any_vars(. <= -1 )) %>%
    column_to_rownames('gene')

nrow(Co_down)

#Find duplicate genes from both up/downregulated
duprows <- rownames(Co_down) %in% rownames(Co_up)

#Combine lists
Co_LFC_thresh <- rbind(Co_up,Co_down[!duprows,]) #1762 of 2915 DEGs

#1762 of 2915 DEGs passed LFC threshold

#Save these as that final DEGs with a LFC threshold
saveRDS(rownames(Co_LFC_thresh),"diffExp/output/impulseResultsCo_FDR-LFC_kallisto.rds")

#############Obtain gene trends in SP

#Exponential Phase TP removal. It happens to be the first column in FullSG_lfc_analysis
CoSG_lfc_analysisSP <- CoSG_lfc_analysis[,-1]

#up-regulated in stationary phase
Co_upSP <- CoSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. > 0 )) %>%
    column_to_rownames('gene')

nrow(Co_upSP)

#downregulated in stationary phase
Co_downSP <-  CoSG_lfc_analysisSP %>%
    rownames_to_column(as.character('gene')) %>%
    filter_if(is.numeric,all_vars(. < 0 )) %>%
    column_to_rownames('gene')

nrow(Co_downSP)

#Filter gene trends that passed LFC threshold
Co_upSP_F <- Co_LFC_thresh[which(rownames(Co_LFC_thresh) %in% rownames(Co_upSP)),]
Co_downSP_F <- Co_LFC_thresh[which(rownames(Co_LFC_thresh) %in% rownames(Co_downSP)),]

saveRDS(Co_upSP_F,"diffExp/output/Co_upSP_kallisto.rds")
saveRDS(Co_downSP_F,"diffExp/output/Co_downSP_kallisto.rds")

###############COG analysis##########################

Cv_Cog_pre <- read.csv("diffExp/initial_files/Cv_COG_EggNogMapper.csv",header=TRUE,sep=",")

#Remove Categories with less that 10 genes
Cv_Cog <- Cv_Cog_pre[!Cv_Cog_pre$COG %in% names(which(table(Cv_Cog_pre$COG) < 10)), ]

#Analyzing coculture conditions COG categories to monoculture- edit to only include signicant genes
Co_upSP <- readRDS("diffExp/output/Co_upSP_kallisto.rds")
Co_downSP <- readRDS("diffExp/output/Co_downSP_kallisto.rds")

library(ImpulseDE2)
#OCvain significant genes of interest
Cogenes <- readRDS("diffExp/output/impulseResultsCo_FDR-LFC_kallisto.rds")

#Confirm selected genes are differentially expressed.
table(is.na(match(rownames(Co_upSP),Cogenes)))
table(is.na(match(rownames(Co_downSP),Cogenes)))

#All match differentially expressed genes

#Plot COG categories by up and down regulated genes compared to monoculture

Cv_FullU <- Cv_Cog[which(Cv_Cog$Locus %in% row.names(Co_upSP)),]
Cv_FullD <- Cv_Cog[which(Cv_Cog$Locus %in% row.names(Co_downSP)),]

Cv_FullU_sum <- as.data.frame(table(Cv_FullU$COG)/table(Cv_Cog$COG)*100)
Cv_FullD_sum <- as.data.frame(table(Cv_FullD$COG)/table(Cv_Cog$COG)*100)

#Merge dataframes
library(tidyverse)
Cv_Results <- list(Cv_FullU_sum, Cv_FullD_sum) %>% reduce(left_join, by = "Var1")

#Change names
names(Cv_Results) <- c("Category","Coculture_Upregulation","Coculture_Downregulation")

#Reshape
library(reshape2)
library(dplyr)
Cv_ResultsF <- melt(Cv_Results)

#Add condition for facet
Cv_ResultsF$Condition <- c(rep("CvBt-CvPs",nrow(Cv_Results)*2))
#Add regulation for facet
Cv_ResultsF$Regulation <- rep(c(rep("Upregulation",nrow(Cv_Results)),rep("Downregulation",nrow(Cv_Results))),1)

#Remove NaN rows
Cv_ResultsF <- Cv_ResultsF[complete.cases(Cv_ResultsF), ]

#Drop levels with no data
Cv_ResultsF$Category <- factor(Cv_ResultsF$Category)

#Let's plot stacked bar
library(ggplot2)

## set the order categories
temp_df_1 <-
  Cv_ResultsF %>%
  arrange(Category)
the_order <- temp_df_1$Category

temp_df_2 <-
  Cv_ResultsF %>%
  arrange(desc(Regulation))
the_order2 <- temp_df_2$Regulation

#Set order of facets
#temp$size_f = factor(temp$size, levels=c('50%','100%','150%','200%'))

#Set order of regulation
#Cv_ResultsF$new = factor(Cv_ResultsF$Regulation, levels=c("Upregulation","Downregulation"), labels=c("Upregulation","Downregulation"))

#Plot
gg <- Cv_ResultsF %>% ggplot(aes(x=Category,weight=value,fill=Regulation)) + geom_bar(position="dodge") +  scale_x_discrete(limits = levels(the_order)) +
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
Cv_ResultsF_new <- Cv_ResultsF

Cv_ResultsF_new$Regulation <- factor(Cv_ResultsF_new$Regulation, levels = c("Upregulation","Downregulation"))
Cv_ResultsF_new$Condition = factor(Cv_ResultsF$Condition, levels=c("CvBt-CvPs"))

#Re-plot
gg <- Cv_ResultsF_new %>% ggplot(aes(x=Category,weight=value,fill=Regulation)) + geom_bar(position="dodge") +  scale_x_discrete(limits = levels(the_order)) +
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


gg <- Cv_ResultsF_new %>% ggplot(aes(x=Category,weight=value,fill=Regulation)) + geom_bar(position="dodge") +  scale_x_discrete(limits = levels(the_order)) +
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

ggsave("CvBt-CvPs_COG_categories.eps",plot=gg_final,device="eps",width=30, units="cm",dpi=300)


```
