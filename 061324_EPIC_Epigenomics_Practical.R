#Student: Pierre Wensel
#UVic-UCC Course: Epigenomcs
#Date: May 31, 2024
#Task: DNA Methylation Epics Array Practical Exercise

#rm(list=ls()) 
#Introduction
#In medicine, hepatopulmonary syndrome is a syndrome of shortness of 
#breath and hypoxemia (low oxygen levels in the blood of the arteries) caused 
#by vasodilation (broadening of the blood vessels) in the lungs of patients with liver 
#disease. Dyspnea and hypoxemia are worse in the upright position (which is 
#called platypnea and orthodeoxia, respectively). 

#Hepatopulmonary syndrome (HPS) is the most common cause
#of respiratory insufficiency in patients with chronic liver disease. 
#It is characterised by a gas exchange abnormality caused
#by intrapulmonary vascular dilatations (IPVD) in patients with
#liver diseases. It occurs in 5-32% of liver transplant candidates.1
#Abnormal oxygenation is defined by elevated alveolar-arterial
#oxygen gradient (>15 mmHg or >20 mmHg in patients >64 years, respectively) 
#while breathing room air in the sitting position at rest. IPVD is usually 
#diagnosed by contrast-enhanced transthoracic echocardiography. 
#The severity of HPS is classified according to the degree of hypoxemia.
#Every patient with advanced liver disease that is undergoing
#evaluation for liver transplantation or suffers from dyspnoea
#should be screened for the presence of HPS. An established
#screening tool is pulse oximetry.In this data set we have 7 patients 
#and 5 controls with an EPIC array for which analysis of differential DNA mythylation  


#The following publication was cited and cross-referenceed to identify potential biomarkers:
#Mendoza N, Rivas E, Rodriguez-Roisin R, Garcia T, Bruguera M, Agusti A, et al. (2021) Liver
#epigenome changes in patients with hepatopulmonary syndrome: A pilot study. PLoS
#ONE 16(2): e0245046. https://doi.org/10.1371/journal.pone.0245046

#None of the 802,688 DNA probes analyzed in the case control comparison achieved a significant False Discovery Rate
#(FDR). WGCNA identified 5 co-methylated gene-modules associated to HPS markers, mainly related to nervous and 
#neuroendocrine system, apoptotic processes, gut bacterial translocation, angiogenesis and 
#vascular remodeling ontologies. To conclude, HPS is associated with nervous/neuroendocrine 
#system and vascular remodeling related liver epigenetic changes

#Workflow


#####INSTALL AND/OR LOAD NECESSARY PACKAGES

#devtools::install_version("dbplyr", version = "2.3.4")
#devtools::install_github("TomKellyGenetics/heatmap.2x", ref="master")
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
#BiocManager::install("minfi")

#BiocManager::install("clusterProfiler")
#BiocManager::install("ExperimentHub")
#if (!require("devtools")) install.packages("devtools")
#library(devtools)
#install_github("husson/FactoMineR")
#install.packages("remotes")
#remotes::install_github("metamaden/recountmethylation")
#BiocManager::install("topGO")
#install.packages("remotes")
#remotes::install_github("perishky/dmrff")
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#library(doParallel)
#registerDoParallel(cores = 3)
library(heatmap.2x) 
library(gplots)
library(DMRcate)
library(dbplyr)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOstats)
library(GO.db)
library(DOSE)
library(annotate)
library(missMethyl)
library(dplyr)
library(ExperimentHub)
library(ggplot2)
library(limma)
library(recountmethylation)
#library(httr)
library(topGO)
library(gplots)
library(DOSE)
library(missMethyl)
library(DMRcate)
library(dmrff)
library(GenomicRanges) 
##############SET DIRECTORIES, RENAME AND READ IN TARGETS FILE AND IDAT DATA FILES:


setwd("C:/Users/User/Desktop/R/Metabolomics/")
getwd()
idat.folder <- "Epigenomics_EPIC/idats"
#change for your path directory
#Note: I renamed the course-provided samples.csv file to "targets.csv
targets <- read.metharray.sheet(base=idat.folder)
class(targets)
targets
#Note: Unlike previous examples involving Illumina440K, there is no "Name" column in "targets.csv"
#Adding a name column with identifiers H and C

#data.met = readRDS(file = "data.met.RDS")
# met matrix
#met <- as.data.frame(SummarizedExperiment::assay(data.met))
## remove probes with NA
#probe.na <- rowSums(is.na(met))
#table(probe.na == 0)
#FALSE   TRUE 
#103553 382024 
# chose those has no NA values in rows
#probe <- probe.na[probe.na == 0]
#met <- met[row.names(met) %in% names(probe), ]

targets$name<-targets$Basename 
targets$name<-paste0(targets$name,c("H","C","C","C","C","H","H","C","H","H","H","H"))

targets
#test1<-test %>% mutate(across(c('Basename'), substr, 24, nchar(Basename)))

targets$Basename <- paste0(c("Class/"),targets$Basename)
#cat("Class/", targets$Basename, sep="")
#targets[] <- lapply(targets, function(x) paste('Class/', x))
#targets[] <- Map(paste, 'Class/', targets)
targets
targets$Basename<-gsub("Class/","Epigenomics_EPIC/idats/",targets$Basename)
targets
#targets$Basename<-gsub("Class/","/Epigenomics_EPIC/idats/",targets$Basename)
#targets$Basename<-gsub("Class/","C:/Users/click/OneDrive/Desktop/R_TRANSCIPTOMICS/Epigenomics_EPIC/idats/",targets$Basename)

###loading data

rgset <- read.metharray.exp(targets = targets,verbose = T)
rgset

####quality control first day

#The class of RGSet is a RGChannelSet object. 
#This is the initial object of a minfi analysis that contains the raw intensities in the green and red channels. 
#Note that this object contains the intensities of the internal control probes as well. 

##################DATA EXPLORATION AND CONVERSION OF OBJECTS ULTIAMTELY FROM rgset to GEnomicRAtioSet:
phenoData <- pData(rgset)

head(phenoData)

#The RGChannelSet stores also a manifest object that contains the probe design information of the array:

manifest <- getManifest(rgset)
manifest

####Probes information

getProbeInfo(manifest, type = "I")

getProbeInfo(manifest, type = "II")


##A MethylSet objects contains only the methylated and unmethylated signals

MSet <- preprocessRaw(rgset)
MSet

#A RatioSet object is a class designed to store Beta values and/or M values instead of the methylated and unmethylated signals. 
#An optional copy number matrix, CN, the sum of the methylated and unmethylated signals, can be also stored. 
#Mapping a MethylSet to a RatioSet may be irreversible, i.e. one cannot be guranteed to retrieve the methylated and unmethylated signals from a RatioSet.
#A RatioSet can be created with the function ratioConvert:

RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
RSet

#The functions getBeta, getM and getCN return respectively the Beta value matrix, M value matrix and the Copy Number matrix.
beta <- getBeta(RSet)
beta
head(beta)

#GenomicRatioSet

#The function mapToGenome applied to a RatioSet object will add genomic coordinates to each probe together with some additional annotation information.
#The output object is a GenomicRatioSet (class holding M or/and Beta values together with associated genomic coordinates). 
#It is possible to merge the manifest object with the genomic locations by setting the option mergeManifest to TRUE.

GRset <- mapToGenome(RSet)
GRset

beta <- getBeta(GRset)
head(beta)

M <- getM(GRset)
class(M)

head(M)

CN <- getCN(GRset)
class(CN)

head(CN)

sampleNames <- sampleNames(GRset)
class(sampleNames)

sampleNames


probeNames <- featureNames(GRset)
class(probeNames)

head(probeNames)


pheno <- pData(GRset)
pheno

gr <- granges(GRset)
head(gr, n= 3)

##Full annotation
annotation <- getAnnotation(GRset)

#head(annotation)

names(annotation)

#NOT DEMONSTRATED IN LECTURE VIDEO BUT INCLUDED IN CODE?
levels(as.factor(unlist(strsplit(annotation$UCSC_RefGene_Group,";"))))
#[1] "1stExon" "3'UTR"   "5'UTR"   "Body"    "ExonBnd" "TSS1500" "TSS200"

################QUALITY CONTROL

#minfi provides a simple quality control plot that uses the log median intensity in both the methylated (M) and unmethylated (U) channels.
#When plotting these two medians against each other, it has been observed that good samples cluster together, 
#while failed samples tend to separate and have lower median intensities. 
#In order to obtain the methylated and unmethylated signals, we need to convert the RGChannelSet to an object containing the methylated and unmethylated signals using the function preprocessRaw. It takes as input a RGChannelSet and converts the red and green intensities to methylated and unmethylated signals according to the special 450K probe design, and returns the converted signals in a new object of class MethylSet. It does not perform any normalization.

#The accessors getMeth and getUnmeth can be used to get the methylated and unmethylated intensities matrices:

head(getMeth(MSet))
head(getUnmeth(MSet))

#The functions getQC and plotQC are designed to extract and plot the quality control information from the MethylSet:
qc <- getQC(MSet)
qc

#png("QC.plot2.png")
plotQC(qc)
#dev.off()

#Above a certain dashed-line threshold, the probes are hybridizing well
#To further explore the quality of the samples, it is useful to look at the Beta value densities of the samples, with the option to color the densities by group:
densityPlot(MSet, sampGroups = targets$Disease_state)
#binomial distribution suggests good quality of hybridization

phenoData2 <- pData(MSet)
densityPlot(MSet, sampGroups = phenoData2$Disease_state)

#The EPIC array contains several internal control probes that can be used to assess the quality control of different sample 
#preparation steps (bisulfite conversion, hybridization, etc.). 
#The values of these control probes are stored in the initial RGChannelSet 
#and can be plotted by using the function controlStripPlot and by specifying the control probe type:
#png("Bis_con.png")
controlStripPlot(rgset, controls="BISULFITE CONVERSION II")

#Sex prediction

#By looking at the median total intensity of the X chromosome-mapped probes, denoted med(X)med(X), 
#and the median total intensity of the Y-chromosome-mapped probes, denoted med(Y)med(Y),
#one can observe two different clusters of points corresponding to which gender the samples belong to.
#To predict the gender, minfi separates the points by using a cutoff on log2med(Y)

#The default cutoff is ???2???2.plotSex expects as input a [Genomic]MethylSet, yet getSex returns a DataFrame. 
#Since the algorithm needs to map probes to the X-chr and to the Y-chr, the input of the function 
#getSex needs to be a GenomicMethylSet or a GenomicRatioSet.

predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
#predictedSex <- getSex(GRset, cutoff = -2)??

predictedSex

phenoData$gender

#Genders match!!!!!!!!!!!

#To choose the cutoff to separate the two gender clusters, one can plot med(Y)med(Y) against med(Y)med(Y) with the function plotSex:
 
GRset<-addSex(GRset)####add sex to GRatio object

#pdf("sex.pdf")
plotSex(GRset)
#plotSex(getSex(GRset, cutoff = -2))
#dev.off()
#Evidently, sex is clustering
getSex
#function (object = NULL, cutoff = -2) 
#{
#  .isGenomicOrStop(object)
#  if (is(object, "GenomicMethylSet")) 
#    CN <- getCN(object)
#  if (is(object, "GenomicRatioSet")) 
#    CN <- getCN(object)
#  xIndex <- which(seqnames(object) == "chrX")
#  yIndex <- which(seqnames(object) == "chrY")
#  .getSex(CN = CN, xIndex = xIndex, yIndex = yIndex, cutoff = cutoff)
#}
#<bytecode: 0x000002792cf22378>
#  <environment: namespace:minfi>

####detection pvals

detP<-detectionP(rgset)
#Plotting the mean detection p-value for each sample will allow us to gauge whether any samples have many failed probes 
#- this will be indicated by a large mean detection p-value. Samples with mean detection p-values exceeding a cutoff 
#such as 0.05 can be excluded from further analysis.

dim(detP)
#[1] 866836     12

head(detP)

# examine mean detection p-values across all samples to identify any failed samples
#pdf("pvals.detect.pdf")
barplot(colMeans(detP),col=factor(targets$Disease_state),las=2,cex.names=0.8,main="Mean detection p-values")
abline(h=0.05,col="red")
#Above 2 lines worked, skipped pdf line and devoff
#barplot(colMeans(detP), las=2, cex.names=0.8, ylab="Mean detection p-values")
#dev.off()
 
###In case we want to remove samples
keep<-colMeans(detP) < 0.01
rgset<-rgset[,keep]
targets = targets[keep,]
detP = detP[,keep]
dim(detP)

###probe
keep <- rowSums(detP < 0.01) == ncol(rgset)
#keep<-rowMeans(detP) < 0.01 #FROM 450K script

table(keep)

not.keep<-names(keep[keep=="FALSE"])

rgset <- rgset[keep,]

dim(rgset)

#################NORMALIZATION

#BLOCK2: Analysis (in order to standardize, normalize by Quantiles). Do not expect 
#huge differences in this case, hence use nominal pvalue to establish differential 
#analysis control instead of FDR.
#1. How many CpgS are differentially methylated without considering covariates 
#(age, sex…)? That means all, promoter, and promoter plus island.

#If there exist global biological methylation differences between your samples, as for instance a dataset 
#with cancer and normal samples, or a dataset with different tissues/cell types, use the preprocessFunnorm function 
#as it is aimed for such datasets. On the other hand, if you do not expect global differences between your samples, 
#for instance a blood dataset, or one-tissue dataset, use the preprocessQuantile function. 
#In our experience, these two normalization procedures perform always better than the functions preprocessRaw, 
#preprocessIllumina and preprocessSWAN discussed below. For convenience, these functions are still implemented 
#in the minfi package.”

#gRatioSet.Illumina <- preprocessIllumina(rgset)##Illumina
#gRatioSet.SWAN <- preprocessSWAN(rgset,verbose = T)##SWAN
#gRatio.Fun<-preprocessFunnorm(rgset,verbose = T)

#The preprocessQuantile function (Input: RGChannelSet,Output: GenomicRatioSet)implements stratified quantile normalization preprocessing. 
#The normalization procedure is applied to the Meth and Unmeth intensities separately. The distribution of type I and type II signals is 
#forced to be the same by first quantile normalizing the type II probes across samples and then interpolating a reference distribution to
#which we normalize the type I probes. Since probe types and probe regions are confounded and we know that DNA methylation varies across 
#regions we stratify the probes by region before applying this interpolation. Note that this algorithm relies on the assumptions necessary 
#for quantile normalization to be applicable and thus is not recommended for cases where global changes are expected such as in cancer-normal 
#comparisons as these would be removed by the normalization.As we are comparing different x types, which are globally relatively similar, 
#we will apply the preprocessQuantile method to our data.

gRatioSet.quantile <- preprocessQuantile(rgset,verbose = T)##SQN
class(gRatioSet.quantile) 
gRatioSet.quantile

#Compare with the unnormalized data to visualize the effect of the normalization. 
#First a comparison of the Beta distributions for the different probe designs. 
#This will give an indication of the effectiveness of the within-array normalization.

par(mfrow=c(1,2))
# Plot distributions prior to normalization for sample 1
plotBetasByType(MSet[,1],main="Raw")
# The normalized object is a GenomicRatioSet which does not contain
# the necessary probe info, we need to extract this from the MethylSet first.
typeI <- getProbeInfo(MSet, type = "I")[, c("Name","nCpG")]
typeII <- getProbeInfo(MSet, type = "II")[, c("Name","nCpG")]
probeTypes <- rbind(typeI, typeII)
probeTypes$Type <- rep(x = c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
# Now plot the distributions of the normalized data for sample 1
plotBetasByType(getBeta(gRatioSet.quantile)[,1], probeTypes = probeTypes, main="Normalized",)
#Evidently, the normalization brought the distributions closer to each other 
#Now let’s see how the between-array normalization worked…

# visualise what the data looks like before and after normalization
par(mfrow=c(1,2))
densityPlot(rgset, sampGroups=targets$Disease_state,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Disease_state)),text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(gRatioSet.quantile), sampGroups=targets$Sample_Group,main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Disease_state)),text.col=brewer.pal(8,"Dark2"))

#After normalization of  data is a good time to look at the similarities and differences between the various samples. 
#One way to do this is by creating a MDS or Multi-Dimensional Scaling plot. This is a method to graphically represent relationships
#between objects (here the different samples) in multidimensional space onto 2 or 3 dimensional space. 
# MDS plots to look at largest sources of variation
# Create color panel

pal <- brewer.pal(8,"Dark2")
# Plot figures
par(mfrow=c(1,2))
plotMDS(getM(gRatioSet.quantile), top=1000, gene.selection="common",
        col=pal[factor(targets$Disease_state)], dim=c(1,2))
legend("top", legend=levels(factor(targets$name)), text.col=pal,
       bg="white", cex=0.7)

#Evidently, 44% of total variance is explained by first 2 Principal Components (PC1=27%, PC2=17%) based on eigenvalues og eigenvectors


############################################################################
densityPlot(getBeta(gRatioSet.quantile), sampGroups = targets$Disease_state)

#UNUSED TEST CODE 
####Correlations
sample.quantile<-colMeans(getBeta(gRatioSet.quantile),na.rm = T)
#sample.illumina<-colMeans(getBeta(gRatioSet.Illumina),na.rm = T)
#sample.swan<-colMeans(getBeta(gRatioSet.SWAN),na.rm = T)


#mod1 <-lm(as.numeric(sample.quantile)~as.numeric(sample.illumina))## 
#modsum<-summary(mod1)

#r2 <- cor.test(as.numeric(sample.quantile),as.numeric(sample.illumina),method="spearman")$estimate
#my.p<-cor.test(as.numeric(sample.quantile),as.numeric(sample.illumina),method="spearman")$p.value
#my.p<-signif(my.p, digits=3)

#png("Quantiles_Illumina.png")
#mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
#plot(as.numeric(sample.quantile)~as.numeric(sample.illumina), main=paste("Correlation", "p-value", my.p, sep=" "), pch=20,col="grey40")

#    ,xlim=c(0,10000),ylim=c(0,10000))
#  plot(as.numeric(cnts[,1])~as.numeric(cnts[,2]), main=paste("Correlation","Pulse Replicates Global", "p-value", my.p, sep=" "), xlab="Ctrl1", ylab="Ctrl2", pch=20)
#abline(mod1, col="red")
# abline(h=0,col="blue1")
#abline(v=0,col="blue1")
#legend('topleft', legend = mylabel, bty = 'n')

#dev.off()

########################################PREPARE GENOME ANNOTATION 

ann850k<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
class(ann850k)

ann850k<-as.data.frame(ann850k)
#head(ann850k)

colnames(ann850k)

############################FILTERING

#Poor performing probes can obscure the biological signals in the data and are generally filtered out prior to 
#differential methylation analysis. As the signal from these probes is unreliable, by removing them we perform fewer
#statistical tests and thus lower the multiple testing penalty. We filter out probes that have failed in one or more 
#samples based on detection p-value.

#UNUSED CODE
#???????????????????START
# ensure probes are in the same order in the mSetSq and detP objects
#detP <- detectionP(rgSet)
#detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

# remove any probes that have failed in one or more samples; this next line
# checks for each row of detP whether the number of values < 0.01 is equal
# to the number of samples (TRUE) or not (FALSE)
#keep <- rowSums(detP < 0.01) == ncol(mSetSq)
#table(keep)
# Subset the GenomicRatioSet
#mSetSqFlt <- mSetSq[keep,]
#mSetSqFlt
#?????????????????????STOP

#Not interested in imprinting or sexual methylation studies. We are instead interested in disease state.etc.
keep <- !(featureNames(gRatioSet.quantile) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])

table(keep)

gRatioSet.quantile <- gRatioSet.quantile[keep,]
gRatioSet.quantile

##Remove probes with SNPs at CpG or SBE site (These are variable among individuals, so we remove the corresponding probes for those Cpg islands)
#####get snp info

snp.info<-getSnpInfo(gRatioSet.quantile)
snp.info

####mapping

gRatioSet.quantile<- dropLociWithSnps(gRatioSet.quantile)
gRatioSet.quantile
densityPlot(getBeta(gRatioSet.quantile), sampGroups = targets$Disease_state)

#UNUSED CODE
#?????????START
#phenoData <- pData(MSet)
#densityPlot(MSet, sampGroups = phenoData$Sample_Group)
#??????????STOP

####get values to proceed with analysis
betas<-getBeta(gRatioSet.quantile)
head(betas)

dim(betas)

betas <- rmSNPandCH(betas, dist=2, mafcut=0.05,rmcrosshyb = T,rmXY = T) ### get out SNPs of my samples
#Probe IDs from EPICv1 or earlier detected. Proceeding...
#Cannot connect to ExperimentHub server, using 'localHub=TRUE' instead
#Using 'localHub=TRUE'
#If offline, please also see BiocManager vignette section on offline use
#snapshotDate(): 2024-05-28
#see ?DMRcatedata and browseVignettes('DMRcatedata') for documentation
#loading from cache
#see ?DMRcatedata and browseVignettes('DMRcatedata') for documentation
#loading from cache
#see ?DMRcatedata and browseVignettes('DMRcatedata') for documentation
#loading from cache

head(betas)

dim(betas)

#Calculate M-values for statistical analysis: as previously mentioned, M-values have nicer statistical properties 
#and are thus better for use in statistical analysis of methylation data
M<-getM(gRatioSet.quantile)
head(M)

dim(M)

#Beta and M dimensions are NO LONGER identical. Is this a problem? There is no
#mention or example of similarly adjusting M in notes or videos. Therefore, I will
#leave this alone

#####
#ann450k<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#M <- getM(gRatioSet.quantile)

#Confiming that beta names and target names are in SAME ORDER:
colnames(betas)
targets$name
#These are in the same order. Therefore, I can safely rename  betas columns
colnames(betas)<-targets$name

#UNUSED CODE
#???????????????START
#ALTERNATIVELY:
###put same order
#betas<-betas[,match(rownames(pData(gRatioSet.quantile)),colnames(betas))]
#colnames(betas)<-pData(gRatioSet.quantile)$name
#betas<-DMRcate::rmSNPandCH(betas, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY=TRUE)
#???????????????STOP

#Remove probes with cross-reactivity

#probes2remove<-read.csv("48639-non-specific-probes-Illumina450k.txt")
#betas<-betas[!(rownames(betas) %in% probes2remove$TargetID), ]
#ann450k<-as.data.frame(ann450k)

#Access lists of cross-reactive probes with the recountmethylation function using get_crossreactive_cpgs() (see ?get_crossreactive_cpgs for details).
#Pidsley et al 2016 cross-reactive EPIC probes. Databases accessed with recountmethylation contain data from GEO (ncbi.nlm.nih.gov/geo/), a live public database where alterations to 
#online records can cause discrepancies with stored data over time. We cannot guarantee the accuracy of stored data, and advise users 
#cross-check their findings with latest available records

crossreactive_cpgs<-get_crossreactive_cpgs("epic")
class(crossreactive_cpgs)
dim(crossreactive_cpgs)

head(crossreactive_cpgs)
length(get_crossreactive_cpgs("epic"))


dim(betas)
#betas<-betas[!(rownames(betas) %in% probes2remove$TargetID), ]
#betas_test<- betas[ -which(row.names(betas) %in% crossreactive_cpgs), ]
betas_test<-betas[!(rownames(betas) %in% crossreactive_cpgs), ]
dim(betas_test)

#Dimensions did not change after removal of cross-reactive probes

#Remove probes with NA
probe.na <- rowSums(is.na(betas))
table(probe.na == 0)

#EVidently, there are no NA probes
#Choose probes that have no NA values in rows

#UNUSED CODE
#????????START
#probe <- probe.na[probe.na == 0]
#betas <- betas[row.names(betas) %in% names(probe), ]
#Dimensions did not change after removal of NA probes
#???????STOP

#Once the data has been filtered and normalized, it is often useful to re-examine the MDS plots to see if the relationship between the samples has changed. 
#From the new MDS plots it is apparent that much of the inter-individual variation has been removed as this is no longer the first principal component,
#likely due to the removal of the SNP-affected CpG probes. However, the samples do still cluster by individual in the second dimension and thus a factor 
#for individual should still be included in the model.
#mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
par(mfrow=c(1,2))
plotMDS(getM(gRatioSet.quantile), top=1000, gene.selection="common",
        col=pal[factor(targets$Disease_state)], cex=0.8)
legend("right", legend=levels(factor(targets$Disease_state)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(getM(gRatioSet.quantile), top=1000, gene.selection="common",
        col=pal[factor(targets$name)])
legend("right", legend=levels(factor(targets$name)), text.col=pal,
       cex=0.7, bg="white")

#dev.off()

#Plotting PCA AGAIN
pca <- FactoMineR::PCA(t(betas), scale.unit = T, graph = F, ncp = 40)
factoextra::fviz_pca_ind(pca, axes = c(1,2), habillage=as.factor(targets$Disease_state), repel = T)

#########PROBE-WISE DIFFERENTIAL METHYLATION
#LIMMA

#NOTE: q-values are the name given to the adjusted p-values found using an optimised FDR approach.
# Since we do not expect huge differences in our analysis, we will use nominal pvalue to establish differential 
#analysis control instead of FDR:

#LIMMA MODEL1 APPROACH (treatment-contrasts parametrization)WITH ONE GROUP (Disease_state)

targets_copy<-targets

#Define the factor of interest
#targets_copy$Disease_state<-as.factor(targets_copy$Disease_state)
#targets_copy$Disease_state<-relevel(targets_copy$Disease_state,ref="Control")
#group<-as.factor(pData(gRatioSet.quantile)$Disease_state)
#group<-relevel(group,ref="Control")
targets_copy$Disease_state<- factor(targets_copy$Disease_state, levels=c("Control","HPS"))
class(targets_copy$Disease_state)
targets_copy$Disease_state

# Set up the design matrix for the Differential Methylation analysis
#design <- model.matrix(~0+cellType+individual, data=targets)
#colnames(design) <- c(levels(cellType),levels(individual)[-1])
#design1 <- model.matrix( ~ diseaseState, data = pData(dat) )
#design2 <- model.matrix( ~ 0 + dat$diseaseState )
design <- model.matrix(~Disease_state, data=targets_copy)
colnames(design) <- c("Control","HPSvsControl")
#colnames(design) <- levels(targets_copy$Disease_state)
#colnames(design)[2]<-"Comp"????
design

####generate adjusted beta matrix
#fit <- lmFit(M, design)
fit <- lmFit(betas, design)
head(fit$coefficients)
#Some second coefficients are negative and positive


#Summary of an lm result includes some readily-understood goodness of fit information that include
#Residual standard error, Multiple R-Squared, Adjusted R-Squared, F-statistic
#effect size (logFC) and 95% confidence intervals (CI.L, CI.R), from which we can back calculate the standard error.
#Standard deviations can be obtained as fit$sigma, or sqrt(fit$s2.post) for values after empirical Bayes shrinkage. 
head(fit$stdev.unscaled*fit$sigma)


#For mean expression within each group, standard errors of the mean are obtained as 
#head(fit$stdev.unscaled * sqrt(fit$s2.post))
#head(sqrt((fit$s2.post)*fit$stdev.unscaled))

# pooling of variance across like-genes
fit2 <- eBayes(fit)
#fit2$F FOR MULTIPLE LEVELS ANOVA not used here
#Extracting significantly methylated probes
#results <- topTable(fit1,coef="basal",adjust="BH",number=Inf,sort.by="P")
#results <- topTable(fit2,coef="Comp",num=dim(fit2)[1],sort.by="P",adjust.method = "BH")??????? 
results = topTable(fit2, coef=ncol(design), sort.by="P",number = nrow(betas), adjust.method = "BH")
dim(results)

head(results)

#Subsetting:

res<-subset(results,results$P.Value<.05)
class(res)
dim(res)
#[1] 40178     6
head(res)

#GARBAGE COLLECTION TO TROUBLESHOOT VOLCANO PLOT ERROR WITH DEPTH FUNCTION

gc()

# Visualization of top 10 most significantly differentially methylated CpGs 
par(mfrow=c(2,5))
sapply(rownames(res)[1:10], function(cpg){
  plotCpg(betas, cpg=cpg, pheno=targets_copy$Disease_state, type="categorical", ylab = "Beta values")
})

#Visualization of volcano plot
#making dataset
dat <- data.frame(foldchange = fit[["coefficients"]][,2], logPvalue =  -log10(fit2[["p.value"]][,2]))
dat$threshold <- as.factor(abs(dat$foldchange) < 0.4)

#Visualization
cols <- c("TRUE" = "grey", "FALSE" = "blue")
ggplot(data=dat, aes(x=foldchange, y = logPvalue, color=threshold)) +
  geom_point(alpha=.6, size=1.2) +
  scale_colour_manual(values = cols) +
  geom_vline(xintercept = 0.4, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = - 0.4, colour="#990000", linetype="dashed") +
  theme(legend.position="none") +
  xlab("Fold Change") +
  ylab("-log10 p value") +
  theme_bw() +
  theme(legend.position = "none")

results <- decideTests(fit2)
results


vennDiagram(results)

#LIMMA MODEL2 APPROACH(group-means parameterization) WITH ONE GROUP (Disease_state)

design <- model.matrix(~0 + Disease_state, data=targets_copy)
design

colnames(design) <- c("Control","HPS")
fit <- lmFit(betas, design)
head(fit$coefficients)

head(fit$stdev.unscaled*fit$sigma)

cont.matrix <- makeContrasts(HPSvsControl=HPS-Control, levels=design)
cont.matrix
#Contrasts
#Levels    HPSvsControl
#Control           -1
#HPS                1
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results = topTable(fit2, sort.by="P", adjust.method = "BH")
dim(results)
 

head(results)
results <- decideTests(fit2)
vennDiagram(results)

#UNUSED CODE START
#THIS CODE IS NOW GENERATING VECTOR . I TRIED getElement(P.value)results$P.value but not working
#Subsetting:
#res<-subset(results,results[['P.Value']]<.05)
#dim(res)
#[1] 10  6 
#head(res)
#STOP

#LIMMA MODEL1 APPROACH (treatment-contrasts parametrization)WITH ONE GROUP (Disease_state) AND ALL COVARIATES
#There's no need to form a contrast because the score variable is already in the model. 
#A continuous variable is its own contrast. 
#Explore covariates from Phenodata or targets_copy dataframe:
str(targets_copy)
#Disease_state=factor
#gender=character
#bni=int
#age=int
#fvc=int
#name=character (Since this is not factorable into meaningful levels as a categorical variable, 
#I will not treat "this "name" as covariate)
#Preparing covariates:
targets_copy$gender<- factor(targets_copy$gender, levels=c("female","male"))
targets_copy$bmi<-as.numeric(targets_copy$bmi)
targets_copy$fvc<-as.numeric(targets_copy$fvc)
targets_copy$age<-as.numeric(targets_copy$age)

design <- model.matrix(~Disease_state+gender+bmi+fvc+age, data=targets_copy)
design

fit <- lmFit(betas, design)
head(fit$coefficients)

head(fit$stdev.unscaled*fit$sigma)

fit2 <- eBayes(fit)
results = topTable(fit2, coef=ncol(design), sort.by="P",number = nrow(betas), adjust.method = "BH")
dim(results)

#SAME # AS BEFORE
head(results)

res<-subset(results,results$P.Value<.05)
dim(res)
##[1] 18770     6 
#THIS IS ~60% FEWER THAN ROWS FROM PREVIOUS APPROACH NOT ACCOUNTING FOR COVARIATES: [1] 40178     6
head(res)


#Visualization of Volcano plot
#making dataset
dat <- data.frame(foldchange = fit[["coefficients"]][,2], logPvalue =  -log10(fit2[["p.value"]][,2]))
dat$threshold <- as.factor(abs(dat$foldchange) < 0.4)


#GARBAGE COLLECTION NEEDED AGAIN error displayed after an "overuse" of ggplot. 

gc()


#Visualization
cols <- c("TRUE" = "grey", "FALSE" = "blue")
ggplot(data=dat, aes(x=foldchange, y = logPvalue, color=threshold)) +
  geom_point(alpha=.6, size=1.2) +
  scale_colour_manual(values = cols) +
  geom_vline(xintercept = 0.4, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = - 0.4, colour="#990000", linetype="dashed") +
  theme(legend.position="none") +
  xlab("Fold Change") +
  ylab("-log10 p value") +
  theme_bw() +
  theme(legend.position = "none")

results <- decideTests(fit2)
results
results <- decideTests(fit2)
#vennDiagram(results)
#Error in vennDiagram(results) : 
#  Can't plot Venn diagram for more than 5 sets


#LIMMA MODEL1 APPROACH (treatment-contrasts parametrization)WITH ONE GROUP (Disease_state) AND "Gender" COVARIATE
group<-as.factor(pData(gRatioSet.quantile)$Disease_state)
group<-relevel(group,ref="Control")
design <- model.matrix(~group+pData(gRatioSet.quantile)$predictedSex)
#design <- model.matrix(~pData(gRatioSet.quantile)$Disease_state) 
colnames(design)<-c("Intercept","Comp","Gender")
design


####generate adjusted beta matrix
fit <- lmFit(betas, design)
head(fit$coefficients)

head(fit$stdev.unscaled*fit$sigma)

fit2 <- eBayes(fit)
results <- topTable(fit2,coef="Comp",num=dim(fit2)[1],sort.by="P",adjust.method = "BH") ####Pval list 1
dim(results)
#[1] 742714      6
res<-subset(results,results$P.Value<.05)
dim(res)
#[1] 41338     6
#This is less than NO COVARIATES considered and more than with ALL COVARIATES considered

###############################################################DMP FINDING

#Confirm Order:
targets$Disease_state
colnames(betas)
targets$name
#ALL IN SAME ORDER

#FOR INITIAL DMP DEMONSTRATION, TO AVOID ERRORS, we will first use linear model for the effect of
#2-level Disease_state factor categorical variable (Control + HPS) with NO COVARIATES
dmp <- dmpFinder(betas, pheno = targets$Disease_state, type = "categorical")
dim(dmp)
#[1] 742714      4
head(dmp)
dmp<-subset(dmp,dmp$pval<.05)
dim(dmp)
#[1] 40891     4
dmp <- merge(dmp,ann850k[,c("UCSC_RefGene_Name","Relation_to_Island","UCSC_RefGene_Group")], by="row.names")
#This provides gene name, relation to island, and position in promoter or not
dmp<-as.data.frame(dmp)
head(dmp)
rownames(dmp)<-dmp$Row.names
#Remove column
dmp<-dmp[,-1]
head(dmp)


#We are particularly interested in methylation in gene promoter methylation as it is related to gene expression control. High methylation results in low gene expression
#Therefore, subsetting to get CpGs in promoters

## Outputs
proms<-dmp[grep("TSS1500|TSS200|5'UTR|1stExon",dmp$UCSC_RefGene_Group),]
dim(proms)

#This is the number of CpGs in promoters



#Specifically, we are further interested in number of CpGs in islands:
proms.island<-proms[grep("Island",proms$Relation_to_Island),]
#proms.island<-proms[proms$Relation_to_Island=="Island",]
dim(proms.island)

#This is the number of CpGs in promoters islands
head(proms.island)


proms.island<-merge(proms.island,betas,by="row.names")
rownames(proms.island)<-proms.island$Row.names
proms.island<-proms.island[,-1]




#We want to visually determine if methylation can be split among samples
write.csv(dmp,file="C.vs.H.csv",row.names=T)

#colnames(betas)<-targets$Disease_state
 
Control <- grep("C01C",colnames(betas),perl=T)
Control

HPS <- grep("C01H",colnames(betas),perl=T)
HPS
#samples<-ifelse(group=="KO","red","blue1")
spcol <- c(rep("blue1",length(Control)),rep("red",length(HPS)))
spcol
#[1] "blue1" "blue1" "blue1" "blue1" "blue1" "red"   "red"   "red"   "red"   "red"   "red"   "red" 


#my.data<-as.data.frame(betas[rownames(betas) %in% rownames(proms.island),])
l <- as.matrix(betas[rownames(betas) %in% rownames(proms.island),c(Control,HPS)])
head(l)
#               200514040133_R04C01C 200517480058_R03C01C 200517480058_R05C01C 200517480058_R06C01C
#cg16352085            0.2552397            0.2756281            0.2893039            0.3283099
#cg27009105            0.1681403            0.1080997            0.1350834            0.1269381
#cg04531633            0.3125911            0.2832128            0.2868397            0.2160284
#cg17368944            0.4528569            0.4023897            0.4464363            0.4461291
#cg14997037            0.7133778            0.7887480            0.7616112            0.7850129
#cg26339484            0.2085316            0.1712462            0.2021708            0.1568247

#Evidently, dataframe split with control first and HPS second

#####scale data
#data <- t(scale(t(my.data))) # z-score normalise each row (feature)
#data <- apply(data, 2, function(x) ifelse(x > 4, 4, x)) # compress data within range
#data <- apply(data, 2, function(x) ifelse(x < -4, -4, x)) # compress data within range
#spcol<-samples
#cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

#heatmap.2x(as.matrix(data[1:2000,]), col=rev(cols), Colv = T,Rowv=TRUE, scale="none", ColSideColors=spcol,
#           trace="none", dendrogram="none", 
#           cexRow=1, cexCol=.7,
#           main="WT vs KO Proms + Islands",
#           # labCol=NA,
#           labRow=NA, 
#          density.info="none",
#          hclust=function(x) hclust(x,method="complete"),
#           distfun=function(x) as.dist((1-cor(t(x)))/2),srtCol=45
           
#)


#pdf("heatmap.diffMeth.group.nocovariate.promsislands.pdf")
heatmap.2(l,
          main="Diffmeth CpG's ranked",
          labRow=NA,
          trace="none",
          na.rm=T,
          col=greenred,
          ColSideColors=spcol,
          distfun=function(x) dist(x,method="euclidean"),
          dendrogram = "column")


pca <- prcomp(t(na.omit(l)))
plot(pca$x[,1],pca$x[,2],col=spcol,cex=.2)
text(pca$x[,1],pca$x[,2],labels=rownames(pca$x),col=spcol,cex=1)

#dev.off()

#####testing candidates

#The following publication was cited and cross-referenceed to identify potential biomarkers:
#Mendoza N, Rivas E, Rodriguez-Roisin R, Garcia T, Bruguera M, Agusti A, et al. (2021) Liver
#epigenome changes in patients with hepatopulmonary syndrome: A pilot study. PLoS
#ONE 16(2): e0245046. https://doi.org/10.1371/journal.pone.0245046
#BAsed on this publication, we tested the following 10 exhibitng high logfold change in methylation between control and HPS
# ORMDL3|VAV3|MYT1|PAX8|CNTNAP2|HIVEP2|ANKRD28|FAM168A|SCAMP3|4 GMP

tested<-proms.island[grep("ORMDL3|VAV3|MYT1|PAX8|CNTNAP2|HIVEP2|ANKRD28|FAM168A|SCAMP3",proms.island$UCSC_RefGene_Name),]

Control <- grep("C01C",colnames(betas),perl=T)
Control

HPS <- grep("C01H",colnames(betas),perl=T)
HPS
#samples<-ifelse(group=="KO","red","blue1")
spcol <- c(rep("blue1",length(Control)),rep("red",length(HPS)))

l <- as.matrix(betas[rownames(betas) %in% rownames(tested),c(Control,HPS)])

heatmap.2(l,
          main="Diffmeth CpG's Biomarkers",
          labRow=NA,
          trace="none",
          na.rm=T,
          col=greenred,
          ColSideColors=spcol,
          distfun=function(x) dist(x,method="manhattan"),
          dendrogram = "column")

#####################WE WILL LATER REPEAT DIFFERENTIAL ANALYSIS FOR MODELLING Disease_state + ALL COVARIATES FOR ALL REGIONS (promoters, islands)

#Consider prom+islands diff methylated Cpgs using covariates adjustment, in 
#terms of Biological pathways, do you think they are related to the disease? 
#Specify which procedure you used it and also, please discuss properly.

##GO enrichment
#####GO enrichment analysis

#WE are first considering the "mapped" genes that are associated with the differentially methylated CpGs in promoter + island
#Note: One CpGs can map to several isoforms or different overlapping genes (MCL1, MCL1)
genesid<-proms.island$UCSC_RefGene_Name
head(genesid)

genesid<- strsplit(as.character(genesid),';')
dim(genesid)
genesid<-unique(unlist(genesid))
dim(genesid)
#NULL
head(genesid)

genesid<-genesid[ genesid != "" ]##remove empty elements from a vector or non-annotated genes
genesid <- genesid[!is.na(genesid)]
#####GO enrichment analysis
head(genesid)


#Need to revmap the gene symbols to NCBI GenBank ENTREZ gene reference accession numbers
#because hypergenetic tests requires processing numerical ENTREZ GeneIDs:

xx2=unlist(mget(as.character(genesid), ifnotfound=NA, revmap(org.Hs.egSYMBOL)))
head(xx2)

#All the genes in the human genome are in the "universe"
univ<-Lkeys(org.Hs.egGO)
#[ reached getOption("max.print") -- omitted 192382 entries ]=1000+192382


#ALTERNATIVELY:
eg<-bitr(genesid, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ego2 <- enrichGO(gene         = eg$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable=T
)

#Note: There are 3 separate ontologies: molecular function, cellular localization, biological process (BP)
param2<-new("GOHyperGParams",geneIds=xx2,universeGeneIds=univ,annotation="org.Hs.eg.db",ontology="BP",pvalueCutoff=0.01,conditional=FALSE,testDirection="over")
param2
#A GOHyperGParams instance
#category: GO 
#annotation: org.Hs.eg 

#Performing the test
hyp=hyperGTest(param2)
hyp

#There are evidently 11475 GEne Ontology biological Process pathways tested, with only 1374
#There are only 4025 genes in gne set size due to some genes not being annotated/NA, since there are more 
#Gene Symbols than Reference GEneIDs

#Get p-values for the test
gGhyp.pv<-pvalues(hyp)
gGhyp.odds<-oddsRatios(hyp)
gGhyp.counts<-geneCounts(hyp)

sigGO.ID<-names(gGhyp.pv[gGhyp.pv<0.01])

#Test the number of counts and move to new matrix
gGhyp.counts<-as.data.frame(gGhyp.counts)
gGhyp.counts$GOterms<-rownames(gGhyp.counts)
gGhyp.counts<-gGhyp.counts[rownames(gGhyp.counts)%in% sigGO.ID,]

#Here only show the significant GO terms of BP (Biological Processes)
#For other categories, like Molecular Function (MF) or Cellular Localization (CL) we follow similar procedure
sigGO.Term<-getGOTerm(sigGO.ID)[["BP"]]
 
#Combining into a total dataframe
results_GO<-cbind(as.data.frame(gGhyp.pv[gGhyp.pv<0.01]),as.data.frame(sigGO.Term),gGhyp.counts)
head(results_GO)

#Navigating the object more deeply at lower level: 

head(results_GO$sigGO.Term)

#Evidently, the gene ontology enrichment is related to many nitrogen-intensive pathways

write.csv(results_GO,"GO_enrichment_BP_CvsH")

#####DOSE
#Specifically Disease Ontology Enrichment

data(DO2EG)

#Check previous gene symbol values for differnetialy mythylaetd CpGs:
head(genesid)
 
eg.ids <- toTable(org.Hs.egALIAS2EG)

gene.set<-eg.ids[eg.ids$alias_symbol %in% genesid,"gene_id"]
#gene.set like universe is vector of NCBI ENtrez GEne IDS

univ <- Lkeys(org.Hs.egGO)
head(univ)
#[ reached getOption("max.print") -- omitted 192382 entries ]
#DO=Diesease Ontology

#NEED TO RE-VIST TO SELECT ONLY NOMINAL P-VALUE<0.05 (NOT ALL AS = 1) WITH NO FDR, ADJUSTED q-VLUE 
#FOR NOW, I AM USING FDR:
x <- enrichDO(gene          = gene.set,
              ont           = "DO", 
              pvalueCutoff  = 1,
              pAdjustMethod = "fdr",
              universe      = univ,
              qvalueCutoff  = .05,
              readable      = FALSE)

results.dose<-as.data.frame(x)

results.dose
#[1] ID          Description GeneRatio   BgRatio     pvalue      p.adjust    qvalue      geneID      Count      
#<0 rows> (or 0-length row.names)
#TOO STRINGENT CUTOFF WITH FDR LEADS TO 0 results
#REDO WITh p-value cutoff < 0.1:
x <- enrichDO(gene= gene.set,ont= "DO", pvalueCutoff  = 0.1, universe = univ,  readable= FALSE)
results.dose<-as.data.frame(x)

dotplot(x)
#The red color are more significant over-representation (hypermethylated) compared to blue (hypomethylated) and the larger circles
#a represent higher number of genes 

############################################KEGG
ego2 <-enrichKEGG(
  eg$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  #universe,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)
dotplot(ego2, title=" Prom islands Cpgs ",showCategory=25)

#############################################

#missmethyl
#enrichment_GO <- gometh(rownames(proms.island),all.cpg = rownames(betas),collection = "GO", 
                        #array.type = "450K",plot.bias = T,prior.prob = T,equiv.cpg = T,anno = ann450k) 

#enrichment_GO<-enrichment_GO[enrichment_GO$ONTOLOGY=="BP",]

#enrichment_GO<-enrichment_GO[enrichment_GO$FDR<.05,]

#enrichment_KEGG <- gometh(rownames(proms.island),all.cpg = rownames(betas),collection = "KEGG", 
                          #array.type = "450K",plot.bias = T,prior.prob = T,equiv.cpg = T,anno = ann450k) 

#enrichment_KEGG<-enrichment_KEGG[enrichment_KEGG$FDR<.05,]
#enrichment_KEGG[enrichment_KEGG$P.DE<.05,]


###########################################DMR
#NOW EXPLORING IDENTIFICATION OF BUMPS (DIFFERENTIALLY METHYLATED REGIONS):
#We avoid using package bumphunter as it is prohibitively memory-intensive
#####bumphunter
#library(bumphunter)
#pheno <- pData(gRatioSet.quantile)$Status
#designMatrix <- model.matrix(~ pheno)
#Run the algorithm with a large number of permutations, say B=1000:

#dmrs <- bumphunter(gRatioSet.quantile, design = designMatrix, cutoff = 0.01, B=1000, type="Beta")

#library(doParallel)
#registerDoParallel(cores = 3)
#The results of bumphunter are stored in a data frame with the rows being the different differentially methylated regions (DMRs):

#GFRA1 and GSTM2,SEPT9

#names(dmrs)
#head(dmrs$table, n=3)

#################################################DMRCate 
#WE USE DMRCate for for DMR due to its less memory-intensive algorithms

#GFRA1 and GSTM2,SEPT9
dmp[grep("GFRA1",dmp$UCSC_RefGene_Name),]

data(dmrcatedata)

myMs<-getM(gRatioSet.quantile)
#myMs <- logit2(myBetas)###convert back the Betas to M
myBetas<-getBeta(gRatioSet.quantile)

#WE WILL USE M VALUES HERE TO AVOID DISPERSION:
#CONFIRM SAME ORDER:
colnames(myMs)

targets$name
targets$Basename
colnames(myMs)<-targets$name

#myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05,rmXY = T,rmcrosshyb = TRUE) ### get out SNPs of my samples
#myBetas <- rmSNPandCH(myBetas, dist=2, mafcut=0.05) ### get out SNPs of my samples

myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05) ### get out SNPs of my samples
myBetas <- rmSNPandCH(myBetas, dist=2, mafcut=0.05) ### get out SNPs of my samples
#Probe IDs from EPICv1 or earlier detected. Proceeding...
#Cannot connect to ExperimentHub server, using 'localHub=TRUE' instead
#Using 'localHub=TRUE'
#If offline, please also see BiocManager vignette section on offline use
#snapshotDate(): 2024-05-28
#see ?DMRcatedata and browseVignettes('DMRcatedata') for documentation
#loading from cache

myMs.noSNPs
#[ reached getOption("max.print") -- omitted 742631 rows ]
###establish the groups
Control <- grep("C01C",colnames(myMs.noSNPs),perl=T)
Control

HPS <- grep("C01H",colnames(myMs.noSNPs),perl=T)
HPS

#AS factor????
pheno <- as.factor(pData(gRatioSet.quantile)$Disease_state)
pheno<-relevel(pheno,ref="Control")
designMatrix <- model.matrix(~pheno)
designMatrix

#Please specify either 'EPICv2' or 'EPICv1' for arraytype. EPICv2 probe IDs have 15 characters, e.g. cg00000029_TC21. EPICv1 probe IDs have 10 characters, e.g. cg00000029.
myannotation <- cpg.annotate("array", myMs.noSNPs, what="M", arraytype = "EPICv1",
                             analysis.type="differential", design=designMatrix, coef=2)

#Your contrast returned no individually significant probes. Try increasing the fdr. 
#Alternatively, set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.
myannotation <- cpg.annotate("array", myMs.noSNPs, what="M", arraytype = "EPICv1",
                             analysis.type="differential", design=designMatrix, fdr=0.9,coef=2)
#Your contrast returned 212666 individually significant probes. We recommend the default setting of pcutoff in dmrcate().
class(myannotation)
#[1] "CpGannotated"
#attr(,"package")
#[1] "DMRcate"

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
#my.dmrs<-dmrcoutput$results
#Error in dmrcoutput$results : $ operator not defined for this S4 class
#my.dmrs<-my.dmrs[my.dmrs$Stouffer<.05,]
#my.dmrs<-my.dmrs[my.dmrs$no.cpgs >= 5,]

#results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
#Cannot connect to ExperimentHub server, using 'localHub=TRUE' instead
#Using 'localHub=TRUE'

#dmrs.res<-as.data.frame(results.ranges)
#Error in h(simpleError(msg, call)) : 
#head(my.dmrs)
#dim(my.dmrs)

###filter those with less than 5 CpGs inside DMR
#my.dmrs<-subset(dmrs.res,dmrs.res$no.cpgs >= 5)
#my.dmrs<-subset(my.dmrs,my.dmrs$Fisher <.05)
#write.csv(my.dmrs,"dmrs.proms.csv")

###select specifical
#head(dmr.res)
#dmrs.res[grep("VAV3",my.dmrs.res$overlapping.genes),]
#dmrs.res[grep("VAV3",my.dmrs.res$overlapping.promoters),]
#3 DMRs overlapping promoters
 
#groups <- c(HPS="magenta", Control="forestgreen")
#type<-pheno
#cols <- groups[as.character(type)]

#DMR.plot(ranges=results.ranges, dmr=14269, CpGs=ilogit2(myMs), phen.col=cols, genome="hg19")

#DMR.plot(ranges=results.ranges, dmr=1, CpGs=myMs.noSNPs,what="M",
         #arraytype = "450K", phen.col=cols, genome="hg19")

#DMR.plot(ranges=results.ranges, dmr=598, CpGs=myBetas, phen.col=cols, genome="hg19")
#dev.off()
#DMR.plot(ranges=results.ranges, dmr=598, CpGs=myBetas, what="Beta", arraytype = "450K",
#         phen.col=cols, genome="hg19")

################PLOTS

####DMRff

#library(dmrff)
 
#group<-as.factor(pData(gRatioSet.quantile)$Status)
#group<-relevel(group,ref="Normal")

#design <- model.matrix(~group)##Cancer up
#colnames(design)<-c("Intercept","Comp")
####generate adjusted beta matrix
#fit <- lmFit(betas, design)

#fit2 <- eBayes(fit)

#stats <- data.frame(estimate=fit2$coefficients[,"Comp"],se=sqrt(fit2$s2.post) * fit2$stdev.unscaled[,"Comp"],p.value=fit2$p.value[,"Comp"])

####add annotation to stats
#stats<-merge(stats,ann450k[,c("chr","pos")],by="row.names")

#dmrs <- dmrff(estimate=stats$estimate,
#              se=stats$se,
#              p.value=stats$p.value,
#              methylation=betas,
#              chr=stats$chr,
#              pos=stats$pos,
#              maxgap=500,
#              verbose=T)


#dmrs <- dmrs[which(dmrs$p.value < 0.05 & dmrs$n > 1),]
#annot<-dmrff.sites(dmrs, stats$chr, stats$pos)


#dmrs<-makeGRangesFromDataFrame(dmrs)

##########################################################################################expression

#expres<-read.csv("../Expression_data.csv")

#genesid<-proms.island$UCSC_RefGene_Name

###split

#genesid<- strsplit(as.character(genesid),';')
#genesid<-unique(unlist(genesid))

#genesid<-genesid[ genesid != "" ]##remove empty elements from a vector
#genesid <- genesid[!is.na(genesid)]

###select the genes of the proms+island
#expres<-expres[expres$Gene %in% genesid,]

##mean
#expres<-aggregate(.~Gene, data=expres[,c(2:13,15)], mean)
#rownames(expres)<-expres$Gene
#expres<-expres[,-1]

#Normal <- grep("\\d+P",colnames(expres),perl=T)
#Cancer <- grep("\\d+T",colnames(expres),perl=T)

#pdf("ANOVA_PLOT_expression_signature_proms.cgi.pdf")
#GFRA1 and GSTM2,SEPT9

#i<-"SEPT9"

#for (i in 1:(dim(expres)[1])){
  
#  Norm<-as.numeric(expres[i,Normal])
#  Canc<-as.numeric(expres[i,Cancer])
  
#  class<-c(rep("Norm",length(Norm)), rep("Canc",length(Canc)))
#  y<-c(Norm,Canc)
  
#  result<-aov(y~as.factor(class))
  
#  res<-summary(result)
  
#  pval<-res[[1]][["Pr(>F)"]][1]
#  pval<-signif(pval, digits=3)
  
#  legend=c("Normal","Cancer")
  
#  boxplot(Norm,Canc,names=legend, col="lightgrey",main=paste(rownames(expres)[i],"p-value",pval, sep=" "), ylab="Expression")
#}

#dev.off()

###methylation

##add beta values to proms+islands results

#proms.island<-merge(proms.island,betas,by="row.names")
#rownames(proms.island)<-proms.island$Row.names
#proms.island<-proms.island[,-1]

#Normal <- grep("\\d+P",colnames(proms.island),perl=T)
#Cancer <- grep("\\d+T",colnames(proms.island),perl=T)

#pdf("ANOVA_PLOT_methylation_signature_proms.cgi.pdf")

#for (i in 1:(dim(proms.island)[1])){
  
#  Norm<-as.numeric(proms.island[i,Normal])
#  Canc<-as.numeric(proms.island[i,Cancer])
  
#  class<-c(rep("Norm",length(Norm)), rep("Canc",length(Canc)))
#  y<-c(Norm,Canc)
  
#  result<-aov(y~as.factor(class))
  
#  res<-summary(result)
  
#  pval<-res[[1]][["Pr(>F)"]][1]
#  pval<-signif(pval, digits=3)
  
#  legend=c("Normal","Cancer")
  
#  boxplot(Norm,Canc,names=legend, col="lightgrey",main=paste(rownames(proms.island)[i],proms.island$UCSC_RefGene_Name[i],"p-value",pval, sep=" "), ylab="Methylation")
#}
#dev.off()

#####################################################################################

#The exercise will be structured in blocks:
#BLOCK1: Quality Control
#1. How many samples passed the quality control? Can you specify which is the 
#procedure you followed? Our overall procedure involved conversion as dollows:
#IDAT Files->PGChannelSet->MethylSet->RAtioSet (and -> Genomic MethylSet in parallel via Ratioconvert)->Genomic RatioSet
#Based on the bimodal distribution for diseased_state suggesting good hybridizations
#all 12 samples passed initial QC.
#2. How many probes pass the quality control? Can you specify which is the 
#procedure you followed? Based on above-the-dashed-line threshhold plot, bisulfite II plot for red/green fluoescence,
#the majority of 866836 probes with detP in rgset survived QC. However, 837623 remained after detP <0.001 test in rgset, 
#835034 remained in GRatioset, 816112 remained after ommiting sexual imprinting or sexual methylation,742714 remained after
#filtering out cross-reactive probes
#3. Can you predict the sex of the samples? Do they correlate with sex in the 
#phenodata? Yes, I was able to predict the sex of the 12 samples and they correlated with sex in phenodata.
#BLOCK2: Analysis (in order to standardize, normalize by Quantiles). Do not expect 
#huge differences in this case, hence use nominal pvalue to establish differential 
#analysis control instead of FDR.
#1. How many CpgS are differentially methylated without considering covariates 
#(age, sex…)? That means all, promoter, and promoter plus island.
#Using first limma, 40178 CpGs were differentially methylated with NO COVARIATES
#Using second DMP, 40891 for all , 14529 for promoter and 5865 for promoter.island were 
#deemed differentially methylated
#2. How many CpGs are differentially methylated considering covariates 
#(look in the Phenodata)? That means all, promoter, and promoter plus island 
#(Hint: use limma design matrix to assign covariates)
#Using limma, 18770 were differentially methylated for ALL COVARIATES and 41338
#were differentially mehylated when accunting for ONLY GENDER COVARIATE

#Question: Which do you think is the best method, considering covariates or 
#not, that give you the best approach to get differential methylated CpGs that 
#explain the disease. Discuss it properly

#As evidenced by limma model, accounting for ALL COVARIATES resulted in a more conservative and lower
#number of differentially methylated CpGs, suggesting multi-factor interactions. Interestingly,
#accounting for ONLY GENDER COVARIATE resulted in higher number of differentially methylated CpGs.
#Accounting for covariates is recommended to elucidate such effects, even though PCA analysis 
#when accountign for NO VARIATES clearly delineated clustering segregation around H-DISEASE vs. CONTROL
#3. Consider prom+islands diff methylated Cpgs using covariates adjustment, in 
#terms of Biological pathways, do you think they are related to the disease? 
#Specify which procedure you used it and also, please discuss properly.
#A HyperG-Test Gene to GO BP (Biological Processes) test for over-representation was performed
#resulting in 11475 GO BP ids tested (1374 have p < 0.01)
#The Selected gene set size: 4025, the Gene universe size: 18888, and the Annotation package: org.Hs.eg 
#There are evidently 11475 GEne Ontology biological Process pathways tested, with only 1374 having p < 0.01 significance,
#including nitrogen-containing nucleic acid metabolism.
#There are only 4025 genes in gne set size due to some genes not being annotated/NA, since there are more 
#Gene Symbols than Reference GEneIDs. There is apparent 
#relation between disease and prom+islands diff methylated Cpgs using covariates adjustment.

#4. Consider to establish a threshod using phenodata FVC column ONLY in HPS 
#patients: i-e less than 90 and more than 90 value. Forced vital 
#capacity (FVC) is the volume of air that can forcibly be blown out after full 
#inspiration, measured in liters. It is an important issue valoring the disease, 
#Do you find dmps in promoter + islands (remember, work with nominal pval)? 
#Are they good candidates for biomarkers? (heatmaps, GO enrichment,..)
#compared to the previous analysis?
  #DMPs in prmoter+islands bsed on non-adjusted nominal p-value were discovered and basd on 
#Based on our depicted CpGs biomarkers heatmap generated, contrsting red vs. green scaled coloration
#suggests that they are good candidates fro biomarkers of the disease.








