#=============================================================#
# R Script for Gene Expression Analysis (with ML/AI)          #
# ------------------------------------------------------------#
# Developer: Dr. GhazalaSultan | Dept. Comp. Sci., AMU, INDIA #
#=============================================================#

# DEG Analysis - Part 1
#-----------------------#
# 1. Reading CEL files
# 2. Expression Normalization (RMA)
# 3. Pre and Post Normalization Expression Visualization


# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("limma")

# Load packages
library(limma)
library(hgu133plus2cdf)
library(affy)
library(IRanges)
library(RColorBrewer)
library(affyPLM)


targets <- readTargets("target_H_DCIS.txt")
targets

#Read CEL Files
data <- ReadAffy(filenames = targets$FileName)
data


#=====================#
# RMA Normalization
#=====================#
eset <- rma(data)
normset <- exprs(eset)
pData(eset)
write.csv(normset, "ExpSet_PostNorm.csv", quote = F)


#========================#
# Probe/Gene Annotation
#========================#
BiocManager::install("hgu133plus2.db")
library("hgu133plus2.db")

# Match probe IDs and retrieve Gene SYMBOLs
probes=row.names(normset)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))

# Combine gene annotations with raw data
normset_anno = cbind(probes,Symbols,normset)
write.csv(normset_anno, "ExpSet_PostNorm_Annotated.csv", quote = F, row.names = F)



#=====================#
#   Box Plot
#=====================#
par(mfrow=c(1,2))

#Boxplot Before Normalization
#colors = c(rep("red",5),rep("blue",9))
tiff(file="Boxplot_Pre-Normalization.tiff", bg="transparent", width=400, height=500)
par(mar = c(12, 4, 6, 2) + 0.1); # This sets the plot margins
boxplot(data,col="red", main="Boxplot Pre-Normalization", las=2, cex.axis=0.74, ylab="Intensities" )
title(xlab = "Sample Array", line = 8); # Add x axis title
#xlab="Sample Array")
dev.off()

#Boxplot After Normalization
tiff(file="Boxplot_Post-Normalization.tiff", bg="transparent", width=400, height=500)
par(mar = c(12, 4, 6, 2) + 0.1); # This sets the plot margins
boxplot(normset,col="blue", main="Boxplot Post-Normalization", las=2, cex.axis=0.74, ylab="Intensities") #, col=colors 
title(xlab = "Sample Array", line = 8); # Add x axis title
dev.off()

