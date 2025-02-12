#=============================================================#
# R Script for Gene Expression Analysis (with ML/AI)          #
# ------------------------------------------------------------#
# Developer: Dr. GhazalaSultan | Dept. Comp. Sci., AMU, INDIA #
#=============================================================#


# DEG Analysis - Part 2
#-----------------------#
# 1. Data Cleaning - Dimensionality Reduction (PCA)
# 2. Differentially Expressed Genes (DEGs) Identification 
# 3. DEGs Visualization (Volcano plot & Heatmap)


#=====================#
#   PCA Plot
#=====================#
#install.packages("factoextra")
library(tidyverse)
library(factoextra)

data <- read.csv("ExpSet_PostNorm.csv")

nrow(data)  # Check the number of rows in your data
length(c(rep("Healthy", 5), rep("DCIS", 9)))  # Check the length of your group labels

# Adjust Group Labels
data_t <- t(data[,-1 ])  # Exclude the first col (gene names) during transpose
data_t <- as.data.frame(data_t)
data_t$Group <- c(rep("Healthy", 5), rep("DCIS", 9))

# Perform PCA
pca_res <- prcomp(data_t[, -ncol(data_t)])

# view all PC scores 
head(pca_res$x)

png("pca.png")
# Generate PCA plot
fviz_pca_ind(pca_res,
             geom.ind = c("point", "text"),
             col.ind = data_t$Group,
             palette = c("blue", "red"),
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Group",
             labelsize = 2
)
dev.off()


#=======================#
#   DEG Identification
#=======================#
# Model Matrix Design
# Let's create a model matrix using the factor() function to represent the condition labels ("Healthy" or "DCIS")
design <- model.matrix(~factor(c("Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "DCIS", "DCIS", "DCIS", "DCIS", "DCIS", "DCIS", "DCIS", "DCIS", "DCIS")))
# Now assign names ("Healthy" and "DCIS") to the columns of the model matrix
colnames(design) <- c("Healthy","DCIS")


# Fits a linear model for each gene based on the given series of arrays
# It estimates the relationship between gene expression and conditions
fit <- lmFit(normset[,1:ncol(normset)], design) 


# Contrast Matrix Design
# Define the specific comparison between conditions you want to analyze.
cont.matrix = makeContrasts(DCIS-Healthy, levels=design)

# Fitting model with Contrasts(2 Groups), so apply the defined contrast to the previously fitted model (fit)
fit2 <- contrasts.fit(fit, cont.matrix)

# Model optimization / Empirical Bayes Moderation
# Improves the estimation of variances for genes with low expression
# Computes moderated t-statistics and log-odds (B-stats) of differential expression by empirical Bayes shrinkage of the standard errors towards a common value
# contrast-specific information from fit2 is incorporated into the fit object, which is then passed to eBayes()
fit2 <- eBayes(fit)  


# Result Top Table
topTable(fit2, coef = 2, adjust.method = "BH") 
x <- topTable(fit2, coef=2, adjust="BH", sort.by="logFC", number=100000);
write.csv(x, "Result_Table_logFCsorted.csv", quote = F)

# Filter & Save DEGs
logFC_1 <- x[x$P.Value < 0.05 & (x$logFC > 1 | x$logFC < -1), ]
write.csv(logFC_1,"pval_0.05_fc_1.csv", quote = F)



#==========================#
# AnnotatE filtered DEGs
#==========================#
BiocManager::install("hgu133plus2.db")
library("hgu133plus2.db")

probes=row.names(x)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))

# Combine gene annotations with raw data
deg_anno = cbind(probes,Symbols, x)
write.csv(deg_anno, "DEGs_Annotated.csv", quote = F, row.names = F)



#================================================#
# DEG Viz (Volcano plot)
#================================================#
#install.packages("gdata")
#install.packages("gplots")
library(gdata)
library(gplots)

png(filename = "VolcanoPlot_FC_1.png")
with(x, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot"))
with(subset(x, P.Value < 0.05 & logFC > 1 ), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(x, P.Value < 0.05 & logFC < -1), points(logFC, -log10(P.Value), pch=20, col="green"))
dev.off()



#====================================#
# DEG Viz (Heatmap DEG Expression)
#====================================#
data <- read.csv(file = "heatmap_expdata.csv", h=T)
rnames <- data[,1]
mat_data <- data.matrix(data[,2:ncol(data)])
rownames(mat_data) <- rnames
my_palette <- colorRampPalette(c("red", "blue"))(n = 299)
col_breaks = c(seq(0,5,length=100), # for red
               seq(5.1,10,length=100),  # for yellow
               seq(10.1,15,length=100)) # for green
tiff("heatmap_exp_deg_cluster.tiff",     
     width = 6*300,        # 5 x 300 pixels
     height = 6*300,
     res = 300,            # 300 pixels per inch
     pointsize = 8)        # smaller font size

heatmap.2(mat_data,
          main = "Heatmap", # heat map title
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks, 
          dendrogram="both",     # only draw a row dendrogram
          Colv="T" ,         # turn off column clustering
          lhei = c(1,7)         # Key size width adjustment
)            
dev.off()

