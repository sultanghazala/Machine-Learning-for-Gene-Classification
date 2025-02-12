#=============================================================#
# R Script for ML Model Development for Gene Expression Data  #
# ------------------------------------------------------------#
# Developer: Dr. GhazalaSultan | Dept. Comp. Sci., AMU, INDIA #
#=============================================================#

# Load the required packages
install.packages("xlsx")
install.packages("caret")      # machine learning library package
install.packages("glmnet")     # computing penalized regression - elasticNet
install.packages("caTools")    # For Logistic regression 
install.packages('pROC')       # For ROC curve to evaluate model 

# Importing required libraries
library(readxl)  
library(tidyverse)
library(caret)  
library(glmnet) 
library(caTools)
library(pROC)  
library(ggplot2)


# Prerequisites 
graphics.off()  # clear all graphs from RStudio
rm(list = ls()) # remove all files from your workspace
# Ctrl+L -> Clear Console


#======================#
## Load Data
#======================#
#df <- read_excel("exprTable.xlsx", sheet = 1)
df <- read_excel("shortTable.xlsx", sheet = 1)


#===============================#
# Data Pre-processing |   EDA  
#===============================#
# Check the structure of the dataset
head(df)
str(df)
dim(df)

# remove rows that contain NA values
df <- df[complete.cases(df), ]
head(df)
dim(df)

#Calculate Mean of duplicate genes
x <- df
x <- data.frame(x)
x <- do.call(rbind,lapply(lapply(split(x,x$Symbols),`[`,2:ncol(x)),colMeans))
dim(x)

#Convert rownames as a 1st column with header Symbols -> which became rownames after previous operation
library(tibble) # from tidyverse
x <- data.frame(x)
x <- tibble::rownames_to_column(x, var="Symbols")
head(x)
dim(x) 

df <- x
# Transpose table 
library(sjmisc)  #install.packages("sjmisc")
df_t <- rotate_df(df, cn=T)
Symbols <- colnames(df[-1])
df_t <- cbind(Symbols, df_t)
write.csv(df_t, "transposed_table.csv", row.names=F)

df_t <- read.csv("transposed_table.csv", h=T)
dim(df_t)
df_t[1]

# N_GSM398104.CEL to N/T
df_t[,1] <- gsub("_.*$", "", df_t[,1])
df_t[1]

# convert N/T from char to factor
df_t[1] <- factor(df_t$Symbols)
str(df_t)

# df_t -> df
df <- df_t # df is df_NT

# view transformed data
str(df)
df <- data.frame(df)


#======================================#
# Visualize Dataset - Figures - Plots 
#======================================#
#####  Box and Whisker Plots  ##### 
# Given that the input variables are numeric, we can create box and whisker plots of each
png("box_and_whisker_plots.png")
par(mfrow=c(1,4))
for(i in 2:9) {
  boxplot(df[,i], main=names(df)[i], col="blue")
}
dev.off()


#####  scatterplot matrix  ##### 
library(ggplot2)
# split input and output
x <- df[,2:ncol(df)]  # x -  inputs attributes
y <- df[,1]     # y -  outputs attributes
y<- as.factor(y)
featurePlot(x=x, y=y, plot="ellipse") # time intensive

# view number of samples/replicates in each condition
plot(y, col="blue")



#=======================#
#   Data Splitting
#=======================#
# create 60%/40% for training and testing dataset
set.seed(101)
split <- createDataPartition(df$Symbols, p=0.60, list=FALSE)
train <- df[split,]
test <- df[-split,]

# view dimensions of dataset, training set, test set
dim(df)
dim(train)
dim(test)


#set cross-validation control for training and evaluation metric
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"


#============================#
#   Build ML Models   
#============================#
# kNN (k Nearest Neighbour)
#----------------------------#
library(ISLR)
library(cowplot)
set.seed(101)
fit.knn <- train(Symbols~., data=train, 
                 method="knn", trControl=control) # using default values of k in kNN
fit.knn

# Feature/Gene Importance
varImp(fit.knn)
plot(varImp(fit.knn, scale=FALSE), top=20) 
# export plot  
tiff("kNN_varimp.tiff")
plot(varImp(fit.knn, main="k-Nearest Neighbors", scale=TRUE), top=30)
dev.off()

# Make Predictions # Confusion Matrix
predict.knn <- predict(fit.knn, newdata = test)
cfm.knn <- confusionMatrix(predict.knn, as.factor(test$Symbols))
cfm.knn

# Plot Confusion Matrix
library(ggplot2)
library(dplyr)
# get confusion metrics table
cfm.knn_tb <- confusionMatrix(predict.knn, test$Symbols, positive = "T")$table
table <- data.frame(cfm.knn_tb)
plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "high", "low")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

ggplot(data = plotTable, mapping = aes(x = Reference, 
                                       y = Prediction, fill = Freq, alpha = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_gradientn(name = "Freq", colors = c("#E0E0E0", "#009194"),
                      values = scales::rescale(c(0, 2, 5))) +
  guides(alpha = "none") +  
  theme_bw() +
  xlim(rev(levels(plotTable$Reference)))+
  coord_fixed(ratio = 1)


#------------------------------#
# SVM (Support Vector Machine)
#------------------------------#
library(kernlab)
set.seed(101)
fit.svm <- train(Symbols~., data=train, trControl=control)
fit.svm

# Feature Importance
varImp(fit.svm)
plot(varImp(fit.svm, scale=FALSE), top=20)
# export plot
png("svm_varimp.png")
plot(varImp(fit.svm, main="Support Vector Machines", scale=TRUE), top=30)
dev.off()

# Make Predictions # Confusion Matrix
predict.svm <- predict(fit.svm, newdata = test)
cfm.svm <- confusionMatrix(predict.svm, as.factor(test$Symbols))
cfm.svm

# Plot Confusion Matrix
cfm.svm_tb <- confusionMatrix(predict.svm, test$Symbols, positive = "T")$table
table <- data.frame(cfm.svm_tb)
plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "high", "low")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))
ggplot(data = plotTable, mapping = aes(x = Reference, 
                                       y = Prediction, fill = Freq, alpha = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_gradientn(name = "Freq", colors = c("#E0E0E0", "#009194"),
                       values = scales::rescale(c(0, 2, 5))) +
  guides(alpha = "none") +  
  theme_bw() +
  xlim(rev(levels(plotTable$Reference)))+
  coord_fixed(ratio = 1)


#------------------------------#
# Logistic regression # E-net
#------------------------------#
set.seed(101)
fit.enet <- train(Symbols ~ ., data = train, method ='glmnet',
                  type.measure="deviation", family="binomial",
                  tuneGrid = expand.grid(alpha = seq(0,1,length=10), lambda = seq(0.0001,0.2,length=5)),
                  trControl = control)
print(fit.enet)    

varImp(fit.enet)
plot(varImp(fit.enet), top=20, main="ENET") #, col="red")

# Make Predictions & Model Evaluation
predict.enet <- predict(fit.enet, newdata = test)
cfm.enet <- confusionMatrix(predict.enet, as.factor(test$Symbols))
cfm.enet


# Plot Confusion Matrix
cfm.enet_tb <- confusionMatrix(predict.enet, test$Symbols, positive = "T")$table
table <- data.frame(cfm.enet_tb)
plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "high", "low")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))
ggplot(data = plotTable, mapping = aes(x = Reference, 
                                       y = Prediction, fill = Freq, alpha = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_gradientn(name = "Freq", colors = c("#E0E0E0", "#009194"),
                       values = scales::rescale(c(0, 2, 5))) +
  guides(alpha = "none") +  
  theme_bw() +
  xlim(rev(levels(plotTable$Reference)))+
  coord_fixed(ratio = 1)