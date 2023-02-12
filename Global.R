##### Install required packages for S:CORT data analysis pipeline

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

if (!require("affy", quietly = TRUE))
  BiocManager::install("affy")
library(affy)

if (!require("limma", quietly = TRUE))
  BiocManager::install("limma")
library(limma)

if (!require("samr", quietly = TRUE))
  BiocManager::install("samr")
library(samr)

if (!require("genefilter", quietly = TRUE))
  BiocManager::install("genefilter")
library(genefilter)

if (!require("ggfortify", quietly = TRUE))
  BiocManager::install("ggfortify")
library(ggfortify)

if (!requireNamespace("sva", quietly = TRUE))
  BiocManager::install("sva")
library(sva)

if (!requireNamespace("MLeval", quietly = TRUE)){
  BiocManager::install("MLeval")}
library(MLeval)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)

if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
library(tidyr)

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)

if (!require("e1071", quietly = TRUE))
  install.packages("e1071")
library(e1071)

if (!require("caret", quietly = TRUE))
  install.packages("caret")
library(caret)

if (!require("glmnet", quietly = TRUE))
  install.packages("glmnet")
library(glmnet)

if (!require("randomForest", quietly = TRUE))
  install.packages("randomForest")
library(randomForest)

if (!require("Boruta", quietly = TRUE))
  install.packages("Boruta")
library(Boruta)

if (!require("pROC", quietly = TRUE))
  install.packages("pROC")
library(pROC)

if (!require("mlbench", quietly = TRUE))
  install.packages("mlbench")
library(mlbench)

if (!require("doMC", quietly = TRUE))
  install.packages("doMC")
library(doMC)
registerDoMC(cores = 8)