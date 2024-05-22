if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.12")

if (!requireNamespace("dplyr", quietly = TRUE)){
  BiocManager::install("dplyr")}
library(dplyr)

if (!requireNamespace("ggplot2", quietly = TRUE)){
  BiocManager::install("ggplot2")}
library(ggplot2)

if (!requireNamespace("shinycssloaders", quietly = TRUE)){
  BiocManager::install("shinycssloaders")}
library(shinycssloaders)

if (!requireNamespace("doMC", quietly = TRUE)){
  BiocManager::install("doMC")}
require(doMC)
registerDoMC(cores=6)
detectCores(all.tests = FALSE, logical = TRUE)

if (!requireNamespace("data.table", quietly = TRUE)){
  BiocManager::install("data.table")}
library(data.table)

if (!requireNamespace("parallel", quietly = TRUE)){
  BiocManager::install("parallel")}
library(parallel)
numCores <- detectCores()
numCores

if (!requireNamespace("tidyverse", quietly = TRUE)){
  BiocManager::install("tidyverse")}
library(tidyverse)

if (!requireNamespace("limma", quietly = TRUE)){
  BiocManager::install("limma")}
library(limma)

if (!requireNamespace("Biobase", quietly = TRUE)){
  BiocManager::install("Biobase")}
library(Biobase)

if (!requireNamespace("raster", quietly = TRUE)){
  BiocManager::install("raster")}
library(raster)

if(!requireNamespace("plotROC", quietly = TRUE)){
  BiocManager::install("plotROC")}
library(plotROC)

if (!requireNamespace("sva", quietly = TRUE)){
  BiocManager::install("sva")}
library(sva)

if (!requireNamespace("e1071", quietly = TRUE)){
  BiocManager::install("e1071")}
library(e1071)

if (!requireNamespace("caret", quietly = TRUE)){
  BiocManager::install("caret")}
library(caret)

if (!requireNamespace("gbm", quietly = TRUE)){
  BiocManager::install("gbm")}
library(gbm)

if (!requireNamespace("Boruta", quietly = TRUE)){
  BiocManager::install("Boruta")}
library(Boruta)

path_to_discovery_data <- "Data/Discovery/"
path_to_validation_data <- "Data/Validation/"


