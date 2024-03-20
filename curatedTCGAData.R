if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("curatedTCGAData")
#BiocManager::install("TCGAutils")

library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)

#browseVignettes("curatedTCGAData")

# check here for 'diseaseCode' and 'assays' inputs
#?curatedTCGAData

#curatedTCGAData("BRCA", version = "2.0.1")

# "Methylation_methyl450","CNVSNP"
brca <- suppressMessages(curatedTCGAData("BRCA", c("RNASeq2GeneNorm",
                                                   "miRNASeqGene"), version = "2.0.1", dry.run = FALSE))

# Filter primary tumor cases (code 01)
tumor <- TCGAprimaryTumors(brca)

# Check for replicates
anyReplicated(tumor)
#   BRCA_miRNASeqGene-20160128 BRCA_RNASeq2GeneNorm-20160128 
#   FALSE                         FALSE 

cd = colData(brca)[which(!is.na(brca$OS.Time)),]
library(survival)
osurv = Surv(cd$OS.Time/365.25, cd$OS.event)
# Plot overall survival curve
plot(survfit(osurv~1), main = "Overall BRCA Survival", xlab = "Years")

# Some cancer datasets contain associated subtype info
# within the clinical datasets
head(getSubtypeMap(brca))

summary(brca)

# common clinical variables
head(getClinicalNames("BRCA"))

colData(brca)[, getClinicalNames("BRCA")][1:5, 1:5]

sampleTables(brca)

# 用来解读上面信息的code
data(sampleTypes, package = "TCGAutils")
sampleTypes

#primaryTumors <- TCGAprimaryTumors(brca)

#sampleTables(primaryTumors)

# Set the working directory to disk D:
setwd("E:/1上课/Capstone")
getwd()

td <- tempdir()
tempd <- file.path(getwd(), "BRCA")
if (!dir.exists(tempd))
  dir.create(tempd)

exportClass(tumor, dir = tempd, fmt = "csv", ext = ".csv")
