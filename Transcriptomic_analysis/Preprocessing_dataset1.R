# Processing Dataset GSE157643

## Dataset 1 
# Effect of iodine on Arabidopsis thaliana (Col-0) trancriptome
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157643
# GEO Accession: GSE157643



library(GEOquery)
library(limma)
library(umap)
library(dplyr)
library(ggpubr)
library(affy)
library(tidyverse)
library(hrbrthemes)
library(ggplot2)
library(tidyr)
library(reshape)
library(arrayQualityMetrics)
library(matrixStats)


### Step 1 - Feature extraction ----------------------


# Get supplementary files
getGEOSuppFiles("GSE157643")
# Untar files
untar("GSE157643/GSE157643_RAW.tar", exdir = "data/")
# Get raw data
raw.data <- ReadAffy(celfile.path="data/")


# There are 22810 genes:
length(featureNames(raw.data))

#And 8 samples:

length(sampleNames(raw.data))


#Phenotype data:
pData(raw.data)

# It is also possible to verify the feature data for more information.
# There are 1609 repeated Gene Symbol entries:

feature.data <- read.delim2("GPL198-17390.txt", comment.char="#")
length(feature.data$Gene.Symbol)-length(unique(feature.data$Gene.Symbol))

#AGI column has 1277

length(feature.data$AGI[feature.data$AGI == ""])

feature.data <- feature.data %>% mutate_all(na_if,"")

# There are 1277 AGI codes missing, which represents about 5% of the total probe number:
sum(is.na(feature.data$AGI))


#### Step 2 - Quality control ------------------------------------------------

df1 <- exprs(raw.data)

#Change colnames:
colnames(df1) <- c("leaves_control", "leaves_KBr", "leaves_NaI", "leaves_KI",
                   "roots_control", "roots_KBr", "roots_NaI", "roots_KI")
colnames(df1)


### Normality verification --------------------------------

Box-and-whisker plot

par(mar=c(7,4,2,1))
title <- paste ("GSE157643")
boxplot(df1, boxwex=0.7, notch=T, main=title, outline=TRUE, las=2,
        col = c("chartreuse4","chartreuse4","chartreuse4","chartreuse4","burlywood4","burlywood4","burlywood4","burlywood4"))


par(mar=c(7,4,2,1))
title <- paste ("GSE157643 - less outliers")
boxplot(df1, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2,ylim = c(0,1500),
        col = c("chartreuse4","chartreuse4","chartreuse4","chartreuse4","burlywood4","burlywood4","burlywood4","burlywood4"))

# Data is not normally distributed according to shapiro.test since the p-value < 0.01 for all samples of the 8 variables in the dataset:

Df1 <- as.data.frame(df1)
pv.orig = sapply(Df1[,1:8], function(x) shapiro.test(sample(x, 5000, replace = FALSE))$p.value)

colnames(df1)[which(pv.orig < 0.01)]



### Step 3 - Normalization---------------------------------------------------------------
# Since data is not normally distributed, it will be normalized using rma (log 2 base scale);
# rma converts an affybatch object into an expression set using the RMA measure.


normalized.data <- rma(raw.data)

# expression estimates
normalized.expr <- as.data.frame(exprs(normalized.data))
colnames(normalized.expr) <- c("leaves_control", "leaves_KBr", "leaves_NaI", "leaves_KI")

par(mar=c(7,4,2,1))
title <- paste ("GSE157643 after normalization")
boxplot(normalized.expr, boxwex=0.7, notch=T, main=title, outline=TRUE, las=2,
        col = c("chartreuse4","chartreuse4","chartreuse4","chartreuse4","burlywood4","burlywood4","burlywood4","burlywood4"))

#Report production using arrayQualityMetrics:

arrayQualityMetrics(expressionset = normalized.data,
                    outdir = "Report_for_dp",
                    force = TRUE)

#Density plot

dp2 <- melt(normalized.expr)
p2 <- ggplot(data=dp2, aes(x=value, group=variable, fill=variable))+
  ggtitle('Density plot')+
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()
p2



### Feature Data manipulation ----------------------------------

# Feature data contains an ID column (column 1) that corresponds to the row names of the data, 
# and a gene symbol column (column 11). Also, column 12 is the NCBI or Entrez Gene ID.
 
sub_features <- feature.data[, c(1,11,12,14)]


sub_features[1:5,]

#There are 1609 entries in GeneSymbol which are duplicates. This number includes NAs:

dim(sub_features[duplicated(sub_features$Gene.Symbol),])

#Remove missing values:

sf <- sub_features %>% drop_na()


#Merge features with gene expression

data_exprs <- normalized.expr %>%
  rownames_to_column(var = "ID") %>%
  inner_join(., sf, by = "ID")


dim(data_exprs)

# Trying the GWENA low variance filter:

library(GWENA)
Gf <- filter_low_var(t(data_exprs[,1:8]))


# Density plot:

df_f <- melt(t(Gf))
colnames(df_f)[2] <- "Treatment"
p2f <- ggplot(data=df_f, aes(x=value, group=Treatment, fill=Treatment))+
  ggtitle('Density plot')+
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()
p2f


# Select one representative row per group:

library(WGCNA)
collapsed_dataV <- collapseRows(datET = data_exprs[,1:8],
                                rowGroup = data_exprs[,9],
                                rowID = rownames(data_exprs),
                                method = "maxRowVariance")


dat1Collapsed=data.frame(collapsed_dataV$group2row, collapsed_dataV$datETcollapsed)




collapsed_data<- collapseRows(datET = data_exprs[,1:8],
                              rowGroup = data_exprs[,9],
                              rowID = rownames(data_exprs))

dat2Collapsed=data.frame(collapsed_data$group2row, collapsed_data$datETcollapsed)



# Check duplicates 

dat2Collapsed$group[duplicated(dat2Collapsed$group)]

# Save data:
  
#  - GE: gene_expression (no Nans)
#  - WCM: WGCNA collapseRows with MaxMean 
#  - WCmV: WGCNA collapseRows with maxRowVariance 
#  - GF: GWENA filter_low_var

write.csv(data_exprs, "GE.csv")
write.csv(dat1Collapsed, "WCM.csv")
write.csv(dat2Collapsed, "WCmV.csv")
write.csv(Gf, "GF.csv")


