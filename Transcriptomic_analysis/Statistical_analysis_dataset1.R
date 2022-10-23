## Dataset 1 
# Statistical analysis 


library(umap)
library(dplyr)
library(ggpubr)
library(affy)
library(tidyverse)
library(affycoretools)
library(hrbrthemes)
library(ggplot2)
library(tidyr)
library(reshape)
library(arrayQualityMetrics)
library(matrixStats)
library(broom)
library(stringr)
library(ggiraphExtra)
library(GWENA)


# Analysis of microarray datasets - part 1

## Dataset 1 
# Effect of iodine on Arabidopsis thaliana (Col-0) trancriptome
# GEO Accession: GSE157643

# This dataset was preprocessed, using rma normalization, and by removing rows that did not have a corresponding Gene symbol. 
# Next, probes were selected in order to maintain the highest variance in data, by using the collapseRows funcion.


WCmV <- read.csv("WCmV.csv", row.names=1)


# Data has the expression information plus 2 columns containing the group of the probe, and the probe ID (selectedRowID):
  
head(WCmV)

# Previously, the dataset had 22810 rows and 8 columns, each corresponding to a treatment. This preprocessing steps led to the exclusion of 1762 rows, thus leaving a total of 21057 rows:
  
dim(WCmV)

# Subseting the expression data only:
  
dfe <- WCmV[,3:10]

# The limma package was studied and tested, however, since there is only one file for each treatment on the 2 different tissues (leaves and roots) it is not possible to obtain robust results (no replicates). 

# Since data was previously normalized it is possible to use parametric tests to study differences between groups. 
# By performing the Kolmorov-Smirnov test on all samples, it is possible to admit that they follow a normal distribution (p>0.05). 

lst.ks <- lapply(1:ncol(dfe), function(i)
  ks.test(dfe[, i], "pexp", 1.0/mean(dfe[, i])))
Dstat <- sapply(lst.ks, function(x) x$statistic)
Dstat


# ANOVA

# Separating leaves and roots into two dataframes: df_leaves and df_roots

df_leaves <- dfe[,1:4]
df_roots <- dfe[,5:8]

df_leaves <- df_leaves %>%
  rownames_to_column(var = "ID")
df_roots <- df_roots %>%
  rownames_to_column(var = "ID")


## Testing the interaction between treatment and tissue:

# Y = u + Tissue*Treatment + error

df_all <- data.frame(df_leaves, df_roots[,2:5])
df_all <- melt(df_all)
df_all[c("Tissue","Treatment")]<-str_split_fixed(df_all$variable,"_",2)

colnames(df_all) <- c("ID","variable","Expression","Tissue","Treatment")
df_all[1:3,]

# Stating the hypothesis:
  
#  Test whether there is an interaction on tissue and treatment affecting gene expression levels.

anova1 <- aov(Expression ~ Treatment*Tissue, data = df_all)

summary(anova1)

# In the table above the "Treatment:Tissue" variable has a p-value=0.0322.
# This means there the interaction is significant and a source of variation on gene expression level. 
# In fact, it is expected that there is a variation on expression levels on different plant tissues. 

tukey1<-TukeyHSD(anova1)

tukey.plot.aov<-aov(Expression ~ Tissue:Treatment, data=df_all)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)

rownames(tukey.plot.test$`Tissue:Treatment`) <- gsub("roots","R",rownames(tukey.plot.test$`Tissue:Treatment`))
rownames(tukey.plot.test$`Tissue:Treatment`) <- gsub("control","C",rownames(tukey.plot.test$`Tissue:Treatment`))
rownames(tukey.plot.test$`Tissue:Treatment`) <- gsub("leaves","L",rownames(tukey.plot.test$`Tissue:Treatment`))

ggHSD(tukey.plot.test)+
  theme(plot.title=element_text(size=rel(2)),
        axis.text.y=element_text(angle = 0, size=rel(1.5),hjust=0.5),
        axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)))

# The Tuckey test confirms that, with 95% confidence, there is a significant difference between gene expression values in roots and leaves (p-value = 3.1e-06). 
# From the Tuckey test it is possible to verify that there are statistically significant differences (p<0.05) between the following groups:
  
# 1. KBr:roots-control:leaves (p-value=0.0000011)
# 2. KBr:roots-KBr:leaves (p-value=0.0000025)
# 3. KBr:roots-KI:leaves (p-value=0.0000000)
# 4. NaI:roots-KI:leaves (p-value=0.0011947)
# 5. KBr:roots-NaI:leaves (p-value=0.0000042)
# 6. NaI:roots-control:leaves (p-value=0.0109511)
# 7. NaI:roots-KBr:leaves (p-value=0.0176926)
# 8. NaI:roots-NaI:leaves (p-value=0.0242985)
# 9. KBr:roots-control:roots (p-value=0.0122500)
# 10. KI:roots-KBr:roots (p-value=0.0001377)


# The first 8 groups show significant differences in treatments in different tissues. 
# The second group shows that there is a significant difference within the same treatment (KBr) between roots and leaves. 
# Groups 9 and 10 show significant differences on the same tissue (roots) between treatments with KBr and the control treatment,
# and that KI and KBr have caused significant differences of expression in the roots. 
# This can mean that KBr is causing significant differences in expression levels in the roots. 
# It is not likely that this difference is attributed to the effect of KI since there are no differences detected here between this treatment and the control on the roots. 
# Other tests will be performed by conducting separate analysis on the leaves and roots. 


## Testing only differences in treatment on Leaves:

DFL <-melt(df_leaves)
head(DFL)


colnames(DFL)<-c("ID","Treatment","Expression")
DFL$Treatment <- gsub("leaves_","",DFL$Treatment)
head(DFL)

# Stating the hypothesis:
# H0: The mean value of expression in the different treatments is the same
# H1: The mean value of expression is different between treatments. 

leaves_anova1 <- aov(Expression ~ Treatment, data = DFL)

summary(leaves_anova1)

# There are 3 degrees of freedom (number of treatments-1). The sum of squares displays the total variation between the means of treatments and the overall mean. 
# P-value = 0.85 (p-value>0.05), therefore, there are no significant differences between treatments in the leaves. 

tukey_leaves1<-TukeyHSD(leaves_anova1)
ggHSD(tukey_leaves1)+
  theme(plot.title=element_text(size=rel(2)),
        axis.text.y=element_text(angle = 0, size=rel(1.5),hjust=0.5),
        axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)))


# Roots:
  
DFR <-melt(df_roots)
head(DFR)

colnames(DFR)<-c("ID","Treatment","Expression")
DFR$Treatment <- gsub("roots_","",DFR$Treatment)


# Stating the hypothesis:
# H0: The mean value of expression on roots in the different treatments is the same
# H1: The mean value of expression on roots is different between treatments. 

roots_anova1 <- aov(Expression ~ Treatment, data = DFR)

summary(roots_anova1)

# We can conclude that there is a significant difference on gene expression levels on the roots between the treatments (p-value<0.05). 


# As a post-hoc test, the Tukey's Honestly Significant Difference (Tuckey's HSD) test can be used to perform pairwise-comparisons between treatments, 
# and therefore discover which groups are statistically different. 

tukey_roots1<-TukeyHSD(roots_anova1)

tukey_roots1

# From the Tukey's HSD test it is possible to observe that there are statistically significant differences (95% confidence) between:
# - roots treated with KBr and the control (p-value = 0.0031584)
# - roots treated with KI and KBr (p-value = 0.0000342)
 
 
# These differences can be clearly observed in the following plot:

ggHSD(tukey_roots1)+
      theme(plot.title=element_text(size=rel(2)),
              axis.text.y=element_text(angle = 0, size=rel(1.5),hjust=0.5),
              axis.text.x=element_text(size=rel(1.5)),
              axis.title.x=element_text(size=rel(1.5)))




# Testing a filter before ANOVA (filter then repeat)
# Genes with constant values may not add relevant information. 


# Using GWENA package

library(GWENA)
library(magrittr)

threads_to_use <- 2

# Transpose data:

tdfe <- t(dfe)


# checking if expression data is correctly defined
is_data_expr(tdfe)
design <- read.delim("design.txt")

# Filter low variance data:
# Using this filtering, 14739 genes remain (6318 genes were removed).

Exp_filter <- filter_low_var(tdfe, pct = 0.7, type = "median")

# Remaining number of genes
ncol(tdfe)
ncol(Exp_filter)


MEF <- melt(t(Exp_filter))


MEF[c("Tissue","Treatment")]<-str_split_fixed(MEF$X2,"_",2)
MEF <- data.frame(MEF$X1,MEF$value,MEF$Tissue,MEF$Treatment)
colnames(MEF)<-c("ID","Expression","Tissue","Treatment")


fanova <- aov(Expression ~ Treatment*Tissue, data = MEF)

summary(fanova)

# The interaction is no longer significative. 


tukey_f<-TukeyHSD(fanova)

tukey_f

DFL2 <- melt(DFL2)
colnames(DFL2)<-c("ID","Treatment","Expression")
DFL2$Treatment <- gsub("leaves_","",DFL2$Treatment)

leaves_anova2 <- aov(Expression ~ Treatment, data = DFL2)
summary(leaves_anova2)

tukey_leaves2<-TukeyHSD(leaves_anova2)
tukey_leaves2
ggHSD(tukey_leaves2)+
        theme(plot.title=element_text(size=rel(2)),
        axis.text.y=element_text(angle = 0, size=rel(1.5),hjust=0.5),
        axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)))



### Performing a second ANOVA test with roots data

DFR2 <- t(Exp_filter)
DFR2<-DFR2[,5:8]
DFR2 <- melt(DFR2)
colnames(DFR2)<-c("ID","Treatment","Expression")
DFR2$Treatment <- gsub("roots_","",DFR2$Treatment)
roots_anova2 <- aov(Expression ~ Treatment, data = DFR2)
summary(roots_anova2)
tukey_roots2<-TukeyHSD(roots_anova2)
tukey_roots2
ggHSD(tukey_roots2)+
        theme(plot.title=element_text(size=rel(2)),
        axis.text.y=element_text(angle = 0, size=rel(1.5),hjust=0.5),
        axis.text.x=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=rel(1.5)))
