# Expression Analysis - Dataset1

library(ggplot2)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(stats)
library(genefilter)
library(GWENA)
library(WGCNA)
library(limma)
library(reshape2)
library(EnhancedVolcano)
library(heatmaply)
library(RColorBrewer) 


# Read data ------------------------------

# collapsedRows data
WCmV <- read.csv("WCmV.csv", row.names=1)
df <- WCmV[,3:10]

# Design matrix
design <- read.delim("design.txt", row.names=1)


# Building the Expression Sets: --------------------
# Building an expression set data with filtered genes:

# Collapsed rows:
exp_data <- as.matrix(df)
pheno_data <- new("AnnotatedDataFrame", data=design)
all_set<- ExpressionSet(assayData=exp_data,
                        phenoData=pheno_data)


# subset with only interest genes from phylogeny--------

igenes <- c("DUR3","AT1G73020", "SULTR4;1", "SULTR3;1", 
            "SULTR2;1", "SULTR4;2", "SULTR1;1", "AST56", 
            "AST91", "SULTR1;3", "SULTR3;4", "SULTR1;2",
            "SULTR3;5", "SULTR3;2", "VPS60.1", "VPS60.2", 
            "VPS2.3", "SNF7.1", "VPS20.1", "VPS2.1","SNF7.2", "ATMRP2", "MRP1", 
            "ATMRP5","ATMRP7", "ATMRP4", "ATMRP10", 
            "ATMRP14", "ATNAP5", "ATMRP8", "ATMRP9", 
            "ATMRP11", "MRP6", "ABCB4", "PGP5", "ATMRP13")

exp_data <- as.matrix(df[igenes,])
pheno_data <- new("AnnotatedDataFrame", data=design)
all_set<- ExpressionSet(assayData=exp_data,
                        phenoData=pheno_data)


# LEAVES AND ROOTS SEPARATION 
leaves <- all_set[, all_set$Tissue =="L"]
#writexl:write_xlsx(as.data.frame(exprs(leaves)), "leaves.xlsx")
roots <- all_set[, all_set$Tissue =="R"]
#writexl:write_xlsx(as.data.frame(exprs(roots)), "roots.xlsx")

colMain <- rev(colorRampPalette(brewer.pal(8, "Spectral"))(25))
heatmaply(exprs(leaves), 
          main = "Expression of genes detected in the Phylogenetic analysis - leaves",
          margins = c(60,100,40,20),
          grid_color = "white",
          col = colMain,
          grid_width = 0.00001,
          label_names = c("Gene", "Sample:", "Value"),
          fontsize_row = 10, fontsize_col = 10,
          labCol = c("Control","KBr","NaI","KI"),
          labRow = rownames(exprs(leaves)),
          heatmap_layers = theme(axis.line=element_blank()), 
          dendrogram = "none"
)


heatmaply(exprs(roots), 
          xlab = "", ylab = "", 
          main = "Expression of genes detected in the Phylogenetic analysis - roots",
          margins = c(60,100,40,20),
          grid_color = "white",
          col = colMain,
          grid_width = 0.00001,
          titleX = FALSE,
          branches_lwd = 0.1,
          label_names = c("Gene", "Sample:", "Value"),
          fontsize_row = 10, fontsize_col = 10,
          labCol = c("Control","KBr","NaI","KI"),
          labRow = rownames(exprs(roots)),
          heatmap_layers = theme(axis.line=element_blank()), 
          dendrogram = "none"
)
