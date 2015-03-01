####################### Pairwise_Correlation.R #############################
# Function: compute the pairwise correlation for all metabolites           #  
# Usage: R --no-save < pairwise_correlation.R --args input outdir cut      #
# Arguments: input =  R data                                               #
# Author: Chelsea Ju                                                       #
# Date: 2014-02-12                                                         #
# Modify: 2015-02-28                                                       #
############################################################################

library(RColorBrewer)
library(lattice)
# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 3){
  stop(paste("Invalid Arguments\n",
             "Usage: R--no-save --slave < Pairwise_Correlation.R --args input outdir\n",
             "\t input = R Data \n",
             "\t outdir = output directory\n",
             "\t cut = number of subtree\n"),
       sep="");
}

#setwd("/Users/Chelsea/Bioinformatics/PingLab/Metabolomics-Analysis-Tools/Chelsea_Analysis/")
data_path <- options[1];
outdir <- options[2];
branches <- as.numeric(options[3]);

dir.create(file.path(outdir), showWarnings = FALSE, recursive = TRUE)


load(data_path);

quant <- Data_Obj_RMRedundancy$quantification;
## remove NA from data by setting them to zero
quant[is.na(quant)] <- 0;

# remove quantification with zero sum
zeroSum <- which(colSums(quant) == 0);
if(length(zeroSum) > 0){
  quant <- quant[,-zeroSum];
}

## run pairwise correlation
correlation <- cor(quant);
rgb.palette <- colorRampPalette(brewer.pal(11, "RdBu"));

# clustering order
d <- dist(correlation)
fit <- hclust(d)  # default method = "complete"
pairwise_reorder <- correlation[fit$order, fit$order];

# cut the trees
trees <- cutree(fit, branches);
tree_label <- fit$label[fit$order];

## plot heatmap
pdf(paste(outdir, "pariwise_heatmap.pdf", sep=""))
levelplot(pairwise_reorder, at = unique(c(seq(-1.2, 1.2, length = 100))), col.regions=rgb.palette(100), cuts=50, 
          scales=list(x=list(cex=.1, rot=60), y=list(cex=.1), tck=c(0,0)),
          xlab = "Metabolites", ylab="Metabolites", main = "Pairwise Correlations")
dev.off()


# heirachical cluster dendrogram
pdf(paste(outdir, "pairwise_dendrogram.pdf", sep=""))
plot(fit, cex = 0.1, main = "Cluster of Protein Pairwise Correlation", hang = -1, xlab = "Metabolites")
rect.hclust(fit, k=branches, border="red")
dev.off()

# plot subgroup
for (i in 1:branches){
  subtree <- tree_label[tree_label %in% rownames(as.matrix(trees[trees == i]))] 
  subtree_data <- quant[,subtree]
  
  pdf(paste(outdir, "pairwise_heatmap_subtree_", i, ".pdf", sep=""))
  print(
    levelplot(cor(subtree_data), at = unique(c(seq(-1.2, 1.2, length = 100))), col.regions=rgb.palette(100), cuts=50,
              scales=list(x=list(cex=.2, rot=60), y=list(cex=.2), tck=c(0,0)),
              xlab = "Metabolites", ylab = "Metabolites", main = "Pairwise Correlations")
  )
  dev.off()
  
  write(subtree, file = paste(outdir,"pairwise_list_subtree_", i, ".txt", sep=""))  
}


