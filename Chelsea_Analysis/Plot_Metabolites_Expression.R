######################## Plot_Metabolites_Expression.R #####################
# Function: for a given list of metabolites, plot the quantification for   #
#   all patients across different time points                              #
# Usage: R --no-save < Plot_Metabolties_Expression.R --args input          # 
#                                                     metabolties outfile  # 
# Arguments: input = R data                                                # 
#            metabolites = list of metabolites                             #
#            outfile = filename of output                                  #
# Author: Chelsea Ju                                                       #
# Date: 2014-02-16                                                         #
# Modify: 2015-02-28                                                       #
############################################################################

library(RColorBrewer)
library(ggplot2)
library(reshape2)

# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 3){
  stop(paste("Invalid Arguments\n",
             "Usage: R--no-save --slave < Plot_Metabolites_Expression.R --args input metabolites outfile\n",
             "\t input = RData \n",
             "\t metabolites = a list of metabolites \n",
             "\t outfile = filename of outfile\n"),
       sep="");
}

#setwd("/Users/Chelsea/Bioinformatics/PingLab/Metabolomics-Analysis-Tools/Chelsea_Analysis/")

data_path <- options[1];
metabolite_file <- options[2];
outfile <- options[3];

# read data
load(data_path);
metabolite_query <- read.table(metabolite_file, sep="\n");
metabolite_query <- metabolite_query$V1;

# retrieve data
quant <- Data_Obj_RMRedundancy$quantification;
quant_query <- melt(quant[, colnames(quant) %in% metabolite_query])

## reformat data
data <- data.frame(
  patients = rep(Data_Obj_RMRedundancy$subjectNumber, length(metabolite_query)),
  days = rep(Data_Obj_RMRedundancy$fromSurgeryDate, length(metabolite_query)),
  metabolites = quant_query$Var2,
  expression = quant_query$value
);


p <- ggplot(data = data, aes(x=days, y= expression, group = patients, colour=patients)) + 
  geom_point() +
  geom_line()

p + facet_wrap( ~ metabolites, ncol = 10, scales="free")
ggsave(outfile, p);
