######################## Plot_Metabolites_FoldChange.R #####################
# Function: for a given list of metabolites, plot the fold change          #
#           quantification for all patients across different time points   #
# Usage: R --no-save < Plot_Metabolties_FoldChange.R --args input          # 
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
quant[is.na(quant)] <- 0;
quant_query <- quant[, colnames(quant) %in% metabolite_query];

T0_quant <- quant_query;
for(i in 1:nrow(T0_quant)){
  patient_id <- as.character(Data_Obj_RMRedundancy$subjectNumber[i]);
#  row_id <- which(Data_Obj_RMRedundancy$subjectNumber==patient_id & Data_Obj_RMRedundancy$timepoint=="T0")
  row_id <- which(Data_Obj_RMRedundancy$subjectNumber==patient_id)[1] #P25 doesn't have T0
  T0_quant[i,] <- quant_query[row_id,];
}

## to avoid 0, add 1 to everything
T0_quant <- T0_quant + 1;
quant_query <- quant_query + 1;
quant_fc <- quant_query / T0_quant;
quant_log_fc <- log(quant_fc, 2);
quant_log_fc <- melt(quant_log_fc);

## reformat data
data <- data.frame(
  patients = rep(Data_Obj_RMRedundancy$subjectNumber, length(metabolite_query)),
  days = rep(Data_Obj_RMRedundancy$fromSurgeryDate, length(metabolite_query)),
  metabolites = quant_log_fc$Var2,
  expression = quant_log_fc$value
);

pdf(file=outfile, width=11, height=11)
p <- ggplot(data = data, aes(x=days, y= expression, group = patients, colour=patients)) + 
  geom_point() +
  geom_line() +
  xlab("Days from Sugery") +
  ylab("log2-based Fold Change vs T0") + 
  facet_wrap( ~ metabolites, ncol = 10, scales="free") + 
  theme(strip.text.x = element_text(size = 6, angle=10)) + 
  theme(axis.text.x = element_text(size = 5, angle=45)) + 
  theme(axis.text.y = element_text(size = 5)) +
  theme(axis.title.x = element_text(size=12)) + 
  theme(axis.title.y = element_text(size=12))  

print(p)
dev.off();
