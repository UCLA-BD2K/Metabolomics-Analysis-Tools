########################### 03_Data_Summary.R ####################################
# Function: Attach the metabolite name and correct for multiple testing          #
# Usage: R --no-save < 03_Mouse_Data_PostProcessing.R infile metabolite outfile  #
# Arguments: infile = result from pylmm                                          #
#            metabolite = list of metabolites                                    #
#            outfile = output file name                                          #
# Author: Chelsea Ju                                                             #
# Date: 2015-03-04                                                               #
##################################################################################

# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 3){
  stop(paste("Invalid Arguments\n",
             "Usage: R--no-save --slave < 03_Data_Summary.R --args infile metabolite outfile\n",
             "\t infile = result from pylmm \n",
             "\t metabolite = list of metabolites \n",
             "\t outfile = output file name\n"),
       sep="");
}

#setwd("/Users/Chelsea/Bioinformatics/PingLab/Metabolomics-Analysis-Tools/Chelsea_LMM/")
infile <- options[1];
metabolite_file <- options[2];
outfile <- options[3];


indata <- read.table(infile, sep="\t", header = TRUE);
metabolite <- read.table(metabolite_file, header = FALSE);

# attach metabolite
rownames(indata) <- metabolite$V1;

# remove NA
na_row <- which(is.na(indata$BETA));
if(length(na_row) > 0){
	indata <- indata[-na_row,];
}

# adjust for multiple correction
indata$ADJ_P_VALUE <- p.adjust(indata$P_VALUE, method = "BH");

# count
P_001 <- nrow(indata[indata$P_VALUE < 0.01,]);
P_005 <- nrow(indata[indata$P_VALUE < 0.05,]);
ADJ_001 <- nrow(indata[indata$ADJ_P_VALUE < 0.01,]);
ADJ_005 <- nrow(indata[indata$ADJ_P_VALUE < 0.05,]);


summary <- cbind(c(P_001, P_005), c(ADJ_001, ADJ_005));
colnames(summary) <- c("P-value", "Adj-P-value");
rownames(summary) <- c("0.01", "0.05");
summary;

# sort by adj p value
indata <- indata[order(indata$ADJ_P_VALUE, indata$P_VALUE),];
write.table(indata, outfile, sep="\t");





