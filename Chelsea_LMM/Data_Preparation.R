######################### Data_Preparation.R ###############################
# Function: format data for running PyLmm                                  #
# Usage: R --no-save < Data_Preparation.R --args input snp outdir          #
# Arguments: input =  RData                                                #
#	           snp = snp information                                         #
#	           outdir = output directory                                     #
# Author: Chelsea Ju                                                       #
# Date: 2015-03-02                                                         #
############################################################################

# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 3){
  stop(paste("Invalid Arguments\n",
             "Usage: R--no-save --slave < Data_Preparation.R --args input snp outdir\n",
             "\t input = R Data \n",
             "\t snp = snp information \n",
             "\t outdir = output directory\n"),
       sep="");
}

#setwd("/Users/Chelsea/Bioinformatics/PingLab/Metabolomics-Analysis-Tools/Chelsea_LMM/")
data_path <- options[1];
snp_path <- options[2];
outdir <- options[3]

dir.create(file.path(outdir), showWarnings = FALSE, recursive = TRUE)

load(data_path);
snp <- read.table(snp_path, header = T);
colnames(snp) <- gsub("\\.", '-', colnames(snp));

# reformat SNP Data
snp_data <- matrix(0, nrow = length(Data_Obj_RMRedundancy$label), ncol = nrow(snp));
for(i in 1:length(Data_Obj_RMRedundancy$label)){
  snp_data[i,] <- t(snp[,Data_Obj_RMRedundancy$label[i]]);
}
colnames(snp_data) <- rownames(snp);
rownames(snp_data) <- Data_Obj_RMRedundancy$label;

# SNP output
snp_outfile <- paste(outdir, "snp.emma", sep="/");
write.table(snp_data, file = snp_outfile, sep="\t");

# phenotype output
metabolite_outfile <- paste(outdir, "phenotype.txt", sep="/");
write.table(Data_Obj_RMRedundancy$quantification, file = metabolite_outfile, sep="\t");
