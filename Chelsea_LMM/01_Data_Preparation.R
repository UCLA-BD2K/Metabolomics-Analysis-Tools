######################### 01_Data_Preparation.R ############################
# Function: format data for running PyLmm                                  #
# Usage: R --no-save < 01_Data_Preparation.R --args input snp outdir       #
# Arguments: input =  RData                                                #
#	           snp = snp information                                   #
#	           outdir = output directory                               #
# Author: Chelsea Ju                                                       #
# Date: 2015-03-02                                                         #
############################################################################

# read in arguments
options <- commandArgs(trailingOnly = TRUE);
if(length(options) != 3){
  stop(paste("Invalid Arguments\n",
             "Usage: R--no-save --slave < 01_Data_Preparation.R --args input snp outdir\n",
             "\t input = R Data \n",
             "\t snp = snp information \n",
             "\t outdir = output directory\n"),
       sep="");
}

#setwd("/Users/Chelsea/Bioinformatics/PingLab/Metabolomics-Analysis-Tools/Chelsea_LMM/")
data_path <- options[1];
snp_path <- options[2];
outdir <- options[3];

dir.create(file.path(outdir), showWarnings = FALSE, recursive = TRUE);

load(data_path);
snp <- read.table(snp_path, header = T);
colnames(snp) <- c("C57", "A-J", "BALB-cJ", "CE-J", "DBA-2J", "FVB"); # hack strain names for now
#colnames(snp) <- gsub("\\.", '-', colnames(snp));

# reformat SNP Data
# snp for emma: individuals on columns, snps on rows
snp_data <- sapply(Data_Obj_RMRedundancy$label, function(x) snp[,x]);
rownames(snp_data) <- rownames(snp);
colnames(snp_data) <- Data_Obj_RMRedundancy$label;

# SNP output
snp_outfile <- paste(outdir, "snp.emma", sep="/");
write.table(snp_data, file = snp_outfile, sep="\t", col.names = F, row.names = F);

# transpose SNP output
#snp_outfile2 <- paste(outdir, "snp_t.emma", sep="/");
#write.table(t(snp_data), file= snp_outfile2, sep="\t", col.names = F, row.names = F);

# processing metabolites data
Data_Obj_RMRedundancy$quantification[is.na(Data_Obj_RMRedundancy$quantification)] <- 0;

# T0 metabolites
# individuals on columns, metabolites on rows
t0_metabolite <- sapply(Data_Obj_RMRedundancy$label, function(x) Data_Obj_RMRedundancy$quantification[which(Data_Obj_RMRedundancy$label==x & Data_Obj_RMRedundancy$timepoint=="0"),]);
t0_outfile <- paste(outdir, "t0_metabolite.emma", sep="/");
write.table(t0_metabolite, file = t0_outfile, sep="\t", col.names = F, row.names = F);

# observation (metabolite expression) output
# individuals on columns, metabolites on rows  (emma format)
metabolite_outfile <- paste(outdir, "metabolite_expression.txt", sep="/");
write.table(t(Data_Obj_RMRedundancy$quantification), file = metabolite_outfile, sep="\t", col.names = F, row.names = F);

# phenotype output
# phenotypes in one row
phenotype_outfile <- paste(outdir, "phenotype.txt", sep="/");
write.table(t(Data_Obj_RMRedundancy$HWBW), file = phenotype_outfile, sep="\t", col.names = F, row.names = F);

# output metabolites
metabolite_list_outfile <- paste(outdir, "metabolite_name.txt", sep="/");
write.table(Data_Obj_RMRedundancy$metabolite, file = metabolite_list_outfile, sep="\t", col.names = F, row.names = F);

# output individual information
sample_outfile <- paste(outdir, "sample_info.txt", sep="/");
write.table(Data_Obj_RMRedundancy$label, file = sample_outfile, sep="\t", col.names = F, row.names = F);

paste("Files written to ", outdir, "/", sep="");

