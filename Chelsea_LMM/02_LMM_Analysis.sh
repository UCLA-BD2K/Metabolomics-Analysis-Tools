#!/bin/bash

cd /u/home/c/chelseaj/project/Metabolomic/Metabolomics-Analysis-Tools/Chelsea_LMM
echo $PWD

KINSHIP="/u/home/c/chelseaj/project/software/pylmm/scripts/pylmmKinship.py"
GWAS="/u/home/c/chelseaj/project/software/pylmm/scripts/pylmmGWAS.py"

INPUT="Mouse_LMM/Input"
OUTPUT="Mouse_LMM/Output"

GENOTYPE_KINSHIP="$OUTPUT/genotype_kinship.kin"
METABOLITE_KINSHIP="$OUTPUT/metabolite_kinship.kin"
T0_METABOLITE_KINSHIP="$OUTPUT/t0_metabolite_kinship.kin"

SNP_FILE="$INPUT/snp.emma"
PHENOTYPE_FILE="$INPUT/phenotype.txt"
METABOLITE_FILE="$INPUT/metabolite_expression.txt"
T0_METABOLITE_FILE="$INPUT/t0_metabolite.emma"

GENOTYPE_CORR="$OUTPUT/genotype_correction.txt"
METABOLITE_CORR="$OUTPUT/metabolite_correction.txt"
T0_METABOLITE_CORR="$OUTPUT/t0_metabolite_correction.txt"

echo "Building Kinship Matrices"
python $KINSHIP --emmaSNP $SNP_FILE --emmaNumSNPs=243128 $GENOTYPE_KINSHIP
python $KINSHIP --emmaSNP $METABOLITE_FILE --emmaNumSNPs=610 $METABOLITE_KINSHIP
python $KINSHIP --emmaSNP $T0_METABOLITE_FILE --emmaNumSNPs=610 $T0_METABOLITE_KINSHIP

echo "Correcting for Genotype"
#echo "python $GWAS --emmaSNP $METABOLITE_FILE --emmaPHENO $PHENOTYPE_FILE --kfile $GENOTYPE_KINSHIP $GENOTYPE_CORR"
python $GWAS --emmaSNP $METABOLITE_FILE --emmaPHENO $PHENOTYPE_FILE --kfile $GENOTYPE_KINSHIP $GENOTYPE_CORR

echo "Correction for Metabolite"
python $GWAS --emmaSNP $METABOLITE_FILE --emmaPHENO $PHENOTYPE_FILE --kfile $METABOLITE_KINSHIP $METABOLITE_CORR

echo "Correction for T0 Metabolite"
python $GWAS --emmaSNP $METABOLITE_FILE --emmaPHENO $PHENOTYPE_FILE --kfile $T0_METABOLITE_KINSHIP $T0_METABOLITE_CORR


echo "Results writting to $OUTPUT"
echo "DONE"
