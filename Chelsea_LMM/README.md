###Download & Install pylmm:
pip install git+https://github.com/nickFurlotte/pylmm --user

###Reference of pylmm:
https://github.com/nickFurlotte/pylmm

###Analysis Steps:
**01: Data Preparation**

```
R --no-save --slave < 01_Data_Preparation.R --args ../../Data/Data_Obj_RMRedundancy_MouseModel.RData ../../Data/df_MySQL_MouseSNP.tsv Mouse_LMM/Input/
```
**02: Running pylmm**

```
bash 02_LMM_Analysis.sh
```
**03: Data Post processing**

```
R --no-save --slave < 03_Data_Summary.R --args Mouse_LMM/Output/genotype_correction.txt Mouse_LMM/Input/metabolite_name.txt Mouse_LMM/Output/genotype_correction_summary.txt
R --no-save --slave < 03_Data_Summary.R --args Mouse_LMM/Output/metabolite_correction.txt Mouse_LMM/Input/metabolite_name.txt Mouse_LMM/Output/metabolite_correction_summary.txt
R --no-save --slave < 03_Data_Summary.R --args Mouse_LMM/Output/t0_metabolite_correction.txt Mouse_LMM/Input/metabolite_name.txt Mouse_LMM/Output/t0_metabolite_correction_summary.txt
```
** Alternatively, pylmm can be run independantly

* Build kinship imatrix
'''
python pylmmKinship.py --emmaSNP Mouse_LMM/Input/snp.emma --emmaNumSNPs=243128 Mouse_LMM/Input/genotype_kinship.kin
python pylmmKinship.py --emmaSNP Mouse_LMM/Input/metabolites.txt --emmaNumSNPs=610 Mouse_LMM/Input/metabolite_kinship.kin
python pylmmKinship.py --emmaSNP Mouse_LMM/Input/t0_metabolite.emma --emmaNumSNPs=610 Mouse_LMM/Input/t0_metabolite_kinship.kin

* Running analysis: correct for genotype
'''
python /u/home/c/chelseaj/project/software/pylmm/scripts/pylmmGWAS.py --emmaSNP Mouse_LMM/Input/metabolite.txt --emmaPHENO Mouse_LMM/Input/phenotype.txt --kfile Mouse_LMM/Input/genotype_kinship.kin -v Mouse_LMM/Output/genotype_correction.txt

* Running analysis: correct for metabolites
'''
python /u/home/c/chelseaj/project/software/pylmm/scripts/pylmmGWAS.py --emmaSNP Mouse_LMM/Input/metabolite.txt --emmaPHENO Mouse_LMM/Input/phenotypes.txt --kfile Mouse_LMM/Input/metabolite_kinship.kin -v Mouse_LMM/Output/metabolite_correction.txt

* Running analysis: correct for T0 metabolites
'''
python /u/home/c/chelseaj/project/software/pylmm/scripts/pylmmGWAS.py --emmaSNP Mouse_LMM/Input/metabolite.txt --emmaPHENO Mouse_LMM/Input/phenotypes.txt --kfile Mouse_LMM/Input/t0_metabolite_kinship.kin -v Mouse_LMM/Output/t0_metabolite_correction.txt
