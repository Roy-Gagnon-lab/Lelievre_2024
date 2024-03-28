###### Genotype Prep for GWAS Project ######

##################################
# Date Started: July 13th 2023
# Current Date: Jan 9th 2023
##################################

indir="/home/shared/CLSA/GeneticDataRelease3/"
#!/bin/bash

###### GRM Calculation for Metabolites GWAS Project ######

#Generate a List of SNPs from autosomes only from the genotyped files to use to make GRM

plink --bfile $indir"clsa_gen_v3"\
 --chr 1-22 --make-bed --out /home/shared/CLSA/vitCmeta/GRM_Calculation/clsa_gen_v3_nosexchr

###### Generate a list of LD pruned SNPs to use in GRM #########
## Generate while excluding high LD regions from LD pruning 

plink --bfile clsa_gen_v3_nosexchr --exclude high-LD-regions-hg38-GRCh38.txt\
 --geno 0.05 --hwe 0.000001 --indep-pairwise 50 5 0.6 \
 --keep all_euro_plink.txt --keep-allele-order --maf 0.001 \
 --out ./LD_Pruned_SNPs/clsa_gen_v3_nosexchr_pruned_no_highLD

########### Calculate GRM ####################

##Calculate GRM with high LD Regions excluded
# To call GCTA, use the shortcut gcta-1.94.1 

gcta-1.94.1 \
 --bfile clsa_gen_v3_nosexchr \
 --autosome \
 --make-grm \
 --exclude ./LD_Pruned_SNPs/clsa_gen_v3_nosexchr_pruned_no_highLD.prune.out \
 --thread-num 4 \
 --out ./GRM_Europ/grm_pruned_no_highLD \
 
########## Make Sparse GRM ###################

## Make Sparse GRM high LD Regions excluded ##
gcta-1.94.1 --grm ./GRM_Calculation/GRM_Europ/grm_pruned_no_highLD --make-bK-sparse 0.05 --out ./GRM_Calculation/GRM_Europ/sparse_grm_pruned_no_highLD

 

### Create Imputed Genotype files with MAF threshold set to 0.01 ###
# Use these imputed files in GWAS analysis 

for i in {1..22}; {
        plink2 --bgen $indir"clsa_imp_${i}_v3.bgen" ref-first \
        --sample $indir"clsa_imp_v3.sample" --maf 0.01 --mach-r2-filter 0.3 2.0 --make-bed \
        --out ./GWAS_FILES/IMP_GENO/clsa_imp_${i}_v3
}
 

  



 


