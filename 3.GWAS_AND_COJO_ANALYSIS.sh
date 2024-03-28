#!/bin/bash
#GWAS_CODE

#GWAS analysis of o-methyl ascorbate and ascorbic acid using GCTA fastGWA 

#Rebecca Lelievre 

indir="/home/shared/CLSA/vitCmeta/GWAS/GWAS_FILES/"

##Using fastGWA-mlm function in GCTA
# Only include participants with european ancestry and no missing sex chromosome data 
# For each phenotype file, outlier metabolite values were excluded, sample size differs for each phenotype 
# Adjusted for: age, sex (in pooled analysis), hours since last meal, first 10 principal components,
# province and genotype batch 

 #Pooled Sex Analysis
for i in {ASC,OMETHYL}; do
gcta-1.94.1 --mbfile $indir"GWAS_List.txt" --grm-sparse $indir"sparse_grm_pruned_no_highLD" \
 --fastGWA-mlm --pheno $indir"Pheno_GWAS_${i}.txt" --qcovar $indir"Qcov_Variables_GWAS_Euro_No_Missing_PC.txt" \
  --covar $indir"COV_SEX_PROV.txt" --threads 10 --out ./GWAS_Results/Pooled_GWAS/Province_Adjusted/${i}_GWAS_Province_PC
done

#Interaction Analysis- include sex as an interaction variable 
for i in {ASC,OMETHYL}; do
gcta-1.94.1 --mbfile $indir"GWAS_List.txt" --grm-sparse $indir"sparse_grm_pruned_no_highLD" \
 --fastGWA-mlm --envir $indir"Sex_Cov.txt" --pheno $indir"Pheno_GWAS_${i}.txt" --maf 0.01 --qcovar $indir"Qcov_Variables_GWAS_Euro_No_Missing_PC.txt" \
  --covar $indir"COV_PROV_GWAS.txt" --threads 10 --out ./GWAS_Results/Interaction_Test/Province_Adjusted/${i}_GWAS_Prov_Ineraction
done
 
#Stratified Analysis 
 for i in {ASCF,ASCM,OMETHYLF,OMETHYLM}; do
gcta-1.94.1 --mbfile $indir"GWAS_List.txt" --grm-sparse $indir"sparse_grm_pruned_no_highLD" \
 --fastGWA-mlm --pheno $indir"Pheno_GWAS_${i}.txt" --qcovar $indir"Qcov_Variables_GWAS_Euro_No_Missing_PC.txt" \
  --covar $indir"COV_PROV_GWAS.txt" --threads 10 --out ./GWAS_Results/Stratified_GWAS/Province_Adjusted/${i}_GWAS_Prov
done

 #COJO Analysis

 for i in {10,16}; do 
gcta-1.94.1 --bfile $dir"clsa_imp_${i}_v3_COJO" --chr ${i} --maf 0.01 --cojo-p 5e-8 --cojo-wind 10000 --cojo-collinear 0.9 --cojo-slct \
 --thread-num 4 --cojo-file <(zcat ./GWAS_Results/Pooled_GWAS/Province_Adjusted/ASC_GWAS_Province_PC.fastGWA.gz | awk 'BEGIN{OFS="\t"; print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"} NR>1 {print $2,$4,$5,$7,$8,$9,$10,$6}') \
 --out ./GWAS_Results/Pooled_GWAS/Province_Adjusted/COJO/ASC_COJO_${i}_Prov
 done 

gcta-1.94.1 --bfile $dir"clsa_imp_22_v3_COJO" --chr 22 --maf 0.01 --cojo-p 5e-8 --cojo-wind 10000 --cojo-collinear 0.9 --cojo-slct \
 --thread-num 4 --cojo-file <(zcat ./GWAS_Results/Pooled_GWAS/Province_Adjusted/OMETHYL_GWAS_Province_PC.fastGWA.gz | awk 'BEGIN{OFS="\t"; print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"} NR>1 {print $2,$4,$5,$7,$8,$9,$10,$6}') \
 --out ./GWAS_Results/Pooled_GWAS/Province_Adjusted/COJO/OMETHYL_COJO_22_Prov