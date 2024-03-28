## Dataset creation and Phenotype/Covariate File Creation ##


#Import Metabolite Dataset (batch-normalized)
library(readr)
CLSA_COMBINED_NORMDATAALL_v1 <- read_csv("~/Research Thesis Etc/Thesis Work/CLSA COMBINED_NORMIMPDATAALL_v1.CSV")
View(CLSA_COMBINED_NORMDATAALL_v1)

library(dplyr)
metabvars <- c("100003258", "100009329", "ADM_METABOLON_COM")
metab_final<-CLSA_COMBINED_NORMDATAALL_v1[metabvars]

#CREATION OF FINAL PHENOTYPE DATSET
##Subset the CLSA dataset with the variables of interest: age, glaucoma, IOP, education, alcohol,smoking, caffeine, fruit and vegetable consumption.
# also include linking keys for metabolic and genetic dataset.

CLSA_Full<- read_csv("~/Research Thesis Etc/Thesis Work/180911_UOttawa_MHRoyGagnon_Baseline_CoPv7_DRU_March2022.csv")
clsa_final<-subset(CLSA_Full, select=c("AGE_NMBR_COM", "ICQ_GLAUC_COM", "ED_HIGH_COM", "ED_OTED_COM", "SDC_CULT_WH_COM", "SDC_CULT_ZH_COM", 
                                       "SDC_CULT_SA_COM", "SDC_CULT_BL_COM", "SDC_CULT_FP_COM",
                                       "SDC_CULT_LA_COM", "SDC_CULT_SE_COM", "SDC_CULT_AR_COM",
                                       "SDC_CULT_WA_COM", "SDC_CULT_JA_COM", "SDC_CULT_KO_COM","SDC_CULT_OT_COM", "SDC_CULT_DK_NA_COM",
                                       "SDC_CULT_REFUSED_COM", "SMK_100CG_COM", "SMK_CURRCG_COM", "ALC_EVER_COM", 
                                       "ALC_FREQ_COM", "ALC_RDWE_NB_COM", "ALC_RDWD_NB_COM", "ALC_WHWD_NB_COM",
                                       "ALC_WHWE_NB_COM", "ALC_BRWD_NB_COM", "ALC_BRWE_NB_COM", "ALC_LQWD_NB_COM", "ALC_LQWE_NB_COM",
                                       "ALC_OTWD_NB_COM", "ALC_OTWE_NB_COM", "TON_IOPCC_R_COM", "TON_IOPCC_L_COM", 
                                       "MEDI_ID_DIN_SP2_1_COM","MEDI_ID_DIN_SP2_2_COM","MEDI_ID_DIN_SP2_3_COM",
                                       "MEDI_ID_DIN_SP2_4_COM","MEDI_ID_DIN_SP2_5_COM","MEDI_ID_DIN_SP2_6_COM","MEDI_ID_DIN_SP2_7_COM",
                                       "MEDI_ID_DIN_SP2_8_COM","MEDI_ID_DIN_SP2_9_COM","MEDI_ID_DIN_SP2_10_COM","MEDI_ID_DIN_SP2_11_COM",
                                       "MEDI_ID_DIN_SP2_12_COM","MEDI_ID_DIN_SP2_13_COM", "MEDI_ID_DIN_SP2_14_COM","MEDI_ID_DIN_SP2_15_COM",
                                       "MEDI_ID_DIN_SP2_16_COM","MEDI_ID_DIN_SP2_17_COM","MEDI_ID_DIN_SP2_18_COM","MEDI_ID_DIN_SP2_19_COM",
                                       "MEDI_ID_DIN_SP2_20_COM","MEDI_ID_DIN_SP2_21_COM","MEDI_ID_DIN_SP2_22_COM","MEDI_ID_DIN_SP2_23_COM",
                                       "MEDI_ID_DIN_SP2_24_COM","MEDI_ID_DIN_SP2_25_COM","MEDI_ID_DIN_SP2_26_COM","MEDI_ID_DIN_SP2_27_COM",
                                       "MEDI_ID_DIN_SP2_28_COM","MEDI_ID_DIN_SP2_29_COM","MEDI_ID_DIN_SP2_30_COM","MEDI_ID_DIN_SP2_31_COM",
                                       "MEDI_ID_DIN_SP2_32_COM","MEDI_ID_DIN_SP2_33_COM","MEDI_ID_DIN_SP2_34_COM","MEDI_ID_DIN_SP2_35_COM",
                                       "MEDI_ID_DIN_SP2_36_COM","MEDI_ID_DIN_SP2_37_COM","MEDI_ID_DIN_SP2_38_COM","MEDI_ID_DIN_SP2_39_COM",
                                       "MEDI_ID_DIN_SP2_40_COM", "NUT_FRUT_NB_COM", "NUT_PURE_NB_COM", "NUT_GREEN_NB_COM", "NUT_CRRT_NB_COM", 
                                       "NUT_VGOT_NB_COM", "SEX_ASK_COM", "BLD_CAFF24_COM", "ADM_METABOLON_COM", "ADM_GWAS3_COM", "WGHTS_PROV_COM"))

#Merge final phenotype dataset with metabolite dataset to create final dataset], total sample is 9992

vitc_final<-merge(metab_final, clsa_final, by="ADM_METABOLON_COM")

colnames(vitc_final)[2]= 'O_METHYL'
colnames(vitc_final)[3]= 'Ascorbic'

#### Subsetting for European Ancestry ######
#Read in file with ancestry information (Sample QC file)
sqc_info<-read.delim("C:/Users/16136/Documents/Research Thesis Etc/Thesis Work/clsa_sqc_v3.txt",sep=" ",header=TRUE)

##Create Dataset with PCA Identifier by Merging with current dataset 

sqc_info2<- subset(sqc_info, select= c("ADM_GWAS_COM", "pca.cluster.id"))

#Rename column to ADM_GWAS3_COM to match with genetic linking key name in phenotypic data
sqc_info2$ADM_GWAS3_COM= sqc_info2$ADM_GWAS_COM

vitc_final<-merge(vitc_final, sqc_info2, by="ADM_GWAS3_COM")
vitc_final<-subset(vitc_final, vitc_final$pca.cluster.id == 4)

#Final dataset has 9052 participants 

#Save a subset of phenotype data with only the variables of interest. 
##At this point, individuals have been subset to those with European ancestry 

metab_GWAS_Variables<- subset(vitc_final, select=c("ADM_METABOLON_COM","ADM_GWAS_COM", "SEX_ASK_COM", "AGE_NMBR_COM", "Ascorbic", "O_METHYL", "WGHTS_PROV_COM"))

#The complete phenotype dataset was read in as "Phenodata"
##Save a subset of variable "Hours since last had food or drink" to add it as a covariate 
metab_hours_since<-subset(CLSA_Full, select=c("ADM_GWAS3_COM", "BLD_FD24_HR_COM"))
colnames(metab_hours_since)<-c("ADM_GWAS_COM", "BLD_FD24_HR_COM")

#Merge hours since meal data with other covariates 
metab_GWAS_Variables<-merge(metab_GWAS_Variables, metab_hours_since, by="ADM_GWAS_COM")

#Incorporate genetic information/add genetic covariates to the file 
## Read in sample quality control file
sqc_info<-read.delim("C:/Users/16136/Documents/Research Thesis Etc/Thesis Work/clsa_sqc_v3.txt",sep=" ",header=TRUE)
## Keep relevant variables 
sqc_variables<- subset(sqc_info, select=c("ADM_GWAS_COM", "batch", "chromosomal.sex", "pca.cluster.id", "in.kinship", "PC1", "PC2", "PC3",
                                          "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

## Merge genetic variables with other covariates of interest
metab_GWAS_Variables<-merge(metab_GWAS_Variables, sqc_variables, by="ADM_GWAS_COM")

#Create a dataset with no missing hours since last meal information (in this case, only European ancestry, but related samples)
data2<- metab_GWAS_Variables %>% filter(BLD_FD24_HR_COM>=0)
write.csv(data2,"C:/Users/16136/Documents/Research Thesis Etc/Metab_GWAS_Variables_No_Missing.csv")

#Create covariate text files for GWAS, stratified by sex 

#Quantitative covariates
Qcov_GWAS <- data2 %>% subset(select=c("ADM_GWAS_COM", "AGE_NMBR_COM", "BLD_FD24_HR_COM"))
Qcov_GWAS$FID = Qcov_GWAS$ADM_GWAS_COM
colnames(Qcov_GWAS)<- c("IID", "AGE", "HOURS_SINCE_MEAL", "FID")
write.table(Qcov_GWAS[,c("FID","IID","AGE","HOURS_SINCE_MEAL")], file="Qcov_Variables_GWAS_Euro_No_Missing.txt", row.names=FALSE, quote=FALSE)

#Categorical Covariates
Cov_GWAS <- data2 %>% subset(select=c("ADM_GWAS_COM", "batch", "WGHTS_PROV_COM" ))
Cov_GWAS$FID = Cov_GWAS$ADM_GWAS_COM
colnames(Cov_GWAS)<- c("IID", "BATCH", "WGHTS_PROV_COM", "FID")
write.table(Cov_GWAS[,c("FID","IID","BATCH", "WGHTS_PROV_COM")], file="COV_PROV_GWAS.txt", row.names=FALSE, quote=FALSE)

## Create File with Principle Components ##
Qcov <- read.table("./Qcov_Variables_GWAS_Euro_No_Missing.txt", header= T)
PC <- subset(metab_GWAS_Variables, select=c("ADM_GWAS_COM", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
names(PC)[names(PC)== 'ADM_GWAS_COM'] <-'FID'
Qcov <- merge(Qcov, PC, by="FID")
write.table(Qcov, file="Qcov_Variables_GWAS_Euro_No_Missing_PC.txt", row.names=FALSE, quote=FALSE)

#Creating Sex Covariate File for Gene-environment interaction test#

Sex_Cov <- subset(metab_GWAS_Variables, select=c("ADM_GWAS_COM", "chromosomal.sex"))
Sex_Cov$FID= Sex_Cov$ADM_GWAS_COM
colnames(Sex_Cov)<-c("IID", "SEX", "FID")
Sex_Cov$SEXN<-with(Sex_Cov, ifelse(SEX==1, 0,
                                   ifelse(SEX==2, 1, NA)))
Sex_Cov2<-na.omit(Sex_Cov)
Sex_Cov2<- subset(Sex_Cov2, select=c("FID", "IID", "SEXN"))
write.table(Sex_Cov2, file="Sex_Cov.txt", row.names=FALSE, quote=FALSE)


#Creating Categorical Covariate Files for pooled sex GWAS #
Cov_GWAS <- data2 %>% subset(select=c("ADM_GWAS_COM", "batch", "chromosomal.sex", "WGHTS_PROV_COM" ))
Cov_GWAS$FID = Cov_GWAS$ADM_GWAS_COM
colnames(Cov_GWAS)<- c("IID", "BATCH","SEX", "WGHTS_PROV_COM", "FID")
Cov_GWAS$SEXN<-with(Cov_GWAS, ifelse(SEX==1, 0,
                                   ifelse(SEX==2, 1, NA)))
write.table(Cov_GWAS[,c("FID","IID","BATCH", "SEXN", "WGHTS_PROV_COM")], file="COV_SEX_PROV.txt", row.names=FALSE, quote=FALSE)



### Phenotype File Creation ###
#Create Phenotype File with European samples, no missing food/drink information, but incl. related samples

#Create phenotype file including only those without missing data 
Pheno_GWAS <- subset(data2, select=c("ADM_METABOLON_COM", "ADM_GWAS_COM", "O_METHYL", "Ascorbic"))

#Create list of metabolites for the loop 
metabolite_list <- unique(colnames(Pheno_GWAS[,3:4]))

#run loop for each metabolite 
for(metabo in metabolite_list) {
  dat_sub<- Pheno_GWAS %>% select(ADM_GWAS_COM, metabo)
  dat_sub$log_metabo<-apply(dat_sub, 1, function(x) log(as.numeric(x[2])))
  mean_metabo <-mean(dat_sub$log_metabo, na.rm=TRUE)
  sd_metabo <-sd(dat_sub$log_metabo, na.rm=TRUE)
  dat_sub_filter<-dat_sub %>% filter(abs((log_metabo-mean_metabo)/sd_metabo)<=3)
  mean_metabo_post <-mean(dat_sub_filter$log_metabo, na.rm=TRUE)
  sd_metabo_post <-sd(dat_sub_filter$log_metabo, na.rm = TRUE)
  dat_sub_filter$std_metabo <-apply(dat_sub_filter, 1, function(x) (as.numeric(x[3])-mean_metabo_post)/sd_metabo_post)
  write.csv(dat_sub_filter, paste0("Pheno_",metabo,"-clsa_europ_.csv"),sep="\t")
}

#Create Phenotype Files by Sex
Pheno_O_METHYL_clsa_europ_ <- read_csv("Pheno_O_METHYL-clsa_europ_.csv")

#O-Methyl Female
Pheno_O_METHYL_clsa_europ_ <- merge(Pheno_O_METHYL_clsa_europ_, data2, by="ADM_GWAS_COM")

Pheno_O_Methyl_GWAS_F<- Pheno_O_METHYL_clsa_europ_ %>% #8916 participants
  filter(chromosomal.sex == 2) %>% #4580 participants
  select(ADM_GWAS_COM, std_metabo)
Pheno_O_Methyl_GWAS_F$FID = Pheno_O_Methyl_GWAS_F$ADM_GWAS_COM

colnames(Pheno_O_Methyl_GWAS_F)<-c("IID", "O-METHYL", "FID")

write.table(Pheno_O_Methyl_GWAS_F[,c("FID","IID","O-METHYL")], file="Pheno_GWAS_OMETHYLF.txt", row.names=FALSE, quote=FALSE)

#O-Methyl Male Sex
Pheno_O_Methyl_GWAS_M<- Pheno_O_METHYL_clsa_europ_ %>% #8916 participants
  filter(chromosomal.sex == 1) %>% #4329 participant
  select(ADM_GWAS_COM, std_metabo)
Pheno_O_Methyl_GWAS_M$FID = Pheno_O_Methyl_GWAS_M$ADM_GWAS_COM

colnames(Pheno_O_Methyl_GWAS_M)<-c("IID", "O-METHYL", "FID")

write.table(Pheno_O_Methyl_GWAS_M[,c("FID","IID","O-METHYL")], file="Pheno_GWAS_OMETHYLM.txt", row.names=FALSE, quote=FALSE)

#7 participants not included in the study with no info on chromosomal sex 

Pheno_asc_clsa_europ_ <- read.csv("Pheno_Ascorbic-clsa_europ_.csv")

#Ascorbic Acid 
## Female sex 
Pheno_asc_clsa_europ_ <- merge(Pheno_asc_clsa_europ_, data2, by="ADM_GWAS_COM")

Pheno_asc_GWAS_F<- Pheno_asc_clsa_europ_ %>% #8835 participants
  filter(chromosomal.sex == 2) %>% #4518 participants 
  select(ADM_GWAS_COM, std_metabo)
Pheno_asc_GWAS_F$FID = Pheno_asc_GWAS_F$ADM_GWAS_COM

colnames(Pheno_asc_GWAS_F)<-c("IID", "ASC", "FID")
write.table(Pheno_asc_GWAS_F[,c("FID","IID","ASC")], file="Pheno_GWAS_ASCF.txt", row.names=FALSE, quote=FALSE)

## Male Sex
Pheno_asc_GWAS_M<- Pheno_asc_clsa_europ_ %>% #8835 participants
  filter(chromosomal.sex == 1) %>% #4310 participants 
  select(ADM_GWAS_COM, std_metabo)
Pheno_asc_GWAS_M$FID = Pheno_asc_GWAS_M$ADM_GWAS_COM

colnames(Pheno_asc_GWAS_M)<-c("IID", "ASC", "FID")
write.table(Pheno_asc_GWAS_M[,c("FID","IID","ASC")], file="Pheno_GWAS_ASCM.txt", row.names=FALSE, quote=FALSE)

## Create Combined Sex Phenotype Files ##

Pheno_asc_clsa_europ_ <- merge(Pheno_asc_clsa_europ_, data2, by="ADM_GWAS_COM")
Pheno_asc_GWAS<- Pheno_asc_clsa_europ_ %>% #8835 participants
  select(ADM_GWAS_COM, std_metabo)
Pheno_asc_GWAS$FID = Pheno_asc_GWAS$ADM_GWAS_COM

colnames(Pheno_asc_GWAS)<-c("IID", "ASC", "FID")
write.table(Pheno_asc_GWAS[,c("FID","IID","ASC")], file="Pheno_GWAS_ASC.txt", row.names=FALSE, quote=FALSE)


Pheno_O_METHYL_clsa_europ_ <- merge(Pheno_O_METHYL_clsa_europ_, data2, by="ADM_GWAS_COM")

Pheno_O_Methyl_GWAS_F<- Pheno_O_METHYL_clsa_europ_ %>% #8916 participants
  select(ADM_GWAS_COM, std_metabo)
Pheno_O_Methyl_GWAS_F$FID = Pheno_O_Methyl_GWAS_F$ADM_GWAS_COM

colnames(Pheno_O_Methyl_GWAS_F)<-c("IID", "O-METHYL", "FID")

write.table(Pheno_O_Methyl_GWAS_F[,c("FID","IID","O-METHYL")], file="Pheno_GWAS_OMETHYL.txt", row.names=FALSE, quote=FALSE)


