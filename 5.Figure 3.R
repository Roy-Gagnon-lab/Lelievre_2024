#Stratified GWAS Effect Size of Lead SNPs Figure 3 

#Create Combined Dataset with stratified effect sizes:

### O-Methyl Ascorbate

#Read in Lead SNPs from Interaction 
### Created separate CSV File with list of SNPs that were lead SNPs based on FUMA criteria

OMET_Inter3<- read.csv("C:/Users/16136/Documents/Research Thesis Etc/Thesis Work/GWAS/GWAS_Results/Interaction Test Results/Lead_Interaction_SNPs_OMET.csv", header=T)
#Read in Female O-Methyl GWAS results 
F_Omet_GWAS <- read.table("C:/Users/16136/Documents/Research Thesis Etc/Thesis Work/GWAS/GWAS_Results/Second Attempt/OMETHYLF_GWAS_Prov.fastGWA", header=T)
#Keep only SNP and beta information to merge with interaction dataset
F_Omet_GWAS <-subset(F_Omet_GWAS, select=c("SNP","BETA","SE"))
OMET_Inter_F <- merge(F_Omet_GWAS,OMET_Inter3, by="SNP")
#rename p-value column to identify as female p-value
names(OMET_Inter_F)[names(OMET_Inter_F)== 'BETA'] <-'BETA_F'
names(OMET_Inter_F)[names(OMET_Inter_F)== 'SE'] <-'SE_F'
#Read in male GWAS results 
M_Omet_GWAS <- read.table("C:/Users/16136/Documents/Research Thesis Etc/Thesis Work/GWAS/GWAS_Results/Second Attempt/OMETHYLM_GWAS_Prov.fastGWA", header=T)
#Subset to only SNP and BETA information and merge with female-interaction dataset by SNP name
M_Omet_GWAS <-subset(M_Omet_GWAS, select=c("SNP","BETA","SE"))
OMET_Inter_F_M <- merge(M_Omet_GWAS,OMET_Inter_F, by="SNP")
names(OMET_Inter_F_M)[names(OMET_Inter_F_M)== 'SE'] <-'SE_M'
names(OMET_Inter_F_M)[names(OMET_Inter_F_M)== 'BETA'] <-'BETA_M'
write.csv(OMET_Inter_F_M, file="OMET_Interaction_F_M_PC.csv", row.names = F)

## In excel, reformatted dataset to be able to create plots, rearranged to have one row per sex-SNP and 
## calculate 95% CI for each SNP
## Final dataset looks like: SNP, Sex, SE, Lower 95% CI, Upper 95% CI
## Read in Dataset 
OMET_Interaction_Plot<-read.csv("C:/Users/16136/Documents/Research Thesis Etc/Thesis Work/GxE/OMET_Interaction_F_M_PC.csv", header = T)

#Create Stratified Effects Plot using GGPLOT

library(ggplot2)

OMET_PLOT<-ggplot(OMET_Interaction_Plot, aes(x=SNP, y=BETA_M, ymin=Lower, ymax=Upper, col=Sex)) +
  geom_linerange(linewidth=2, position=position_dodge(width=0.5)) +
  geom_hline(yintercept = 0, lty=2) +
  geom_point(size=2)+theme_bw() +
  scale_x_discrete(name="SNPs") +
  scale_y_continuous(name="Effect Size", limits=c(-0.6,0.6)) +
  coord_flip()

###Ascorbic Acid

#Read in Lead SNPs from Interaction 
### Created separate CSV File with list of SNPs that were lead SNPs based on FUMA criteria
ASC_Inter3<- read.csv("C:/Users/16136/Documents/Research Thesis Etc/Thesis Work/GWAS/GWAS_Results/Interaction Test Results/Lead_Interaction_SNPs_ASC.csv", header=T)
#Read in Female O-Methyl GWAS results 
F_ASC_GWAS <- read.table("C:/Users/16136/Documents/Research Thesis Etc/Thesis Work/GWAS/GWAS_Results/Second Attempt/ASCF_GWAS_Prov.fastGWA", header=T)
#Keep only SNP and beta information to merge with interaction dataset
F_ASC_GWAS <-subset(F_ASC_GWAS, select=c("SNP","BETA","SE"))
ASC_Inter_F <- merge(F_ASC_GWAS,ASC_Inter3, by="SNP")
#rename p-value column to identify as female p-value
names(ASC_Inter_F)[names(ASC_Inter_F)== 'BETA'] <-'BETA_F'
names(ASC_Inter_F)[names(ASC_Inter_F)== 'SE'] <-'SE_F'
#Read in male GWAS results 
M_ASC_GWAS <- read.table("C:/Users/16136/Documents/Research Thesis Etc/Thesis Work/GWAS/GWAS_Results/Second Attempt/ASCM_GWAS_Prov.fastGWA", header=T)
#Subset to only SNP and BETA information and merge with female-interaction dataset by SNP name
M_ASC_GWAS <-subset(M_ASC_GWAS, select=c("SNP","BETA","SE"))
ASc_Inter_F_M <- merge(M_ASC_GWAS,ASC_Inter_F, by="SNP")
names(ASc_Inter_F_M)[names(ASc_Inter_F_M)== 'SE'] <-'SE_M'
names(ASc_Inter_F_M)[names(ASc_Inter_F_M)== 'BETA'] <-'BETA_M'
write.csv(ASc_Inter_F_M, file="ASc_Interaction_F_M_PC2.csv", row.names = F)

## In excel, reformatted dataset to be able to create plots, rearranged to have one row per sex-SNP and 
## calculate 95% CI for each SNP
## Final dataset looks like: SNP, Sex, SE, Lower 95% CI, Upper 95% CI
## Read in Dataset 
ASC_Interaction_Plot<-read.csv("C:/Users/16136/Documents/Research Thesis Etc/Thesis Work/GxE/ASc_Interaction_F_M_PC.csv", header = T)

ASC_PLOT<-ggplot(ASC_Interaction_Plot, aes(x=SNP, y=BETA_M, ymin=Lower, ymax=Upper, col=Sex)) +
  geom_linerange(size=2, position=position_dodge(width=0.5)) +
  geom_hline(yintercept = 0, lty=2) +
  geom_point(size=2)+theme_bw() +
  scale_x_discrete(name="SNPs") +
  scale_y_continuous(name="Effect Size", limits=c(-0.7,0.7)) +
  coord_flip()

