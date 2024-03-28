#-------------------Figure 1  ----------------#

##Pooled Sex GWAS Results Manhattan Plots and Suggestive/Signficiant SNP Results Tables 
library(qqman)
#O Methyl
OMET <- read.table("./Pooled_GWAS/Province_Adjusted/OMETHYL_GWAS_Province_PC.fastGWA.gz", header= T)
OMET_SIG <- OMET[OMET$P < 1*10**-5,]
OMET_SIG2 <- OMET[OMET$P < 5*10**-8,]
write.table(OMET_SIG, file="./Pooled_GWAS/Province_Adjusted/OMET_SUG_Prov.txt", row.names=F)
write.table(OMET_SIG2, file="./Pooled_GWAS/Province_Adjusted/OMET_SIG_Prov.txt", row.names=F)
manhattan(OMET, main= "Overall O-Methylascorbate GWAS Manhattan Plot", chr="CHR",  bp="POS", p="P", snp="SNP")

ASC <- read.table("./Pooled_GWAS/Province_Adjusted/ASC_GWAS_Province_PC.fastGWA.gz", header= T)
ASC_SIG <- ASC[ASC$P < 1*10**-5,]
ASC_SIG2 <- ASC[ASC$P < 5*10**-8,]
write.table(ASC_SIG, file="./Pooled_GWAS/Province_Adjusted/ASC_SUG_Prov.txt", row.names=F)
write.table(ASC_SIG2, file="./Pooled_GWAS/Province_Adjusted/ASC_SIG_Prov.txt", row.names=F)
manhattan(ASC, main= "Overall Ascorbic Acid 2 Sulfate GWAS Manhattan Plot", chr="CHR",  bp="POS", p="P", snp="SNP")

#-------------------Figure 2  ----------------#
##Interaction GWAS Results Manhattan Plots and Suggestive/Signficiant SNP Results Tables
#O Methyl
OMET2 <- read.table("./Interaction_Test/Province_Adjusted/OMETHYL_GWAS_Prov_Ineraction.fastGWA.gz", header= T)
OMET_SIG22 <- OMET2[OMET2$P_G_by_E < 1*10**-5,]
OMET_SIG222 <- OMET2[OMET2$P_G_by_E < 5*10**-8,]
write.table(OMET_SIG22, file="./Interaction_Test/Province_Adjusted/OMET_SUG_Prov_Int.txt", row.names=F)
write.table(OMET_SIG222, file="./Interaction_Test/Province_Adjusted/OMET_SIG_Prov_Int.txt", row.names=F)
manhattan(OMET2, main= "O-Methylascorbate Gene-Sex Interaction Manhattan Plot", chr="CHR",  bp="POS", p="P_G_by_E", snp="SNP")

#Ascorbic
ASC2 <- read.table("./Interaction_Test/Province_Adjusted/ASC_GWAS_Prov_Ineraction.fastGWA.gz", header= T)
ASC_SIG22 <- ASC2[ASC2$P_G_by_E < 1*10**-5,]
ASC_SIG222 <- ASC2[ASC2$P_G_by_E < 5*10**-8,]
write.table(ASC_SIG22, file="./Interaction_Test/Province_Adjusted/ASC_SUG_Prov_Int.txt", row.names=F)
write.table(ASC_SIG222, file="./Interaction_Test/Province_Adjusted/ASC_SIG_Prov_Int.txt", row.names=F)
manhattan(ASC2, main= "Ascorbic Acid 2 Sulfate Gene-Sex Interaction Manhattan Plot", chr="CHR",  bp="POS", p="P_G_by_E", snp="SNP")