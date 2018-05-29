###############################################################################
###############################################################################

#setwd("/home/tmtareq/05.tRNA_check/05.Codon_Bias/10.Analysis_Base/")

#setwd("G:/STUDY/10.tRNA_idea/10.Analysis_Base/")

setwd("/media/molgen/New_Volume_E/20.trna_idea/10.Analysis_Base/")

library(ggplot2)
library(reshape2)
library(dplyr)

################################################################
################################################################
#             Loding the reformed data

#need to do it in the office laptop

# changing the file name in the folder to separate these from previous attempts
#
# phage-tRNA is denoted as 40.01
# phage-bac-codon is denoted as 40.02
# loading and merging the data

ph_trna = read.table("40.01.phage_tRNA_codon", header = T, sep = "\t")
ph_bac_codon = read.table("40.02.host_phage_codon", header = T, sep = "\t")


ph_bac_trna_codon = merge(ph_trna, ph_bac_codon, by = c("taxid","phage","AA","codon","pFraction","pFreq","pNumber"))
dat0 = ph_bac_trna_codon

# writing to a file and uploading in the google drive
#write.csv(dat0, file = "40.03.Phage_Bac_Codon_tRNA.csv", quote = F, row.names = F)

# let's start the actual work
# Starting with part one

###########################################################################
# THE PHAGE TRNA EXISTS TO CONPENSATE THE PHAGE AND ITS HOST BACTERIA CODON
# OR AMINO ACID USAGE BIAS
########################################################################### 

# to simplify, using only one genome per bacteria

dat1 = subset(dat0, dat0$bGenome == 1, select = c("taxid","phage","AA","codon","pFraction","pFreq","pNumber","ptRNA","bacteria","bGenome","bFraction","bFreq","bNumber","number"))

#using dplyr to cumulate the data for AA
#x %>% group_by(Category) %>%  summarise(Frequency = sum(Frequency))

dat1.01 = dat1 %>% group_by(bacteria, AA) %>% 
        summarize(sum(bNumber), sum(pNumber))
colnames(dat1.01) =c("bacteria","AA","AAbNumber","AApNumber")

dat1.02 = dat1.01 %>% group_by(bacteria) %>%
        summarise(sum(AAbNumber), sum(AApNumber))
colnames(dat1.02) = c("bacteria","AAbSum","AApSum")

dat1.03 = merge(dat1.01, dat1.02, by = "bacteria")
dat1.03$AAbFreq = dat1.03$AAbNumber/dat1.03$AAbSum*1000
dat1.03$AApFreq = dat1.03$AApNumber/dat1.03$AApSum*1000

#adding the AA fraction with the previous dataset
dat1.04 = merge(dat1, dat1.03, by = c("bacteria","AA"))

#calculating the difference of the usage
dat1.04$diff_CFreq = dat1.04$pFreq - dat1.04$bFreq
dat1.04$diff_AAFreq = dat1.04$AApFreq - dat1.04$AAbFreq
dat1.04$diff_CFraction = dat1.04$pFraction - dat1.04$bFraction
# separating the tRNA containing phage and removing stop codon from the data

#removing the stop codons
dat1.10 = subset(dat1.04, dat1.04$codon != "*")

#creating the list of the phages that contain at least one tRNA
dat1.11 = unique(subset(dat1.10, dat1.10$ptRNA == "Presence", select = c("phage")))

#the dataset contain all the phages with at least one tRNA
dat1.12 = merge(dat1.11, dat1.10, by = "phage")

#dat1.15 contain the phages with no tRNA
dat1.13 = dat1.11
dat1.13$missing = "tRNA"
dat1.14 = merge(dat1.13, dat1.10, by = "phage", all = T)
dat1.15 = subset(dat1.14, is.na(missing))
dat1.15$ptRNA = "None"
dat1.15$missing = NULL


#reattaching the data
dat1.16 = rbind(dat1.15,dat1.12)

# preparing the data for t-test
#dat1.12 %>% group_by(ptRNA) %>% summarise(mean(diff_AAFreq))

t.test(dat1.12$diff_CFreq ~ dat1.12$ptRNA, alternative = "less", var.equal = TRUE)
t.test(dat1.12$diff_AAFreq ~ dat1.12$ptRNA, alternative = "less", var.equal = TRUE)
t.test(dat1.12$diff_CFraction ~ dat1.12$ptRNA, alternative = "less", var.equal = TRUE)




#for codon Frequency
df_C_Freq_a = as.vector(subset(dat1.16, dat1.16$ptRNA == "Absence", select = "diff_CFreq"))
df_C_Freq_p = as.vector(subset(dat1.16, dat1.16$ptRNA == "Presence", select = "diff_CFreq"))
df_C_Freq_n = as.vector(subset(dat1.16, dat1.16$ptRNA == "None", select = "diff_CFreq"))
t.test(x = df_C_Freq_p, y = df_C_Freq_a, alternative = "g", var.equal = T)
t.test(x = df_C_Freq_p, y = df_C_Freq_n, alternative = "g", var.equal = T)

#for Codon Fraction
df_C_Frac_a = as.vector(subset(dat1.16, dat1.16$ptRNA == "Absence", select = "diff_CFraction"))
df_C_Frac_p = as.vector(subset(dat1.16, dat1.16$ptRNA == "Presence", select = "diff_CFraction"))
df_C_Frac_n = as.vector(subset(dat1.16, dat1.16$ptRNA == "None", select = "diff_CFraction"))
t.test(x = df_C_Frac_p, y = df_C_Frac_a, alternative = "g", var.equal = T)
t.test(x = df_C_Frac_p, y = df_C_Frac_n, alternative = "g", var.equal = T)


#for AA Freq
df_AA_Freq_a = as.vector(subset(dat1.16, dat1.16$ptRNA == "Absence", select = "diff_AAFreq"))
df_AA_Freq_p = as.vector(subset(dat1.16, dat1.16$ptRNA == "Presence", select = "diff_AAFreq"))
df_AA_Freq_n = as.vector(subset(dat1.16, dat1.16$ptRNA == "None", select = "diff_AAFreq"))
t.test(x = df_AA_Freq_p, y = df_AA_Freq_a, alternative = "g", var.equal = T)
t.test(x = df_AA_Freq_p, y = df_AA_Freq_n, alternative = "g", var.equal = T)


#generating some graphs

ggplot(data = dat1.16, aes(x = codon, y = diff_CFraction, fill = ptRNA)) + 
        geom_boxplot(outlier.shape = NA)

ggplot(data = dat1.16, aes(x = codon, y = diff_CFreq, fill = ptRNA)) + 
        geom_boxplot(outlier.shape = NA)

ggplot(data = dat1.16, aes(x = codon, y = diff_AAFreq, fill = ptRNA)) + 
        geom_boxplot(outlier.shape = NA)


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# it looks like the significance is much lower for codon
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


###########################################################################
# THE NUMBER OF TRNA CORRELATES WITH THE DIFFERENCE
########################################################################### 

# for this, need to do the ANOVA
# to determine whether the number of the tRNA correlates with the 
# diiference in AA or codon frequency

#only keeping the 1st genome and keeping the phage codon, 
#with at least one tRNA in their genome
dat2.00 = subset(dat1.16, dat1.16$ptRNA == "Presence")

#basic correlation to check the possibility

cor(x = dat2.00$number, y = dat2.00$diff_CFraction)
cor(x = dat2.00$number, y = dat2.00$diff_CFreq)
cor(x = dat2.00$number, y = dat2.00$diff_AAFreq)

#maybe there is not any connection, gotta do the ANOVA anyway
















