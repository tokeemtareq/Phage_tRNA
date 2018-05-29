#bin/R

setwd("/media/molgen/scan/10.experiment/05.tRNA_Check/05.Codon_Bias/04.Compiled/")

setwd("/home/tmtareq/05.tRNA_check/05.Codon_Bias/04.Compiled/")

setwd("/media/molgen/New_Volume_E/20.trna_idea/04.Compiled/")

phage_codon = read.table("02.Codon_Composition_Phage.dat", header = F)
colnames(phage_codon) = c("taxid","phage","codon","AA","Fraction","Freq","Number")

host_codon = read.table("03.Codon_Composition_Bacteria_NoAnn.dat",header = F)
colnames(host_codon) = c("taxid","host","genome_no","codon","AA",
                         "Fraction","Freq","Number")


library(ggplot2)
library(reshape2)
theme_set(theme_gray(base_size = 18, base_family = "ubuntu"))

ggplot(data = host_codon, 
       aes(x = codon, y = Fraction, fill = AA)) + 
        geom_boxplot() + 
        xlab("Host_Amino_Acid")

ggplot(data = phage_codon, 
       aes(x = codon, y = Fraction, fill = AA)) + 
        geom_boxplot() + 
        xlab("phage_Amino_Acid")


trna = read.table(file = "01.raw_Phage_tRNA.dat", sep = "\t")
colnames(trna) = c("taxid","phage","amino_acid","anti_codon")

ac = trna$anti_codon
ac = toupper(ac)
ac
c = chartr("ATGC", "TACG", ac)

trna$codon = toupper(c)

codon_to_aa = as.data.frame(table(phage_codon$codon,phage_codon$AA))
#df[df == 0] <- NA
codon_to_aa[codon_to_aa == 0] = NA
codon_to_aa = na.omit(codon_to_aa)
codon_to_aa$Freq = NULL

colnames(codon_to_aa) = c("codon","AA")
trna_phage = merge(trna, codon_to_aa, by = "codon")
trna_phage$amino_acid = NULL
trna_phage$anti_codon = NULL
trna_phage$tRNA = "Presence"


phage_codon_trna = merge(phage_codon, trna_phage, 
                         by = c("AA","codon","taxid","phage"), all = T)
#d[is.na(d)] <- 0
phage_codon_trna[is.na(phage_codon_trna)] = "Absence"



ggplot(data = subset(phage_codon_trna, phage_codon_trna$tRNA == "Presence"), 
       aes(x = phage, y = AA, col = codon, size = Fraction)) + 
        geom_point(alpha = 0.5)

phage_codon_trna2 = subset(phage_codon_trna, 
                          select = c("taxid","phage","AA","codon","Fraction","tRNA"))

colnames(phage_codon_trna2) = c("taxid","phage","AA","codon","phage_fraction","tRNA")


library(dplyr)
#group_by(symbol, date_by_month) %>% summarize

host_codon_mean = host_codon %>% group_by(taxid, host, codon) %>% 
        summarise(mean(Fraction))
#remove(new_dat)
colnames(host_codon_mean) = c("taxid","host","codon","host_fraction")


phage_host_codon_trna = merge(host_codon_mean, phage_codon_trna2, by = c("taxid","codon"))

phct = phage_host_codon_trna

phct$codon_use_diff = abs(phct$host_fraction - phct$phage_fraction)

#women<-prof[which(prof$Sex=="F"),] 
ph_trna = phct[which(phct$tRNA == "Presence"),]
ph_ntrna = phct[which(phct$tRNA == "Absence"),]

summary(ph_trna$codon_use_diff)
summary(ph_ntrna$codon_use_diff)

t.test(ph_ntrna$codon_use_diff, ph_trna$codon_use_diff)

t.test(x = ph_ntrna$codon_use_diff, y = ph_trna$codon_use_diff, 
       alternative = "less")

t.test(x = ph_ntrna$codon_use_diff, y = ph_trna$codon_use_diff, 
       alternative = "greater")

phct_gg = subset(phct, select = c("codon","AA","codon_use_diff","tRNA"))


ggplot(phct_gg, aes(codon, codon_use_diff, col = AA)) + 
        geom_boxplot()

ggplot(phct_gg, aes(codon, codon_use_diff, fill = tRNA, col = tRNA)) + 
        geom_boxplot()

ggplot(phct_gg, aes(codon, codon_use_diff, fill = tRNA, col = AA)) + 
        geom_boxplot() + theme_classic() 

phct = phage_host_codon_trna

phct$codon_use_diff = phct$phage_fraction - phct$host_fraction

#women<-prof[which(prof$Sex=="F"),] 
ph_trna = phct[which(phct$tRNA == "Presence"),]
ph_ntrna = phct[which(phct$tRNA == "Absence"),]

summary(ph_trna$codon_use_diff)
summary(ph_ntrna$codon_use_diff)

t.test(ph_ntrna$codon_use_diff, ph_trna$codon_use_diff)

t.test(x = ph_ntrna$codon_use_diff, y = ph_trna$codon_use_diff, 
       alternative = "less")

t.test(x = ph_ntrna$codon_use_diff, y = ph_trna$codon_use_diff, 
       alternative = "greater")

phct_gg = subset(phct, select = c("codon","AA","codon_use_diff","tRNA"))


ggplot(phct_gg, aes(codon, codon_use_diff, col = AA)) + 
        geom_boxplot(outlier.shape = NA)

ggplot(phct_gg, aes(codon, codon_use_diff, fill = tRNA, col = tRNA)) + 
        geom_boxplot(outlier.shape = NA)

ggplot(phct_gg, aes(codon, codon_use_diff, fill = tRNA, col = AA)) + 
        geom_boxplot(outlier.shape = NA) + theme_classic() 


trna_cod = as.data.frame(table(ph_trna$codon))
ntrna_cod = as.data.frame(table(ph_ntrna$codon))
trna_phage = as.data.frame(table(ph_trna$phage))

trna_phage[trna_phage == 0] = NA
trna_phage = na.omit(trna_phage)

#comparing which trnas are frequently present

colnames(trna_cod) = c("codon","trna_present")
colnames(ntrna_cod) = c("codon","trna_absent")

amino_acid = unique(subset(trna, select = c("amino_acid","codon")))
colnames(amino_acid) = c("a_acid","codon")
amino_acid[amino_acid == ""] = NA
amino_acid = na.omit(amino_acid)

cod_trna = merge(trna_cod, ntrna_cod, by = "codon")
cod_to_trna = merge(cod_trna, amino_acid, by = "codon")

summary(cod_to_trna$a_acid)

#targeting the amino acid use - not codon





