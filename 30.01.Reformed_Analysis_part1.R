###############################################################################
###############################################################################

setwd("/home/tmtareq/05.tRNA_check/05.Codon_Bias/10.Analysis_Base/")

setwd("G:/STUDY/10.tRNA_idea/10.Analysis_Base/")

library(ggplot2)
library(reshape2)
library(dplyr)

# loading the data

phage_codon = read.table("02.Codon_Composition_Phage.dat", header = F)
colnames(phage_codon) = c("taxid","phage","codon","AA",
                          "pFraction","pFreq","pNumber")

host_codon = read.table("11.reformed_Codon_Composition_all_Bacteria_Ann.dat",header = F)
colnames(host_codon) = c("host","genome_no","codon","AA",
                         "hFraction","hFreq","hNumber")

host_taxid = read.table("07.BacHost_taxid.list", header = F)
colnames(host_taxid) = c("taxid","host")
host_taxid$taxid = as.factor(host_taxid$taxid)

host_merged = merge(host_codon, host_taxid, by = "host")
host_merged$taxid = as.factor(host_merged$taxid)
phage_codon$taxid = as.factor(phage_codon$taxid)

host_phage = merge(host_merged, phage_codon, by = c("taxid","codon","AA"))

ptrna = read.table("01.raw_Phage_tRNA.dat", header = F, sep = "\t")
colnames(ptrna) = c("taxid","phage","amino_acid","anti_codon")
ptrna$codon = toupper(chartr("ATGC", "TACG", toupper(ptrna$anti_codon)))


#htrna = read.table("09.bac_host_tRNA.dat", header = F, sep = "\t")
#colnames(htrna) = c("taxid","host","genome_no","amino_acid","anti_codon")
#htrna$codon = toupper(chartr("ATGC", "TACG", toupper(htrna$anti_codon)))

codon_to_aa = as.data.frame(table(phage_codon$codon,phage_codon$AA))
codon_to_aa[codon_to_aa == 0] = NA
codon_to_aa = na.omit(codon_to_aa)
codon_to_aa$Freq = NULL
colnames(codon_to_aa) = c("codon","AA")

trna_phage = merge(ptrna, codon_to_aa, by = "codon")
trna_phage$amino_acid = NULL
trna_phage$anti_codon = NULL
tphage = as.data.frame(table(trna_phage$codon, trna_phage$taxid,
                             trna_phage$phage, trna_phage$AA))
tphage[tphage == 0] = NA
tphage = na.omit(tphage)
colnames(tphage) = c("codon","taxid","phage","AA","number")
tphage$ptRNA = "Presence"




#trna_host = merge(htrna, codon_to_aa, by = "codon")
#trna_host$amino_acid = NULL
#trna_host$anti_codon = NULL
#trna_host$htRNA = "Presence"

#compiled1 = merge(trna_host, host_phage, by = c("taxid","host","genome_no",
#                                                "AA","codon"))


#phage_codon_trna[is.na(phage_codon_trna)] = "Absence"

ptrna_codon = merge(tphage, phage_codon, by = c("taxid","phage","AA","codon"), all = T)

#ptrna_codon[is.na(ptrna_codon$ptRNA)] = "Absence"


################################################################
################################################################
#             Reloding the reformed data




