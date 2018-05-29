#bin/R

#targeting the amino acid use - not codon
#
#setwd("/media/molgen/scan/10.experiment/05.tRNA_Check/05.Codon_Bias/04.Compiled/")
#working on my laptop
#changing the work_dir

setwd("/home/tmtareq/05.tRNA_check/05.Codon_Bias/04.Compiled/")


phage_codon = read.table("02.Codon_Composition_Phage.dat", header = F)
colnames(phage_codon) = c("taxid","phage","codon","AA","Fraction","Freq","Number")


library(dplyr)
library(moments)
library(ggplot2)
library(reshape2)
#host_codon = read.table("03.Codon_Composition_Bacteria_NoAnn.dat",header = F)
#colnames(host_codon) = c("taxid","host","genome_no","codon","AA",
#                         "Fraction","Freq","Number")


#       filtering host to reduce the biasness of AA usage
#       as the number of genome is increasing sd
#       reducing the the number to 1

#working on all the bacteria host to find out the relative bias in amino acid use


host_codon = read.table("11.reformed_Codon_Composition_all_Bacteria_Ann.dat")
colnames(host_codon) = c("host_name","host_genome_no", "codon", "aa", 
                         "host_codon_freq", "host_aa_fraction", "host_aa_count")



ggplot(data = host_codon, aes(x = host_aa_fraction, y = host_aa_count)) + 
  geom_point()


host_codon1 = host_codon %>% 
  group_by(host_name, host_genome_no) %>%
  summarise(sum(host_aa_fraction), sum(host_aa_count))
remove(host_codon1)


ggplot(data = host_codon, aes(aa, host_aa_fraction, fill = aa, col = codon)) + 
  geom_boxplot() 

#group_by(symbol, date_by_month) %>% summarize





host_aa_mean = host_codon %>% group_by(host, AA) %>% 
        summarise(mean(Number), max(Number), min(Number))

host_aa_total = host_codon %>% 
        group_by(taxid, host) %>%
        summarise(sum(Freq)) 

host_freq = as.data.frame(table(host_codon$host))
host_freq[host_freq == 0] = NA
host_freq = na.omit(host_freq)
colnames(host_freq) = c("host","occurance")
host_freq$count = host_freq$occurance/64

host_check = merge(host_aa_total, host_freq, by = "host")
host_check$norm_total_freq = host_check$`sum(Freq)`/host_check$count


#remove(new_dat)
#colnames(host_codon_mean) = c("taxid","host","codon","host_fraction")


host_codon_1 = subset(host_codon, genome_no == 1, 
                    select = c("taxid","host","AA","Freq"))

host_AA_total = host_codon_1 %>%
        group_by(taxid, host, AA) %>%
        summarise(sum(Freq))

View(phage_codon)

phage_AA_total = phage_codon %>%
        group_by(taxid, phage, AA) %>%
        summarise(sum(Freq))

#       loading the trna data
#

trna = read.table(file = "01.raw_Phage_tRNA.dat", sep = "\t")
colnames(trna) = c("taxid","phage","amino_acid","anti_codon")

trna$anti_codon1 = toupper(trna$anti_codon)
trna$codon = chartr("ATGC", "TACG", ac)

codon_to_aa = as.data.frame(table(phage_codon$codon,phage_codon$AA))
#df[df == 0] <- NA
codon_to_aa[codon_to_aa == 0] = NA
codon_to_aa = na.omit(codon_to_aa)
codon_to_aa$Freq = NULL
colnames(codon_to_aa) = c("codon","AA")

codon_trna = merge(trna, codon_to_aa, by = "codon")
codon_trna_1 = codon_trna %>%
        transmute(taxid, phage, AA, codon)
codon_trna_1$tRNA = "Presence"


codon_trna_sum = as.data.frame(table(codon_trna_1$taxid, 
                                     codon_trna_1$phage,
                                     codon_trna_1$AA))
codon_trna_sum[codon_trna_sum == 0] = NA
codon_trna_sum = na.omit(codon_trna_sum)


 