###################
# Kraken2 Results #
###################

library(ggplot2)

# load in table from phyloflash
Ktable <- read.table(file = "I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/kraken_domain_summary_adj.tsv", sep = "\t", header = T) 

Ktable <- Ktable[1:24,-6]

Kmelt <- melt(Ktable)

ggplot(Kmelt, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c(Bacteria ="#00d5a0", Eukaryota ="#1389b1", Viruses ="#ee476e", Archaea="#FDA32F", Unclassified_adj="#083a4b")) +
  theme_bw() +  scale_y_continuous(expand = c(0,0)) + ylab("Fraction of mRNA Reads")
