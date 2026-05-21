#!/usr/bin/env Rscript

library(dplyr)
library(reshape2)
library(stringr)
library(tidyr)
library("optparse")

option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# read in compiled table
o1<-read.table(opt$file, sep="\t", header=F)
colnames(o1) <- c("tax","ID","count")

# make table wide
o2<-dcast(o1, tax ~ ID)
o2[is.na(o2)] <- 0

# separate taxonomic levels
o3<- o2 %>% separate(tax,
           into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";")

# subset by domain
a<- o3 %>% filter(Domain == "Archaea")
e<- o3 %>% filter(Domain == "Eukaryota")
b<- o3 %>% filter(Domain == "Bacteria")

# write to directory
write.table(a, "archaea.ALL_count_table_L7.txt", sep="\t", quote = FALSE)
write.table(b, "bacteria.ALL_count_table_L7.txt", sep="\t", quote = FALSE)
write.table(e, "eukaryota.ALL_table_L7.txt", sep="\t", quote = FALSE)
