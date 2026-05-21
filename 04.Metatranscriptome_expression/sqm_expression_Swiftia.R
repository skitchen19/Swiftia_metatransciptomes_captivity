#####################
## SQMTools + DEG  ##
#####################

# load R packages
#if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager")}
#BiocManager::install("SQMtools")

library(SQMtools)
library(DESeq2)
library(stringr)
library("tibble")
library("ggplot2")
library(ggrepel)
library(PairedData)
library(variancePartition)
library("edgeR")
library("BiocParallel")
library(EnhancedVolcano)
library(dplyr)
library(ggpubr)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library("sva")
library(fgsea)
library("vegan")
library(viridis)
library(data.table)
library(pathview)
library(KEGGREST)

###############################
# load in SqueezeMeta results #
###############################

SQM <- loadSQM("I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/PlaceSamplesHere.zip")

# get taxa per transcript orf
tax <-as.data.frame(SQM$orfs$tax) %>% rownames_to_column()

# average coverage for each transcript by year
cov<-as.data.frame(SQM$orfs$cov) %>% rownames_to_column() %>% summarize(y0=rowMeans(.[,2:13]), y1=rowMeans(.[,14:25]), rowname=rowname)

# average abundance for each transcript by year
abun<-as.data.frame(SQM$orfs$abun) %>% rownames_to_column() %>% summarize(y0=rowMeans(.[,c(2:8,10:13)]), y1=rowMeans(.[,c(14:20,22:25)]), rowname=rowname)

# table with GC percent, join above tables, remove rRNA and tRNA counts
tab <-as.data.frame(SQM$orfs$table) %>% filter(Method != "barrnap") %>% filter(Method != "Aragorn") %>% rownames_to_column() %>%
  left_join(tax,by="rowname") %>% left_join(cov,by="rowname") 

# just rRNA counts
taxC <-as.data.frame(SQM$contigs$tax) %>% rownames_to_column()
tabrRNA <-as.data.frame(SQM$orfs$table) %>% filter(Method == "barrnap") %>% rownames_to_column() %>%
  left_join(taxC,by=c("Contig ID"="rowname")) %>% left_join(cov,by="rowname") 


#####################################
# percentage of transcripts by taxa #
#####################################

taxPercent <-tab %>% select(superkingdom, phylum, `Contig ID`) %>% filter(superkingdom != "Unclassified") %>% 
  mutate(superkingdom = ifelse(phylum == "Cnidaria", paste0(superkingdom,"-Cnidaria"), superkingdom)) %>% group_by(superkingdom) %>% summarise(percent = n() / nrow(.) * 100, count= n())

# pie chart of transcript taxonomy
ggplot(taxPercent, aes(x = "", y = percent, fill = superkingdom)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Removes background grid and axes
  scale_fill_manual(values=c(Bacteria ="#00d5a0", `Eukaryota-Cnidaria` ="#1389b1", `Eukaryota` ="#1389b150",Viruses ="#ee476e", Archaea="#FDA32F")) +
  theme_classic() +
  geom_text(aes(label = paste0(round(percent), "%")), position = position_stack(vjust = 0.5)) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))+
  labs(x = NULL, y = NULL, fill = NULL)

pie(taxPercent$percent, labels = taxPercent$superkingdom,
    col = rainbow(length(top_ten_descend$Population)), )
  

####################
# load in metadata #
####################

metadata <- read.table("I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/metadata_swiftia.txt",sep = "\t", header = T)
metadata$Sample_ID<-str_replace(metadata$Sample_ID, pattern="_rep", replacement="")

# remove sample R due to failed year 1 results
metadata2 <- metadata %>% dplyr::filter(Alt_ID !="R") %>% dplyr::filter(sample_type =="sample") %>% column_to_rownames("Sample_ID")

# set the different variables as factors
metadata2$year<-as.factor(metadata2$year)
metadata2$Alt_ID<-as.factor(metadata2$Alt_ID)
metadata2$location<-as.factor(metadata2$location)

# check the structure of the metadata file
str(metadata2)

#################
# apicomplexans #
#################

project.api <- subsetTax(SQM, "phylum", "Apicomplexa")
plotTaxonomy(project.api,"class",count='percent',rescale = F)

summary(project.api)

options(scipen=999)

# extract counts table
api<-as.data.frame(project.api$orfs$abund) %>% summarize(across(where(is.numeric), sum)) 

api<-reshape2::melt(api) %>% left_join(metadata, by=c("variable"="Sample_ID")) %>%
  dplyr::filter(Alt_ID !="R")
api$year <- as.factor(api$year)

a_my_comparisons <- list( c("0", "1"))

ggplot(api, aes(x=year, y=value, fill=year))+
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitter(w=0.1, h=0), shape=21) + # Jitter points for visibility
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif") +
  ylab("Plastid Read Count") +
  theme_classic() +scale_fill_manual(values=c(`0`="#FED068",`1`="#8540A0")) 

############################
# bacterial community data #
############################

project.bac <- subsetTax(SQM, "superkingdom", "Bacteria")

# extract counts table for dream
bacCount<-as.data.frame(project.bac$orfs$abund)

dim(bacCount)

# remove sample "R"
bacCount<-bacCount[,c(-8,-20)]

# remove rRNA read counts
btab <- tab %>% filter(superkingdom == "Bacteria")

bacCount <- bacCount %>% filter(rownames(bacCount) %in% btab$rowname)
dim(bacCount) # number of genes left

# minimum 10 read counts per row in three samples
min_reads = 10
min_samples = 3

# filter transcripts
keep_genes <- rowSums(bacCount >= min_reads) >= min_samples

bacCount <- bacCount[keep_genes,]
dim(bacCount) # number of genes left

# extract table of taxonomy
btax<-as.data.frame(project.bac$orfs$tax) %>% rownames_to_column()

bfun<-as.data.frame(project.bac$orfs$table) %>% rownames_to_column() %>% 
  mutate(`KEGG ID`=str_replace(`KEGG ID`, pattern="\\*", replacement=""),
         `COG ID`=str_replace(`COG ID`, pattern="\\*", replacement=""))


###################################
## Principal Components Analysis ##
###################################

# Standard usage of limma/voom
dge <- DGEList(bacCount)

dge <- calcNormFactors(dge)

# log counts per million
cpm_log <- cpm(dge, log = TRUE)

# top 500 genes, used by default in deseq2
rv <- rowVars(cpm_log)
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]

# run PCA
pca <- prcomp(t(cpm_log[select, ]), scale. = TRUE)
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2),1)

# set up table for plotting PCA
intgroup.df <- metadata2[,c(2,4)]
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                intgroup.df, libsize=dge$samples$norm.factors, names=metadata2$Alt_ID)

# plot PCA
ggplot(d, aes(PC1, PC2, color=year, shape=location, fill=year)) + 
  theme_bw()+
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c(`0`="#FED068",`1`="#8540A0"))+
  scale_fill_manual(values=c(`0`="#FED06880",`1`="#8540A080")) +
  scale_shape_manual(values=c(24,21,22))+
  geom_text_repel(d,mapping = aes(label = names), size = 4) +
  stat_ellipse(linetype = 2,alpha=0.2, geom = "polygon", aes(group=year))

# This reduces the data to 5 uncorrelated variables that explain most variance
pc_scores <- pca$x[, 1:5]

# Perform PERMANOVA using Euclidean distance on PCA scores
permanova_res <- adonis2(pc_scores ~ d$year , method = "euclidean")

print(permanova_res)

absDat <- read.table("I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/absolute.abundance.csv",sep = ",", header = T) %>%
  left_join(metadata, by=c("customer_label"="Sample_ID")) %>% filter( year != "N")

metadata3 <- cbind(metadata2, genome_copies=absDat$DNA_ng_per_ul.[c(-8,-20)])

#########################
## Analysis with DREAM ##
#########################

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(2, "SOCK", progressbar = TRUE, exportglobals = F)

# The variable to be tested must be a fixed effect
form <- ~ year + genome_copies + (1 | Alt_ID) + (1 | location)

# estimate weights using linear mixed model of dream
#vobjDream <- voomWithDreamWeights(dge, form, metadata2)
vobjDream <- voomWithDreamWeights(dge, form, metadata3, BPPARAM =param)

# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation
fitmm <- dream(vobjDream, form, metadata3)
fitmm <- eBayes(fitmm)

dim(dge)
dim(fitmm)

# Get results of hypothesis test on coefficients of interest
tabBac<-topTable(fitmm, coef = "year1", number = nrow(fitmm)) %>% rownames_to_column() %>% 
  left_join(bfun[,c(1,10:17)], by=c("rowname"="rowname")) %>% left_join(btax, by=c("rowname"="rowname")) %>%
  mutate(direction=ifelse(`adj.P.Val` < 0.05 & logFC > 0, "captive",
                          ifelse(`adj.P.Val` < 0.05 & logFC < 0, "wild", "ns")))

#write.table(tabBac,"I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/dream_Year_bacterial_20260520.txt",sep="\t", quote=F)

tabBac2 <- tabBac

tabBac2$genus<-ifelse(tabBac2$genus %in% c("BD1-7 clade", "Mycoplasma", "Unclassified Mycoplasmatales","Unclassified Mycoplasmatales","Unclassified Mycoplasmatota",
                                           "Unclassified Mollicutes","Unclassified Mycoplasmataceae",
                                           "Vibrio", "Unclassified Vibrionaceae",
                                           "Mycobacterium", "Methylobacterium", "FS140-16B-02 marine group",
                                            "Neptuniibacter","Pseudoalteromonas","Endozoicomonas", "Marivibrio",
                                           "Unclassified Oceanospirillaceae", "Unclassified Oceanospirillales",
                                           "Cutibacterium","Unclassified Propionibacteriaceae","Unclassified Actinomycetes"), tabBac2$genus, "other")

# volcano plot
p<-ggplot(tabBac2 %>% filter(genus == "other"), aes(x=logFC, y=-log(adj.P.Val)))+ theme_bw() +
  geom_point(aes(fill=direction),alpha = 0.4, shape=21, size =1, color="lightgrey") +
  xlab('Captive/Wild Expression (Log2)') + ylab(expression(-log[10]*" (P-value)")) +
  geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash') +
  geom_hline(yintercept = -log(0.05), colour = 'black', linetype = 'longdash') 
 
p + geom_point(data = tabBac2 %>% filter(genus != "other") %>% filter(adj.P.Val <= 0.05),
               aes(x=logFC, y=-log(adj.P.Val), fill=genus), alpha=0.9, shape=21, size=2, inherit.aes = F, color="grey25") +
  scale_fill_manual(values = c(`ns`="lightgrey",`BD1-7 clade`="#14655E", 
                                `Mycoplasma`="#B9DFBC", `Unclassified Mycoplasmatales`="#B9DFBC",
                               `Unclassified Mycoplasmatales`="#B9DFBC",`Unclassified Mycoplasmatota`="#B9DFBC",
                               `Unclassified Mollicutes`="#B9DFBC",`Unclassified Mycoplasmataceae`="#B9DFBC",
                                Mycobacterium="#1A936F", Pseudoalteromonas="#49AB7C",
                               Neptuniibacter="#3BC95A", `Unclassified Oceanospirillaceae`="#3BC95A",`Unclassified Oceanospirillales`="#3BC95A",
                               Vibrio="#14847A",`Unclassified Vibrionaceae`="#14847A",
                               "FS140-16B-02 marine group"="#BFEDC9", 
                               Methylobacterium= "#DDE2C9", Endozoicomonas="#F3E9D2",  
                               Marivibrio="#F7F3EC", `other`="darkgrey", `Unclassified Actinomycetes`="darkorange",
                               Cutibacterium="darkorange",`Unclassified Propionibacteriaceae`="darkorange")) 

# by environment
ggplot(tabBac2 %>% filter(genus == "other"), aes(x=logFC, y=-log(adj.P.Val)))+ theme_bw() +
  geom_point(aes(color=direction),alpha = 0.4, shape=19, size =2) +
  xlab('Captive/Wild Expression (Log2)') + ylab(expression(-log[10]*" (P-value)")) +
  geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash') +
  geom_hline(yintercept = -log(0.05), colour = 'black', linetype = 'longdash') +
  scale_color_manual(values=c(`captive`="#8540A0",`wild`="#FED068", `ns`="lightgrey"))

################################
# KEGG overrepresentation test #
################################

# KEGG over-representation in wild
kk <- enrichKEGG(gene         = tabBac$`KEGG ID`[tabBac$direction == "wild"],
                 organism     = 'ko',
                 pvalueCutoff = 0.05)
head(kk)


# KEGG over-representation in captive
kk <- enrichKEGG(gene         = tabBac$`KEGG ID`[tabBac$direction == "captive"],
                 organism     = 'ko',
                 pvalueCutoff = 0.05)
head(kk)

as.data.frame(kk)


###############################
# BD1-7-like genome reference #
###############################

# differential expression of bacteria
bdT<- read.table("I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/bacterial_genome_assemblies/13.bd17_expression.orftable",sep = "\t", header = T, quote="") %>%
  filter(Molecule=="CDS")
rownames(bdT) <- bdT$ORF.ID
colnames(bdT) <- str_replace(colnames(bdT), pattern="Raw.read.count.", replacement="")

# subset to column with read counts
bdT2<-bdT[,c(65:88)]
dim(bdT2)

# remove sample "R"
bdT2<-bdT2[,c(-8,-20)]

# minimum 10 read counts per row
min_reads = 10
min_samples = 1

keep_genes <- rowSums(bdT2 >= min_reads) >= min_samples
#keep_genes <- rowSums(bdT2 >= min_reads) 

bdT2 <- bdT2[keep_genes,]
dim(bdT2) # number of genes left

# number of genes with at least one read
colSums(bdT2>=10)

###################################
## Principal Components Analysis ##
###################################
dge <- DGEList(bdT2)

dge <- calcNormFactors(dge)

# log counts per million
cpm_log_BD <- cpm(dge, log = TRUE)

# top 500 genes, used by default in deseq2
rv <- rowVars(cpm_log_BD)
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(cpm_log_BD[select, ]),scale. = TRUE)
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2),1)

intgroup.df <- metadata2[,c(2,4)]
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                intgroup.df, names = metadata2$Alt_ID, libsize=dge$samples$norm.factors)

# plot PCA
ggplot(d, aes(PC1, PC2, color=year, shape=location, fill=year)) + 
  theme_bw()+
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  scale_color_manual(values=c(`0`="#FED068",`1`="#8540A0"))+
  scale_fill_manual(values=c(`0`="#FED06880",`1`="#8540A080")) +
  scale_shape_manual(values=c(24,21,22))+
  geom_text_repel(d,mapping = aes(label = names), size = 4)+
  stat_ellipse( linetype = 2,alpha=0.2, geom = "polygon",aes(fill=year, group=year))

#########################
## Analysis with DREAM ##
#########################

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(2, "SOCK", progressbar = TRUE, exportglobals = F)

# The variable to be tested must be a fixed effect
form <- ~ year + (1 | Alt_ID) + (1 | location)

# estimate weights using linear mixed model of dream
#vobjDream <- voomWithDreamWeights(dge, form, metadata2)
vobjDream <- voomWithDreamWeights(dge, form, metadata2, BPPARAM =param)

# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation
fitmm <- dream(vobjDream, form, metadata2)
fitmm <- eBayes(fitmm)

dim(dge)
dim(fitmm)


# Get results of hypothesis test on coefficients of interest
tabBD<-topTable(fitmm, coef = "year1", number = nrow(fitmm)) %>% rownames_to_column() %>% 
  left_join(bdT[,c(1,8:16)], by=c("rowname"="ORF.ID")) %>% 
  mutate(direction=ifelse(`adj.P.Val` < 0.05 & logFC > 0, "captive",
                          ifelse(`adj.P.Val` < 0.05 & logFC < 0, "wild", "ns")))%>% 
  mutate(`KEGG.ID`=str_replace(`KEGG.ID`, pattern="\\*", replacement=""),
         `COG.ID`=str_replace(`COG.ID`, pattern="\\*", replacement=""))

#write.table(tabBD,"dream_Year_BD17_20260421.txt",sep="\t", quote=F)

# volcano plot
ggplot(tabBD, aes(x=logFC, y=-log(adj.P.Val)))+ theme_bw() +
  geom_point(aes(color=direction),alpha = 0.7, shape=19) +
  scale_color_manual(values=c(`captive`="#FED068",`wild`="#8540A0", `ns`="grey"))+
  xlab('Captive/Wild Expression (Log2)') + ylab(expression(-log[10]*" (P-value)")) +
  geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash') +
  geom_hline(yintercept = -log(0.05), colour = 'black', linetype = 'longdash')



##########################
# Mycoplasma transcripts #
##########################

# differential expression of bacteria
myT <- read.table("I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/bacterial_genome_assemblies/13.mycoplasma_expression_v2.orftable",sep = "\t", header = T, quote="") %>%
  filter(Molecule=="CDS")

myT.rRNA <- read.table("I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/bacterial_genome_assemblies/13.mycoplasma_expression_v2.orftable",sep = "\t", header = T, quote="") %>%
  filter(Molecule=="rRNA")

rownames(myT) <- myT$ORF.ID
colnames(myT) <- str_replace(colnames(myT), pattern="Raw.read.count.", replacement="")

# subset to column with read counts
myT2<-myT[,c(65:88)]
dim(myT2)

# keep samples that have Mycoplasma present
myT2<-myT2[,c(2,3,5,9,10,14,15,17,21,22)]

# minimum 10 read counts per row
min_reads = 10
min_samples = 1

keep_genes <- rowSums(myT2 >= min_reads) >= min_samples

myT2 <- myT2[keep_genes,]
dim(myT2) # number of genes left

# number of genes with at least one read
colSums(myT2>=10)


###################################
## Principal Components Analysis ##
###################################
dge <- DGEList(myT2)

dge <- calcNormFactors(dge)

# log counts per million
cpm_log_MY <- cpm(dge, log = TRUE)

# top 500 genes, used by default in deseq2
rv <- rowVars(cpm_log_MY)
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(cpm_log_MY[select, ]),scale. = TRUE)
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2),1)

intgroup.df <- metadata2[c(2,3,5,8,9,13,14,16,19,20),c(2,4)]
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                intgroup.df, names = metadata2$Alt_ID[c(2,3,5,8,9,13,14,16,19,20)], libsize=dge$samples$norm.factors)

# plot PCA
ggplot(d, aes(PC1, PC2, color=year, shape=location, fill=year)) + 
  theme_bw()+
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values=c(`0`="#FED068",`1`="#8540A0"))+
  scale_fill_manual(values=c(`0`="#FED06880",`1`="#8540A080")) +
  scale_shape_manual(values=c(21,22))+
  geom_text_repel(d,mapping = aes(label = names), size = 4)+
  stat_ellipse( linetype = 2,alpha=0.2, geom = "polygon",aes(fill=year, group=year))

#########################
## Analysis with DREAM ##
#########################

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(2, "SOCK", progressbar = TRUE, exportglobals = F)

# The variable to be tested must be a fixed effect
form <- ~ year + (1 | Alt_ID) + (1 | location)

# estimate weights using linear mixed model of dream
#vobjDream <- voomWithDreamWeights(dge, form, metadata2)
vobjDream <- voomWithDreamWeights(dge, form, metadata2[c(2,3,5,8,9,13,14,16,19,20),], BPPARAM =param)

# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation
fitmm <- dream(vobjDream, form, metadata2[c(2,3,5,8,9,13,14,16,19,20),])
fitmm <- eBayes(fitmm)

dim(dge)
dim(fitmm)

# Get results of hypothesis test on coefficients of interest
tabMY<-topTable(fitmm, coef = "year1", number = nrow(fitmm)) %>% rownames_to_column() %>% 
  left_join(myT[,c(1,8:16)], by=c("rowname"="ORF.ID")) %>% 
  mutate(direction=ifelse(`adj.P.Val` < 0.05 & logFC > 0, "captive",
                          ifelse(`adj.P.Val` < 0.05 & logFC < 0, "wild", "ns"))) %>% 
  mutate(`KEGG.ID`=str_replace(`KEGG.ID`, pattern="\\*", replacement=""),
         `COG.ID`=str_replace(`COG.ID`, pattern="\\*", replacement=""))

write.table(tabMY,"dream_Year_mycoplasma_v2_20260521.txt",sep="\t", quote=F)

################################
# KEGG overrepresentation test #
################################

# KEGG over-representation in wild
kk <- enrichKEGG(gene         = tabMY$`KEGG.ID`[tabMY$direction == "wild"],
                 organism     = 'ko',
                 pvalueCutoff = 0.05)
head(kk)


# KEGG over-representation in captive
kk <- enrichKEGG(gene         = tabMY$`KEGG.ID`[tabMY$direction == "captive"],
                 organism     = 'ko',
                 pvalueCutoff = 0.05)
head(kk)

as.data.frame(kk)

#######################################
# percent of 16S reads to total reads #
#######################################

## % of total reads
cntMY <- read.csv("read_count_mycoplasma.csv",header=T) %>% 
  left_join(absDat %>% select(Alt_ID,year,location, customer_label, DNA_ng_per_ul., genome_copies_per_ul.), by=c("Sample"="customer_label"))%>%
  mutate(abs=percent_total*genome_copies_per_ul.)

cntMY_sub <- cntMY %>% filter(Alt_ID %in% c("C","D","F","U","AB"))

cntMY_sub  %>%
  group_by(Gen) %>%
  reframe(
    compare_means(abs ~ year, data = cur_data(), method = "wilcox.test", paired=TRUE)
  )

ggplot(cntMY_sub %>% filter(Gen == "PlaceSamplesHere_273674_42-1812") , aes(x=year, y=abs, fill=year)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(w=0.1, h=0), shape=21) + # Jitter points for visibility
  #geom_text_repel(aes(label=rownames(plot_data))) + # Add sample labels
  theme_bw() +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons , label = "p.signif",paired=TRUE) +
  theme(plot.title = element_text(hjust=0.5))+
  scale_fill_manual(values=c(`0`="#FED068",`1`="#8540A0")) 


########################################
## combined heat map of microbial seqs #
########################################

# extract counts table for dream
bacCount<-as.data.frame(project.bac$orfs$abund)

# remove sample "R"
bacCount<-bacCount[,c(-8,-20)]

# remove rRNA read counts
btab <- tab %>% filter(superkingdom == "Bacteria")

bacCount <- bacCount %>% filter(rownames(bacCount) %in% btab$rowname)

dge <- DGEList(bacCount)

dge <- calcNormFactors(dge)

# log counts per million
cpm_log <- cpm(dge, log = TRUE)

cpm <- as.data.frame(cpm_log) %>% rownames_to_column() %>% summarize(y0=rowMeans(.[,c(2:12)]), y1=rowMeans(.[,c(13:23)]), rowname=rowname) %>%
  left_join(bfun[,c(1,10:17)], by=c("rowname"="rowname")) %>% left_join(btax, by=c("rowname"="rowname")) %>%
  left_join(tabBac %>% select(rowname, logFC,adj.P.Val, direction)) %>% filter(`KEGG ID` != "") %>% mutate (group = "community") %>%
  select(rowname, `KEGG ID`, KEGGFUN,logFC,adj.P.Val, `Gene name`, genus, direction, group, y0, y1) %>%
  mutate(s.y0 = as.numeric(scale(y0)), s.y1 = as.numeric(scale(y1)))

colnames(cpm) <-c( "rowname", "KEGG.ID" ,"KEGGFUN","logFC" ,"adj.P.Val","Gene.name",  "genus" ,"direction","group" ,"y0", "y1", "s.y0", "s.y1" )

# Mycoplasma
myT2<-myT[,c(65:88)]

# keep samples that have Mycoplasma present
myT2<-myT2[,c(2,3,5,9,10,14,15,17,21,22)]

dgeMY <- DGEList(myT2)

dgeMY <- calcNormFactors(dgeMY)

# log counts per million
cpm_log_MY <- cpm(dgeMY, log = TRUE)

cpmMY <- as.data.frame(cpm_log_MY) %>% rownames_to_column() %>% summarize(y0=rowMeans(.[,c(2:6)]), y1=rowMeans(.[,c(7:11)]), rowname=rowname)%>%
  left_join(myT[,c(1,8,10:16)], by=c("rowname"="ORF.ID")) %>%
  left_join(tabMY %>% select(rowname,`KEGG.ID`, KEGGFUN,logFC, adj.P.Val, `Gene.name`, direction), by = c("rowname")) %>%
  filter(`KEGG.ID.x` != "") %>% mutate (group = "Mycoplasma", genus = "Mycoplasma") %>%
  mutate(`KEGG.ID.x`=str_replace(`KEGG.ID.x`, pattern="\\*", replacement="")) %>%
  select(rowname, `KEGG.ID.x`, KEGGFUN.x, logFC, adj.P.Val, `Gene.name.x`, genus, direction, group, y0, y1) %>%
  mutate(s.y0 = as.numeric(scale(y0)), s.y1 = as.numeric(scale(y1)))

colnames(cpmMY) <-c( "rowname", "KEGG.ID" ,"KEGGFUN","logFC" ,"adj.P.Val","Gene.name",  "genus" ,"direction","group" ,"y0", "y1" , "s.y0", "s.y1")

# BD1-7 
bdT2<-bdT[,c(65:88)]

# remove sample "R"
bdT2<-bdT2[,c(-8,-20)]

dgeBD <- DGEList(bdT2)

dgeBD <- calcNormFactors(dgeBD)

# log counts per million
cpm_log_BD <- cpm(dgeBD, log = TRUE)

cpmBD <- as.data.frame(cpm_log_BD) %>% rownames_to_column() %>% summarize(y0=rowMeans(.[,c(2:6)]), y1=rowMeans(.[,c(7:11)]), rowname=rowname)%>%
  left_join(bdT[,c(1,8,10:16)], by=c("rowname"="ORF.ID")) %>%
  left_join(tabBD %>% select(rowname,`KEGG.ID`,KEGGFUN,logFC,adj.P.Val, `Gene.name`, direction), by="rowname") %>%
  filter(`KEGG.ID.x` != "") %>% mutate (group = "BD1-7", genus = "BD1-7 clade") %>%
  mutate(`KEGG.ID.x`=str_replace(`KEGG.ID.x`, pattern="\\*", replacement="")) %>%
  select(rowname, `KEGG.ID.x`, KEGGFUN.x,logFC,adj.P.Val, `Gene.name.x`, genus, direction, group, y0, y1) %>%
  mutate(s.y0 = as.numeric(scale(y0)), s.y1 = as.numeric(scale(y1)))


colnames(cpmBD) <-c( "rowname", "KEGG.ID" ,"KEGGFUN","logFC" ,"adj.P.Val","Gene.name",  "genus" ,"direction","group" ,"y0", "y1" , "s.y0", "s.y1")

# join all tables together
all_bac<-rbind(cpm, cpmMY,cpmBD)

genusList <- c("BD1-7 clade", "Mycoplasma", "Unclassified Mycoplasmatales","Unclassified Mycoplasmatales","Unclassified Mycoplasmatota",
               "Unclassified Mollicutes","Unclassified Mycoplasmataceae",
               "Vibrio", "Unclassified Vibrionaceae",
               "Mycobacterium", "Methylobacterium", "FS140-16B-02 marine group",
               "Neptuniibacter","Pseudoalteromonas","Endozoicomonas", "Marivibrio",
               "Unclassified Oceanospirillaceae", "Unclassified Oceanospirillales")

all_bac_sig <- all_bac %>% 
  mutate(genus= case_when(group == "community" ~ ifelse(genus %in% c("Mycoplasma", "Unclassified Mycoplasmatales","Unclassified Mycoplasmatales","Unclassified Mycoplasmatota",
                                                                     "Unclassified Mollicutes","Unclassified Mycoplasmataceae"), "Mycoplasmatota",genus), TRUE ~ genus )) %>%
  mutate(genus= case_when(group == "community" ~ ifelse(genus %in% c("Unclassified Oceanospirillaceae", "Unclassified Oceanospirillales",
                                                                     "Neptuniibacter"), "Oceanospirillaceae",genus), TRUE ~ genus )) %>%
  mutate(genus= case_when(group == "community" ~ ifelse(genus %in% c("Vibrio", "Unclassified Vibrionaceae"), "Vibrionaceae",genus), TRUE ~ genus )) %>%
  mutate(genus= case_when(group == "community" ~ ifelse(genus %in% c("Cutibacterium", "Unclassified Propionibacteriaceae",
                                                                     "Unclassified Actinomycetes"), "Actinomycetes",genus), TRUE ~ genus )) %>%
  group_by(group) %>%
  filter(genus %in% c("Vibrionaceae", "Mycoplasma","BD1-7 clade","Mycoplasmatota","Oceanospirillaceae","Pseudoalteromonas","Actinomycetes")) %>%
  ungroup() %>%
  group_by(KEGG.ID) %>%
  filter(any(direction == "captive" | direction == "wild")) %>%
  ungroup() %>%
  group_by(KEGG.ID, genus, group) %>%
  arrange(adj.P.Val) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(diffExp = y1-y0) %>%
  mutate(KEGG_grp = ifelse(KEGG.ID %in% c("K00611","K00926", "K01478"), "arginine biosyn",
                           ifelse(KEGG.ID %in% c("K00114","K00134", "K00138","K01689"), "glycolysis", 
                                  ifelse(KEGG.ID %in% c("K03475","K02822","K02821","K03077", "K03079", "K02800", "K02768;K02769;K02770","K02799;K02800"), "Ascorbate_PTS",
                                         ifelse(KEGG.ID %in% c("K00239","K00404","K02108","K02109", "K02110","K02111","K02112","K02113","K02115","K00864","K01903",
                                                               "K00111","K00347","K02133"), "Oxphos",
                                                ifelse(KEGG.ID %in% c("K02863","K02864","K02867","K02871","K02874","K02876","K02878","K02881",
                                                                      "K02884","K02886","K02888","K02890","K02892","K02895","K02899","K02902",
                                                                      "K02904","K02906","K02907","K02913","K02926","K02931","K02933","K02935",
                                                                      "K02945","K02946","K02948","K02950","K02952","K02954","K02961","K02965",
                                                                      "K02967","K02982","K02986","K02988","K02992","K02994","K02996","K02358","K02357",
                                                                      "K02519","K00554","K02355","K02860"),"Translation", "other"))))))


all_bac_sig$group<-factor(all_bac_sig$group, levels=c("community", "Mycoplasma", "BD1-7"))
all_bac_sig$genus <-factor(all_bac_sig$genus, levels=c("Actinomycetes","Pseudoalteromonas","Vibrionaceae","Oceanospirillaceae", "Mycoplasmatota","Mycoplasma","BD1-7 clade"))


bsMelt <- melt(all_bac_sig %>% select(genus, KEGG.ID, KEGG_grp,s.y0, s.y1, group))

# final heat map of all taxa from the different data sets
ggplot(data = all_bac_sig %>% filter(KEGG_grp != "other"), mapping = aes(x = genus, y = KEGG.ID, fill = diffExp)) +
  geom_tile() +
  xlab(label = "Sample") +
  facet_grid(KEGG_grp~  group , scales = "free", space = "free",switch = "y")+
  scale_fill_viridis_c(option = "mako", direction = 1) +
  theme(strip.placement = "outside")+
  theme_bw()

####################
# Swifita analysis #
####################

project.cnid <- subsetTax(SQM, "phylum", "Cnidaria")

# extract counts table for dream
cnidCount<-as.data.frame(project.cnid$orfs$abund) 
dim(cnidCount) # number of genes left

# remove sample "R"
cnidCount<-cnidCount[,c(-8,-20)]

# remove rRNA read counts
ctab <- tab %>% filter(phylum == "Cnidaria")

cnidCount <- cnidCount %>% filter(rownames(cnidCount) %in% ctab$rowname)
dim(cnidCount) # number of genes left

# minimum 10 read counts per row
min_reads = 10
min_samples = 11

keep_genes <- rowSums(cnidCount >= min_reads) >= min_samples

cnidCount <- cnidCount[keep_genes,]
dim(cnidCount) # number of genes left

# extract table of functions
cfun<-as.data.frame(project.cnid$orfs$table) %>% rownames_to_column() %>% 
  mutate(`KEGG ID`=str_replace(`KEGG ID`, pattern="\\*", replacement=""),
         `COG ID`=str_replace(`COG ID`, pattern="\\*", replacement=""))

########################################
## Principal Components Analysis: ALL ##
########################################

# load in counts table
dge <- DGEList(cnidCount)

# normalization factor for each library
dge <- calcNormFactors(dge)

# extract log counts per million
cpm_log <- cpm(dge, log = TRUE)

# top 500 genes, used by default in deseq2
rv <- rowVars(cpm_log)
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(cpm_log[select, ]),scale. = TRUE)
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2),1)

# create groupings
intgroup.df <- metadata2[,c(2,4)]
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                intgroup.df, names = metadata2$Alt_ID)

# plot PCA
ggplot(d, aes(PC1, PC2, color=year, shape=location, fill=year)) + 
  theme_bw()+
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  scale_color_manual(values=c(`0`="#FED068",`1`="#8540A0"))+
  scale_fill_manual(values=c(`0`="#FED06880",`1`="#8540A080")) +
  scale_shape_manual(values=c(24,21,22))+
  geom_text_repel(d,mapping = aes(label = names), size = 4)

####################
## DREAM analysis ##
####################

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(2, "SOCK", progressbar = TRUE, exportglobals = F)

# The variable to be tested must be a fixed effect
form <- ~ year + (1 | Alt_ID) + (1 | location)

# estimate weights using linear mixed model of dream
#vobjDream <- voomWithDreamWeights(dge, form, metadata2)
vobjDream <- voomWithDreamWeights(dge, form, metadata2, BPPARAM = param)

# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation
fitmm <- dream(vobjDream, form, metadata2)
fitmm <- eBayes(fitmm)

# check contrasts tested
fitmm$contrasts

# Get results of hypothesis test on coefficients of interest
tabhost<-topTable(fitmm, coef = "year1", number = nrow(fitmm)) %>% rownames_to_column() %>% 
  left_join(cfun[,c(1,10:17)], by=c("rowname"="rowname")) %>%
  mutate(direction=ifelse(`adj.P.Val` < 0.05 & logFC > 0, "captive",
                          ifelse(`adj.P.Val` < 0.05 & logFC < 0, "wild", "ns")))

str(tabhost)

#write.table(tabhost,"G:/Shared drives/MDBC_data/NOAA_Swiftia_grafting/metatranscriptomics/expression/dream_Year_host_20250729.txt",sep="\t", quote=F, row.names=F)

tabhost<- read.table("G:/Shared drives/MDBC_data/NOAA_Swiftia_grafting/metatranscriptomics/expression/dream_Year_host_20250729.txt",sep="\t", quote="", header=T) %>%
  mutate(direction2=ifelse(`adj.P.Val` < 0.05 & logFC > 1, "captive",
                          ifelse(`adj.P.Val` < 0.05 & logFC < -1, "wild", "ns")))

########################################
## Principal Components Analysis: DEG ##
########################################

cnidCountSig <- cnidCount %>% filter(rownames(cnidCount) %in% tabhost$rowname[tabhost$direction2 != "ns"])
dim(cnidCountSig) # number of genes left

# load in counts table
dge <- DGEList(cnidCountSig)

# normalization factor for each library
dge <- calcNormFactors(dge)

# extract log counts per million
cpm_log <- cpm(dge, log = TRUE)

# top 500 genes, used by default in deseq2
rv <- rowVars(cpm_log)
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(cpm_log[select, ]),scale. = TRUE)
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2),1)

# create groupings
intgroup.df <- metadata2[,c(2,4)]
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                intgroup.df, names = metadata2$Alt_ID)

# plot PCA
ggplot(d, aes(PC1, PC2, color=year, shape=location, fill=year)) + 
  theme_bw()+
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  scale_color_manual(values=c(`0`="#FED068",`1`="#8540A0"))+
  scale_fill_manual(values=c(`0`="#FED06880",`1`="#8540A080")) +
  scale_shape_manual(values=c(24,21,22))+
  geom_text_repel(d,mapping = aes(label = names), size = 4)+
  stat_ellipse( linetype = 2,alpha=0.2, geom = "polygon",aes(fill=year, group=year))

#####################################
# KEGG gene set enrichment analysis #
#####################################

# KEGG over-representation in wild
tabhostSig <- tabhost %>% filter(direction2 != "ns") %>% filter( KEGG.ID != "")

list1 <- tabhostSig %>% filter(direction2 == "wild") %>% select(rowname, KEGG.ID) 
list2 <- tabhostSig %>% filter(direction2 == "captive") %>% select(rowname, KEGG.ID) 

# load in human kegg pathways to exclude
human <- read.table("I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/human_KEGG.txt",sep="\t", quote="", header=T)

#
input <- tabhost %>% select(KEGG.ID, logFC,adj.P.Val,z.std) %>% filter(KEGG.ID != "") %>% 
  mutate(rank= -log10(adj.P.Val)*sign(logFC)) %>% group_by(KEGG.ID) %>% summarize(z.stdAvg = mean(z.std)) %>% 
  ungroup() %>% dplyr::arrange(-z.stdAvg) 

original_gene_list <- input$z.stdAvg
names(original_gene_list) <- input$KEGG.ID
gene_list = sort(original_gene_list, decreasing = TRUE)

# load in human kegg pathways to exclude
human <- read.table("I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/human_KEGG.txt",sep="\t", quote="", header=T)

# kegg pathways
ko.pathways <- download_KEGG(species="ko")
term2gene <- ko.pathways$KEGGPATHID2EXTID %>% filter(!from %in% human$ID[human$group == "human"])
grouped_list <- split(term2gene$to, term2gene$from)
term2name <- ko.pathways$KEGGPATHID2NAME %>% filter(!from %in% human$ID[human$group == "human"])

# run GSEA 
fgseaRes <- fgseaMultilevel(pathways = grouped_list, 
                            stats    = gene_list,
                            eps      = 0.0,
                            minSize  = 10,
                            maxSize  = 200,
                            nPermSimple = 10000)

# download generic pathways
pathways <- keggList("pathway", "ko")

get_kegg_category <- function(path_id) {
  info <- keggGet(path_id)[[1]]
  
  if (!is.null(info$CLASS)) {
    # CLASS looks like: "Metabolism; Carbohydrate metabolism"
    classes <- unlist(strsplit(info$CLASS, "; "))
    
    data.frame(
      pathway = path_id,
      category = classes[1],
      subcategory = ifelse(length(classes) > 1, classes[2], NA),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      pathway = path_id,
      category = NA,
      subcategory = NA,
      stringsAsFactors = FALSE
    )
  }
}

# add categories to the pathways
kegg_df <- do.call(rbind, lapply(fgseaResFilt$pathway, get_kegg_category))
kegg_df

# remove pathways that are human, diseaes or otherwise not related to cnidarians
fgseaResFilt <- fgseaRes %>% left_join(term2name, by= c("pathway"="from")) %>% 
  filter(padj < 0.05) %>%
  filter(!str_detect(to, "yeast")) %>%
  filter(!str_detect(to, "Microbial")) %>%
  filter(!pathway %in% human$ID[human$group == "bacteria"]) %>%
  left_join(kegg_df, by=c("pathway")) %>%
  mutate(category=ifelse(is.na(category),"Metabolism",category))

fgseaResFilt$category <- as.factor(fgseaResFilt$category )

# plot enriched pathways
ggplot(fgseaResFilt, aes(x=NES, y= reorder(to, -padj+NES), fill= NES))+ 
  geom_bar(stat = "identity") +
  scale_fill_gradient(low= "#FED068", high="#8540A0") +
  theme_bw() +
  xlab("Normalized Enrichment Score") + ylab("") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y")+
  theme(
    strip.background = element_blank(),
    panel.spacing = unit(0, "lines")
  )

# write the nested table results to file
fwrite(fgseaResFilt, file = "KEGG_gsea.txt", sep = "\t", quote = FALSE, sep2 = c("", " ", ""))

# Produce the native KEGG plot (PNG)
p1 <- pathview(gene.data= gene_list, pathway.id="04144", species = "ko", low= "#FED068", high="#8540A0", 
                node.sum="mean", limit = list(gene = 3), gene.idtype="kegg",same.layer = TRUE)
