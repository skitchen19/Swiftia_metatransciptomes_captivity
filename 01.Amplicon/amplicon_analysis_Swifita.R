#######################
## 16S Amplicon Data ##
#######################

# load R packages
library(dplyr)
library(ggplot2)
library("ggpubr")
library(PairedData)
library("tibble")
library("phyloseq")
library("decontam")
library(reshape2)
library(vegan)
library(ggplot2)
library(ggrepel)
library(stringr)
library(DESeq2)
library("MicEco")
library("ANCOMBC")

##########################
# create phyloseq object #
##########################

#phyloseq_object_old <- readRDS("C:/Users/Sheila's Comp/Desktop/swiftia_phyloseq.rds") #12.3.25
phyloseq_object <- readRDS("C:/Users/Sheila's Comp/Desktop/phyloseq_all_16Sdata.rds") #3.18.25

# change column names to correct taxonomic levels
colnames(tax_table(phyloseq_object)) <- c("Kingdom","Phylum","Class", "Order","Family","Genus", "Species")

# check top 10 lines
head(tax_table(phyloseq_object))

# load in metadata
metadata <- read.table("I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/metadata_swiftia.txt",sep = "\t", header = T)

# make names consistent
metadata$Sample_ID<-str_replace(metadata$Sample_ID, pattern="_rep", replacement="")

# remove controls
metadata2 <-metadata %>% dplyr::filter(Sample_ID !="Tank") %>% dplyr::filter(year !="N") %>% column_to_rownames("Sample_ID")

# set variables as factors
metadata2$year<-as.factor(metadata2$year)
metadata2$Alt_ID<-as.factor(metadata2$Alt_ID)
metadata2$location<-as.factor(metadata2$location)

# structure of data
str(metadata2)

# sample table for phyloseq
SAMPLE <- sample_data(metadata2)

# apicomplexan ASVs
apicomplexan_asvs <- c("ASV_1122","ASV_982","ASV_8","ASV_811")

# remove host/contaminant mitochondrial ASVs
mito_asvs <- c("ASV_700")

# Keep only bacterial OTUs, remove chloroplast, mitochondrial and apicomplexan ASVs
taxTAB <- as.data.frame(phyloseq_object@tax_table) %>% 
  filter(Kingdom == "Bacteria") %>%
  filter(Order != "Chloroplast") %>%
  dplyr::filter(.,!(rownames(.) %in% apicomplexan_asvs)) %>%
  dplyr::filter(.,!(rownames(.) %in% mito_asvs))

# make taxonomy table
TAX = tax_table(as.matrix(taxTAB))

# remove contaminants from ASV table
otuTAB<-as.data.frame(phyloseq_object@otu_table)
colnames(otuTAB) <- metadata$Sample_ID

otuTAB <-  otuTAB[rownames(otuTAB) %in% rownames(taxTAB), ] 

# rarefaction curve
rarecurve(t(otuTAB), step=50, cex=0.5)

# make otu table
OTU <- otu_table(otuTAB[,1:24],taxa_are_rows=TRUE)

# create phyloseq object
physeq <- phyloseq(OTU, TAX, SAMPLE)
physeq

####################
# check everything #
####################

sample_names(physeq)
rank_names(physeq)
sample_variables(physeq)

###########################
# contamination screening #
###########################
# frequency of contaminants
contamdf.freq <- isContaminant(physeq.prune, method="frequency", conc="DNA_conc")
head(contamdf.freq)

table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

plot_frequency(physeq, taxa_names(physeq.prune)[c(12)], conc="DNA_conc") + 
  xlab("DNA Concentration")

# prevalence of contaminants
sample_data(physeq)$is.neg <- sample_data(physeq)$year == "N"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg",threshold=0.5)
table(contamdf.prev$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$year == "0", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$year == "1", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

########################################
# remove ASV with only one observation #
########################################

physeq.prune = prune_taxa(taxa_sums(physeq) > 1, physeq)
physeq.prune

#############################
# alpha diversity estimates #
#############################

# plot alpha diversity
rich <- estimate_richness(physeq.prune, split=TRUE, measures=c("Chao1", "Shannon")) %>% 
  mutate(year=metadata2$year[1:24]) %>% filter (year !="N")

# kruskal wallis on diversity indicies by site (reef vs. crevice)
sppdiv_kw <- kruskal.test(Shannon ~ year, data = rich)
sppdiv_kw # significantly different

a_my_comparisons <- list( c("0", "1"))

a <-ggplot(rich, aes(x=year, y=Chao1, fill=year)) +
  geom_boxplot(aes(fill = year), show.legend = FALSE,outlier.shape = NA) +
  geom_point(position=position_jitter(w=0.1, h=0), shape=21) + # Jitter points for visibility
  theme_classic() +
  scale_fill_manual(values=c(`0`="#FED068",`1`="#8540A0")) +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons , label = "p.signif") +
  ylab("Estimated ASV Richness (Chao1)") +
  theme(axis.title.y = element_text(face = "bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.text.x = element_text(face = "bold", size=14),
        strip.background = element_blank())

b<-ggplot(rich, aes(x=year, y=Shannon, fill=year)) +
  geom_boxplot(aes(fill = year), show.legend = FALSE,outlier.shape = NA) +
  geom_point(position=position_jitter(w=0.1, h=0), shape=21) + # Jitter points for visibility
  theme_classic() +
  scale_fill_manual(values=c(`0`="#FED068",`1`="#8540A0")) +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons , label = "p.signif") +
  ylab("Shannon Index") +
  theme(axis.title.y = element_text(face = "bold", size=14),
        axis.title.x = element_text(face = "bold", size=14),
        axis.text.x = element_text(face = "bold", size=14),
        strip.background = element_blank())

gridExtra::grid.arrange(a,b, ncol=2)

###################################
# identify top 10 abundant genera #
###################################

# 1. Agglomerate taxa to the Genus level
# This sums up the counts for all ASVs/OTUs belonging to the same genus
physeq_genus <- tax_glom(physeq.prune, taxrank = "Genus")

# 2. Calculate the total abundance for each genus across all samples
# taxa_sums() calculates the total abundance per taxon
genus_abundances <- taxa_sums(physeq_genus)

# 3. Sort the abundances in decreasing order
sorted_abundances <- sort(genus_abundances, decreasing = TRUE)

# 4. Select the top 10 most abundant genera
top_10_genera_names <- names(head(sorted_abundances, 10))

# 5. Extract the taxonomy table and filter for the top 10 names
top_10_tax_table <- tax_table(physeq_genus)[top_10_genera_names, ]

# 6. Convert the taxonomy table to a data frame for better viewing
top_10_df <- as.data.frame(top_10_tax_table)

# 7. Add the abundance values to the data frame
top_10_df$TotalAbundance <- sorted_abundances[top_10_genera_names]

# 8. Display the top 10 abundant genera
print(top_10_df)

# 9. Calculate average by group
cntGenera <-as.data.frame(otu_table(physeq_genus))

cntGeneraSub <- cntGenera[top_10_genera_names,]

te<-cntGenera[top_10_genera_names,] %>% filter(rownames(.) %in% top_10_genera_names) %>% rowwise() %>%
  mutate(avg_0=mean(c_across(MDBC01:MDBC30),na.rm = TRUE),
         avg_1=mean(c_across(MDBC32:MDBC46),na.rm = TRUE))

######################################################
# normalize and transform data for relative analysis #
######################################################

# normalize reads using median sequencing depth
total <- median(sample_sums(physeq.prune)) #36456
standf <- function(x, t=total) round(t * (x / sum(x)))
physeq.norm <- transform_sample_counts(physeq.prune, standf)

# transform to relative abundance
physeq2 <- transform_sample_counts(physeq.norm, function(x) x / sum(x) )
physeq2

ps2 <- psmelt(physeq2) %>% filter (year !="N")
ps2$Genus <- as.character(ps2$Genus) #convert to character
ps2$Order <- as.character(ps2$Order) #convert to character

# bar plot by genus
ps2$Genus<-ifelse(ps2$Genus %in% c("BD1-7 clade", "Mycoplasma", "Vibrio", "Mycobacterium", "Methylobacterium", "FS140-16B-02 marine group",
                                   "Neptuniibacter","Pseudoalteromonas","Endozoicomonas", "Marivibrio"), ps2$Genus, "other")
ps2$Genus<-factor(ps2$Genus, levels=c("BD1-7 clade", "Mycoplasma", "Mycobacterium", "Pseudoalteromonas",
                                      "Neptuniibacter","Vibrio", "FS140-16B-02 marine group","Methylobacterium","Endozoicomonas", "Marivibrio", "other"))


p <- ggplot(data=ps2, aes(x=year, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c(`BD1-7 clade`="#14655E", `Mycoplasma`="#B9DFBC",Mycobacterium="#1A936F", 
                               Pseudoalteromonas="#49AB7C",
                               Neptuniibacter="#3BC95A", Vibrio="#14847A","FS140-16B-02 marine group"="#BFEDC9", 
                               Methylobacterium= "#DDE2C9", Endozoicomonas="#F3E9D2",  
                               Marivibrio="#F7F3EC", `other`="lightgrey")) +
  theme_classic() + ylab( "Relative Abundance %") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +facet_grid(~Longitude*Alt_ID, scales="free") +
  scale_y_continuous(expand = c(0,0))


###################
# beta-diversity  #
###################

# create an ordination plot
ord <- ordinate(physeq.norm,"NMDS","bray")

#stressplot(ord)
plot_ordination(physeq.norm,ord,type="samples",color="year", shape="location") + geom_point(size=5) +
  stat_ellipse( linetype = 2,alpha=0.2, geom = "polygon",aes(fill=year, group=year)) +
  scale_color_manual(values=c(`0`="#FED068",`1`="#8540A0","grey")) +
  scale_fill_manual(values=c(`0`="#FED068",`1`="#8540A0","grey")) +
  theme_bw() +
  geom_text_repel(mapping = aes(label = Alt_ID), size = 4)


## PERMANOVA
library(pairwiseAdonis)

pdist<- phyloseq::distance(physeq.norm, method="bray")

perm <- adonis2(pdist ~ year, data = metadata2, permutations = 9999)
perm
pairwise.adonis(pdist, metadata2$year)

# how much variation in data
dispersion_result <- betadisper(pdist, metadata2$year)

plot(dispersion_result)
boxplot(dispersion_result)

########################################
## Differential abundance with DEseq2 ##
########################################

ps.taxa <- tax_glom(physeq.prune, taxrank="Genus")
ps.taxa

# convert to deseq object
ps_ds <- phyloseq_to_deseq2(ps.taxa, ~ year)

# use alternative estimator on a condition of "every gene contains a sample with a zero"
ds <- estimateSizeFactors(ps_ds, type="poscounts")

# run differential abundance test
ds = DESeq(ds, test="Wald", fitType="parametric")
alpha = 0.05
res <- results(ds, alpha=alpha)
res <- as.data.frame(res)
res.filt <- res %>% filter(padj < 0.05) 

# subset to most abundant ASVs
taxa_sig <- rownames(res.filt) # significant taxa plus two dominant asv
ps.taxa.sig <- prune_taxa(top_10_genera_names, ps.taxa)
taxtab <- as.data.frame(tax_table(ps.taxa.sig)) %>% rownames_to_column("ASV")

res.filt2 <- res %>% rownames_to_column("ASV") %>% filter(ASV %in% top_10_genera_names) %>% arrange(desc(log2FoldChange)) %>% 
  left_join(taxtab, by=c("ASV" = "ASV")) %>% left_join(top_10_df %>% dplyr::select(Genus, TotalAbundance), by=c("Genus"))

# plot differential abundance results
ggplot(res.filt2, aes(x=log2FoldChange, y=reorder(Genus, log2FoldChange), fill=log2FoldChange, size=log(TotalAbundance,10))) +
  geom_errorbarh(mapping=aes(y=Genus, x=log2FoldChange+lfcSE, xmin=log2FoldChange+lfcSE, xmax=log2FoldChange-lfcSE), 
                 height=0, linewidth=1, color="darkgrey") +
  geom_point(shape=21) + 
  geom_vline(xintercept=0, linetype=2) +
  theme_bw() +
  scale_size_continuous(limits = c(2,6)) +
  scale_fill_gradient2(low = "#FED068", mid = "white", high = "#8540A0")+
  ylab("") 

# Get the data for a specific gene, returning data
ASV= "ASV_49"
plot_data <- plotCounts(ds, gene=ASV, intgroup="year", returnData=TRUE)

# Plot with ggplot2
ggplot(plot_data, aes(x=year, y=count, fill=year)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(w=0.1, h=0), shape=21) + # Jitter points for visibility
  #geom_text_repel(aes(label=rownames(plot_data))) + # Add sample labels
  theme_bw() +
  ggtitle(ASV) +
  theme(plot.title = element_text(hjust=0.5))+
  scale_fill_manual(values=c(`0`="#FED068",`1`="#8540A0"))


#########################
an <-ANCOMBC::ancombc2(data = ps.taxa, 
              tax_level = NULL,
              fix_formula = "year",
              rand_formula = "(1 | Alt_ID)",
              p_adj_method = "fdr",
              prv_cut = 0.2,
              group="year",
              struc_zero = TRUE,
              neg_lb = TRUE,
              global = TRUE)

#rand_formula = "(1 | Alt_ID)",

an$res

results_ancom_bc = data.frame(ASV = an$res$taxon,
                              lcf = an$res$lfc_year1, 
                              W = an$res$W_year1, 
                              p_val = an$res$p_year1, 
                              q_value = an$res$q_year1, 
                              Diff_ab =  an$res$diff_year1)

results_ancom_bc$lcf = results_ancom_bc$lcf * -1

############################# 
## Absolute quantification  #
#############################

absDat <- read.table("I:/My Drive/Carballosa_Swiftia_Msphere_2026/Code_Data/absolute.abundance.csv",sep = ",", header = T) %>%
  left_join(metadata, by=c("customer_label"="Sample_ID")) %>% filter( year != "N")

absDat %>% 
  summarize(med=median(genome_copies_per_ul.),sd=sd(genome_copies_per_ul.) )

ggplot(absDat, aes(x=year, y=`DNA_ng_per_ul.`, fill=year)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(w=0.1, h=0), shape=21) + # Jitter points for visibility
  #geom_text_repel(aes(label=rownames(plot_data))) + # Add sample labels
  theme_bw() +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons , label = "p.signif") +
  theme(plot.title = element_text(hjust=0.5))+
  scale_fill_manual(values=c(`0`="#FED068",`1`="#8540A0")) + ylim (0,0.2)

compare_means(`DNA_ng_per_ul.` ~ year, absDat,paired = TRUE)
compare_means(gene_copies_per_ul ~ year, absDat, paired = TRUE, method = "t.test")

# Subset Bacteria data before treatment
wild <- subset(absDat,  year == "0", DNA_ng_per_ul.,
               drop = TRUE)
# subset Bacteria data after treatment
captive <- subset(absDat ,  year == "1",  DNA_ng_per_ul.,
                  drop = TRUE)

pd <- paired(wild, captive)

p<-ggpaired(pd, cond1 = "wild", cond2 = "captive",
            fill = "condition", palette = c("#FED068","#8540A0"),
            ylab = "% Bacterial rRNA", xlab = "", ylim=c(0, 0.2),
            font.label = list(size = 21, color = "black"))

p
p + stat_compare_means( paired = TRUE)
