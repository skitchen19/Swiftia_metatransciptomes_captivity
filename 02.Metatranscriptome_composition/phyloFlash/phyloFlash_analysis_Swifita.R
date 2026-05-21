#####################
# phyloFlash output #
#####################

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

####################################
# create phyloseq object: Bacteria #
####################################

## Bacteria
# import your files
table <- read.table(file = "bacteria.ALL_count_table.txt", sep = "\t", header = T) %>%
  mutate(OTU = paste0("NTU",OTU))
table.noMito <- table %>% column_to_rownames(var = "OTU") %>% filter(Family != "Mitochondria")
table.noApi <- table.noMito %>% filter(Order != "Chloroplast")
table.polish <- table.noApi %>% mutate(across(everything(), ~replace(., str_detect(., "\\(.*\\)"), NA))) %>%
  mutate(across(everything(), ~replace(., str_detect(., "uncultured\\.*"), NA))) %>%
  mutate(across(everything(), ~replace(., str_detect(., "metagenome"), NA)))

# split for OTUs
otu <-table.polish[,c(7:29)]
OTU <- otu_table(otu,taxa_are_rows=TRUE)

# split for taxa of OTUs
tax <- table.polish[,c(1:6)]
TAX <- tax_table(as.matrix(tax))

# load in metadata
metadata <- read.table("../metadata_swiftia.txt",sep = "\t", header = T)
#metadata$Sample_ID<-str_replace(metadata$Sample_ID, pattern="_rep", replacement="")
metadata2 <-metadata %>% dplyr::filter(Alt_ID !="R") %>% column_to_rownames("Sample_ID") %>%
  filter(year != "N")
metadata2$year<-as.factor(metadata2$year)
metadata2$Alt_ID<-as.factor(metadata2$Alt_ID)
metadata2$location<-as.factor(metadata2$location)

# structure of data
str(metadata2)

# sample table for phyloseq
SAMPLE <- sample_data(metadata2)

# create phyloseq object
physeq <- phyloseq(OTU, TAX, SAMPLE)
physeq

########################################
# remove ASV with only one observation #
########################################

physeq.prune <- prune_taxa(taxa_sums(physeq) > 1, physeq)
physeq.prune

#############################
# alpha diversity estimates #
#############################

rich <- estimate_richness(physeq.prune, split=TRUE, measures=c("Chao1","Shannon")) %>% 
  mutate(year=metadata2$year) 

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


######################################################
# normalize and transform data for relative analysis #
######################################################

total <- median(sample_sums(physeq.prune)) #253477
standf <- function(x, t=total) round(t * (x / sum(x)))
physeq.norm <- transform_sample_counts(physeq.prune, standf)

# transform to relative abundance
physeq2 <- transform_sample_counts(physeq.norm, function(x) x / sum(x) )
physeq2

ps2 <- psmelt(physeq2)
ps2$Genus <- as.character(ps2$Genus) #convert to character
ps2$Order <- as.character(ps2$Order) #convert to character

# Genus
ps2$Genus<-ifelse(ps2$Genus %in% c("BD1-7 clade", "Mycoplasma", "Vibrio", "Mycobacterium", "Methylobacterium-Methylorubrum", "FS140-16B-02 marine group",
                                   "Neptuniibacter","Pseudoalteromonas","Endozoicomonas", "Marivibrio"), ps2$Genus, "other")
ps2$Genus<-factor(ps2$Genus, levels=c("BD1-7 clade", "Mycoplasma", "Mycobacterium", "Pseudoalteromonas",
                                      "Neptuniibacter","Vibrio", "FS140-16B-02 marine group","Methylobacterium-Methylorubrum","Endozoicomonas", "Marivibrio", "other"))

p <- ggplot(data=ps2, aes(x=year, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c(`BD1-7 clade`="#14655E", `Mycoplasma`="#A7D7AC",Mycobacterium="#1A936F", Pseudoalteromonas="#51B484",
                               Neptuniibacter="#88D498", Vibrio="#166F67","FS140-16B-02 marine group"="#C6DABF", 
                               `Methylobacterium-Methylorubrum`= "#DDE2C9", Endozoicomonas="#F3E9D2",  
                               Marivibrio="#F4EBD6", `other`="lightgrey")) +
  theme_classic() + ylab( "Relative Abundance %") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +facet_grid(~Longitude*Alt_ID, scales="free") +
  scale_y_continuous(expand = c(0,0))


###################
# beta-diversity  #
###################

# create an ordination plot
ord <- ordinate(physeq.norm,"NMDS","bray")
plot_ordination(physeq,ord,type="samples",color="year", shape="location") + geom_point(size=5) +
  stat_ellipse( linetype = 2,alpha=0.2, geom = "polygon",aes(fill=year, group=year)) +
  scale_color_manual(values=c(`0`="#FED068",`1`="#8540A0")) +
  scale_fill_manual(values=c(`0`="#FED068",`1`="#8540A0")) +
  theme_bw() +
  geom_text_repel(mapping = aes(label = Alt_ID), size = 4)

## PERMANOVA
pdist<- phyloseq::distance(physeq.norm, method="bray")

perm <- adonis2(pdist ~ year, data = metadata2, permutations = 999)
perm
pairwise.adonis(pdist, metadata2$year)

# how much variation in data
dispersion_result <- betadisper(pdist, metadata2$year)

plot(dispersion_result)
boxplot(dispersion_result)

# anova
anova(dispersion_result)

# Tukey's Honest Significant Differences
T_HSD <- TukeyHSD(dispersion_result)
T_HSD

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
#res.filt <- head(res, n=15)

#taxa_sig <- rownames(res.filt) # significant taxa plus two dominant asv
top_10_ntus<-c("NTU2379","NTU999","NTU93","NTU2138","NTU2335", "NTU2161","NTU1429", "NTU2279","NTU1613")

ps.taxa.sig <- prune_taxa(top_10_ntus, ps.taxa)
taxtab <- as.data.frame(tax_table(ps.taxa.sig)) %>% rownames_to_column("ASV")

# read counts per genus
genus_abundancesMT <- data.frame(`TotalAbundance`= taxa_sums(ps.taxa)) %>% rownames_to_column(var = "ASV")

res.filt2 <- res %>% rownames_to_column("ASV") %>% filter(ASV %in% top_10_ntus) %>% arrange(desc(log2FoldChange)) %>% 
  left_join(taxtab, by=c("ASV" = "ASV")) %>% left_join(genus_abundancesMT, by=c("ASV"))


# Plot with ggplot2
ggplot(res.filt2, aes(x=log2FoldChange, y=reorder(Genus, log2FoldChange), fill=log2FoldChange, size=log(TotalAbundance,10))) +
  geom_errorbarh(mapping=aes(y=Genus, x=log2FoldChange+lfcSE, xmin=log2FoldChange+lfcSE, xmax=log2FoldChange-lfcSE), 
                 height=0, linewidth=1, color="grey") +
  geom_point(shape=21) + 
  geom_vline(xintercept=0, linetype=2) +
  theme_bw() +
  scale_fill_gradient2(low = "#FED068", mid = "white", high = "#8540A0") +  ylab("") 


# Get the data for a specific gene, returning data
NTU= "NTU93"
plot_data <- plotCounts(ds, gene=NTU, intgroup="year", returnData=TRUE)

# Plot with ggplot2
ggplot(plot_data, aes(x=year, y=count, color=year)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(w=0.1, h=0)) + # Jitter points for visibility
  #geom_text_repel(aes(label=rownames(plot_data))) + # Add sample labels
  theme_bw() +
  ggtitle(res.filt2$Genus[res.filt2$ASV ==NTU]) +
  theme(plot.title = element_text(hjust=0.5))+
  scale_color_manual(values=c(`0`="#FED068",`1`="#8540A0"))

#########################
## Apicomplexa 16S Data #
#########################

# Keep only apicomplexan ASVs
taxTAB <-as.data.frame(phyloseq_object@tax_table) %>% 
  subset(.,(rownames(.) %in% apicomplexan_asvs)) 

# make taxonomy table
TAX = tax_table(as.matrix(taxTAB))

# remove contaminants from ASV table
otuTAB<-as.data.frame(phyloseq_object@otu_table)
colnames(otuTAB) <- metadata$Sample_ID

otuTAB <-  otuTAB[rownames(otuTAB) %in% rownames(taxTAB), ] 

# make otu table
OTU <- otu_table(otuTAB[,1:24],taxa_are_rows=TRUE)

# create phyloseq object
physeqApi <- phyloseq(OTU, TAX, SAMPLE)
physeqApi

ps.api <- tax_glom(physeqApi, taxrank="Order")

as<-as.data.frame(t(sample_sums(ps.api))) %>% rowwise() %>%
  mutate(avg_0=mean(c_across(MDBC01:MDBC30),na.rm = TRUE),
         avg_1=mean(c_across(MDBC32:MDBC46),na.rm = TRUE))

# plot and test for difference in 16S
psApi <- psmelt(ps.api)
a_my_comparisons <- list( c("0", "1"))
ggplot(psApi, aes(year, (Abundance), fill = year)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitter(w=0.1, h=0), shape=21) + # Jitter points for visibility
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif") +
  ylab("Plastid Read Count") +
  theme_classic() +scale_fill_manual(values=c(`0`="#FED068",`1`="#8540A0")) 

compare_means(Abundance ~ year, method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif",data = psApi)

p <- ggplot(data=psApi, aes(x=year, y=Abundance, fill=OTU)) +
  geom_bar(aes(), stat="identity", position="stack", fill="#041D25") +
  theme_classic() + ylab( "Plastid Counts") +
  theme(legend.position="none") + guides(fill=guide_legend(nrow=2)) +facet_grid(~Longitude*Alt_ID, scales="free") +
  scale_y_continuous(expand = c(0,0), limits=c(0,3500))

p


########################
## Eukaryotic 18S data #
########################

# load in table from phyloflash
table <- read.table(file = "eukaryota.ALL_table_L7.txt", sep = "\t", header = T) %>%
  mutate(OTU = paste0("NTU",OTU)) %>% column_to_rownames(var = "OTU")
table.polish <- table %>% 
  mutate(across(everything(), ~replace(., str_detect(., "uncultured\\.*"), NA))) %>%
  mutate(across(everything(), ~replace(., str_detect(., "metagenome"), NA))) %>%
  rename_with(~str_remove(., "_rep"))

# split for OTUs
otu <-table.polish[,c(9:31)]
OTU <- otu_table(otu,taxa_are_rows=TRUE)


# split for taxa of OTUs
tax <- table.polish[,c(1:8)]
TAX <- tax_table(as.matrix(tax))

# sample table for phyloseq
SAMPLE <- sample_data(metadata2)

# create phyloseq object
physeq <- phyloseq(OTU, TAX, SAMPLE)
physeq

########################################
# remove ASV with only one observation #
########################################

physeq.prune <- prune_taxa(taxa_sums(physeq) > 1, physeq)
physeq.prune

#collapse to Class level
ps.class <- tax_glom(physeq.prune, taxrank="Class")
ps.class

# normalize reads using median sequencing depth
total <- median(sample_sums(ps.class)) #56630.5
standf <- function(x, t=total) round(t * (x / sum(x)))
physeq.norm <- transform_sample_counts(ps.class, standf)

# transform to relative abundance
physeq2 <- transform_sample_counts(physeq.norm, function(x) x / sum(x) )
physeq2

ps2 <- psmelt(ps.class)
ps2$Class <- as.character(ps2$Class) #convert to character

ps2$Class[ps2$Abundance < 0.05] <- "<5% abund."

ps2$Class<-factor(ps2$Class, levels=c("Octocorallia", "Conoidasida","Dinophyceae", "Spirotrichea",
                                      "Globothalamea","(Dikarya)","<5% abund."))

p <- ggplot(data=ps2, aes(x=year, y=Abundance, fill=Class))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c(`Octocorallia`="#1389b1", `Conoidasida`="#041D25", 
                               `Dinophyceae`="#9EC8D3",`Spirotrichea`="#EAF5F7", `Globothalamea`="#0E6481", `(Dikarya)`="#2BA8CF", `<5% abund.`= "lightgrey")) +
  theme_classic() + ylab( "Relative Abundance %") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +facet_grid(~Longitude*Alt_ID, scales="free") +
  scale_y_continuous(expand = c(0,0))

p <- ggplot(data=ps2, aes(x=year, y=Abundance, fill=Class))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c(`Opisthokonta`="grey30", `Apicomplexa`="lightblue", 
                               `Dinoflagellata`="grey95",`Ciliophora`="black", `Retaria`="grey10", `<5% abund.`= "lightgrey")) +
  theme_classic() + ylab( "Relative Abundance %") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  facet_grid(~location*Alt_ID, scales="free") +
  scale_y_continuous(expand = c(0,0))

ps265 <- ps2 %>% filter(OTU == "NTU265")

compare_means(Abundance ~ year, method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif",data = ps265)
