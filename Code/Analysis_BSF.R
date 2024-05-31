####################################
# SCRIPT 1: CREATE A MASTERTABLE
####################################
####################################
# 0. Packages: install packages and load them wherever needed
####################################
#install.packages("FSA")
library(lme4)
library (microbiome)
library(FSA)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#install.packages("plyr")
library(plyr)
#install.packages("reshape")
library(reshape)
#install.packages("reshape2")
library(reshape2)
#install.packages("scales")
library(scales)
#install.packages("viridis")
library(viridis)
#install.packages("readr")
library(readr)
#install.packages("WriteXLS")
library(WriteXLS)
#install.packages("readxl")
library(readxl)
#install.packages("remotes")
library(remotes)
#install.packages("devtools") # I have the problem that Mac Big sure blocks access to all kind of system on the mac and had to disable system integratie protection
library(devtools)
#install.packages("BiocManager") # # before downlading from github you need to load devtools
library(BiocManager)
#BiocManager::install("ComplexHeatmap")
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
#source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local = TRUE) ## like this it works!!
library(phyloseq)
#BiocManager::install("ALDEx2")
library(ALDEx2)
#install.packages("RSQLite")
library(RSQLite)
#install.packages("https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.9.900.3.0.tar.gz", repos=NULL, type="source")
library(RcppArmadillo)
#BiocManager::install("DESeq2")
library(DESeq2)
#remotes::install_github("yiluheihei/microbiomeMarker")
library(microbiomeMarker)
#install.packages("vegan")
library(vegan)
#install.packages("dendsort")
library(dendsort)
#install.packages("ape")
library(ape)
#install.packages("seriation")
library(seriation)
#install.packages("taxonomizr")
library(taxonomizr)
#install.packages("esquisse")
library(esquisse)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("grid")
library(grid)
#install.packages("coin")
library(coin)
#install.packages("rstatix")
library(rstatix)
#install.packages("PMCMR")
library(PMCMR)
#remotes::install_github("tapj/biotyper")
library(biotyper)
#BiocManager::install("DirichletMultinomial")
#install.packages("DirichletMultinomial")
library(DirichletMultinomial)
#install.packages("shiny")
library(shiny)
#install.packages("fpc")
library (fpc)
#install.packages("clValid")
library (clValid)
#install.packages("cluster")
library(cluster) 
#install.packages("philentropy")
library(philentropy)
####################################
# 1. Create abundance table 
####################################
## A) set directory
setwd("/Users/daan/Desktop/Transfer/Abundances")
getwd()
dir()

## B) read the abundance text files and store in a list
list <- list()
list_txt <- dir(pattern = "*.magnitudes.txt", full.names = FALSE)
for (k in 1: length(list_txt)){
  list[[k]] <- read.delim(list_txt[k])
}

## C) create names file
names <- sub(pattern = "*.magnitudes.txt", "", list_txt)

## D) evaluate lists
head(list[[1]]) ## perform this command to check if the first header is an 'X', needed to merge in next command

## E) merge lists in a table by matching headers
table <- merge(list[[1]], list[[2]], by="magnitudes", all = TRUE) 
for (i in 1:((length(list))-2)){
  table <- merge(table, list [[i+2]], by = "magnitudes", all = TRUE)
}

## F) create a dataframe, discard first row and make contigs the row names,
table <- data.frame(table[,-1], row.names = table[,1])

## G) change all colnames by names vector
for (j in 1:length(names)){
  colnames(table) <- c(names)
}

## H) convert NA to 0 and dots with underscores
table[is.na(table)] <- 0
rownames(table) <- gsub("\\.{1}", "_", rownames(table))

## I) Remove sample with no reads mapped - 51 sampels present
table <- table[, colSums(table != 0) > 0]
ncol(table)
####################################
# 2. Create vector of contigs
####################################
setwd("/Users/daan/Desktop/Transfer/")
getwd()
dir()

non_redundant_contig_names <- read.table("names_contigs.txt", quote = "\"", comment.char = "")
colnames(non_redundant_contig_names) <- "contigs_orginal_name.txt"
rownames(non_redundant_contig_names) <- non_redundant_contig_names$contigs_orginal_name.txt
rownames(non_redundant_contig_names) <- gsub("\\.{1}", "_", rownames(non_redundant_contig_names))
contigs_project <- non_redundant_contig_names$contigs_orginal_name.txt

#view(non_redundant_contig_names) # 'view' is performed with tibble
#View(contigs_project)           # 'View' is performed with utils
####################################
# 3. Viral identification and classification
####################################
# 3.1 Genomead
####################################
setwd("/Users/daan/Desktop/Transfer/Taxonomy")
getwd()
dir()

AllSamples_classificationtaxonomy_genomead <- as.data.frame(read_delim("Genomead3_final.txt", "\t", escape_double = FALSE, col_names = TRUE, trim_ws = FALSE))
AllSamples_classificationtaxonomy_genomead[is.na(AllSamples_classificationtaxonomy_genomead)] <- "Unclassified"
AllSamples_classificationtaxonomy_genomead <- AllSamples_classificationtaxonomy_genomead[AllSamples_classificationtaxonomy_genomead$NODE %in% contigs_project,]
rownames(AllSamples_classificationtaxonomy_genomead) <- AllSamples_classificationtaxonomy_genomead$NODE
AllSamples_classificationtaxonomy_genomead$NODE <- NULL
rownames(AllSamples_classificationtaxonomy_genomead) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_genomead))

####################################
# 3.2 Insert Completeness
####################################
AllSamples_completeness_Checkv <- as.data.frame(read_delim("quality_summary_R.tsv", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_completeness_Checkv) <- c('contig', 'total genes', 'viral genes', 'completeness_quality', 'completeness')
rownames(AllSamples_completeness_Checkv) <- AllSamples_completeness_Checkv$contig
AllSamples_completeness_Checkv$contig <- NULL
AllSamples_completeness_Checkv$`total genes` <- NULL
AllSamples_completeness_Checkv$`viral genes` <- NULL
#AllSamples_completeness_Checkv$`total genes` <- NULL
#AllSamples_completeness_Checkv$`viral genes` <- NULL
rownames(AllSamples_completeness_Checkv) <- gsub("\\.{1}", "_", rownames(AllSamples_completeness_Checkv))
AllSamples_completeness_Checkv$completeness[is.na(AllSamples_completeness_Checkv$completeness)] <- 0
AllSamples_completeness_Checkv <- AllSamples_completeness_Checkv[!AllSamples_completeness_Checkv$completeness==0,]
# Check Blastn.txt file for genus/species taxonomy.
####################################
# 4. In silico host prediction
####################################
# 4.1 In silico host Phyla
####################################
setwd("/Users/daan/Desktop/Transfer/Host")
getwd() 
dir()

AllSamples_completeness_Host_phyla <- as.data.frame(read_delim("Host_phylum.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_completeness_Host_phyla) <- c('contig', 'Host_phyla')
rownames(AllSamples_completeness_Host_phyla) <- AllSamples_completeness_Host_phyla$contig
AllSamples_completeness_Host_phyla$contig <- NULL
rownames(AllSamples_completeness_Host_phyla) <- gsub("\\.{1}", "_", rownames(AllSamples_completeness_Host_phyla))
#View(AllSamples_completeness_Host_phyla)
####################################
# 4.2 In silico host Genus
####################################
getwd() 
dir()

AllSamples_completeness_Host_genus <- as.data.frame(read_delim("Host_genus.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_completeness_Host_genus) <- c('contig', 'Host_genus')
rownames(AllSamples_completeness_Host_genus) <- AllSamples_completeness_Host_genus$contig
AllSamples_completeness_Host_genus$contig <- NULL
rownames(AllSamples_completeness_Host_genus) <- gsub("\\.{1}", "_", rownames(AllSamples_completeness_Host_genus))
#View(AllSamples_completeness_Host_genus)
####################################
# 5. Prophage prediction
####################################
setwd("/Users/daan/Desktop/Transfer/Function")
getwd() 
dir()

AllSamples_lyso <- as.data.frame(read_delim("All_predicted_protein_lysogenic_nodes.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllSamples_lyso) <- c('contig', 'temperate')
rownames(AllSamples_lyso) <- AllSamples_lyso$contig
AllSamples_lyso$contig <- NULL
rownames(AllSamples_lyso) <- gsub("\\.{1}", "_", rownames(AllSamples_lyso))
#View(AllSamples_lyso)
####################################
# 6. Add genus level
####################################
setwd("/Users/daan/Desktop/Transfer/Taxonomy")
getwd() 
dir()

AllsamplesGenus <- as.data.frame(read_delim("genus_clusters2.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllsamplesGenus) <- c('contig', 'Genus')
rownames(AllsamplesGenus) <- AllsamplesGenus$contig
AllsamplesGenus$contig <- NULL
rownames(AllsamplesGenus) <- gsub("\\.{1}", "_", rownames(AllsamplesGenus))
####################################
# 7. Create mastertable: merge scaffolds, abundances and annotations
####################################
COMBO_1 <- merge(table, AllSamples_classificationtaxonomy_genomead, by = 0, all = F)
rownames(COMBO_1) <- COMBO_1$Row.names
COMBO_1$Row.names <- NULL
columns_to_remove <- grep("^NC", names(COMBO_1))
COMBO_1 <- COMBO_1[, -columns_to_remove]
COMBO_1 <- COMBO_1 %>% select(-starts_with("NC"))

COMBO_2 <- merge(COMBO_1, AllSamples_completeness_Checkv, by = 0, all = F)
rownames(COMBO_2) <- COMBO_2$Row.names
COMBO_2$Row.names <- NULL
dim(COMBO_2)

COMBO_3 <- merge(COMBO_2, AllSamples_completeness_Host_phyla, by=0, all=T)
rownames(COMBO_3) <- COMBO_3$Row.names
COMBO_3$Row.names <- NULL
dim(COMBO_3)

COMBO_4 <- merge(COMBO_3,AllSamples_completeness_Host_genus, by=0, all=T)
rownames(COMBO_4) <- COMBO_4$Row.names
COMBO_4$Row.names <- NULL
dim(COMBO_4)

COMBO_5 <- merge(COMBO_4,AllsamplesGenus, by=0, all=T)
rownames(COMBO_5) <- COMBO_5$Row.names
COMBO_5$Row.names <- NULL
dim(COMBO_5)

Mastertable <- merge(COMBO_5,AllSamples_lyso, by=0, all=T)
rownames(Mastertable) <- Mastertable$Row.names
Mastertable$Row.names <- NULL

Mastertable$temperate[is.na(Mastertable$temperate)] <- 'Virulent'
Mastertable$temperate[!Mastertable$temperate == "Virulent"] <- 'Temperate'
####################################
# 6.1 Additional information to Mastertable
####################################
# 6.1.1 Insert total number of trimmed reads
####################################
Mastertable[is.na(Mastertable)] <- 0
Mastertable$Totalnumberofreads <- as.numeric(rowSums(Mastertable[,colnames(Mastertable) %in% names]))
Mastertable <- Mastertable[!(Mastertable$Totalnumberofreads == "0"),] ## remove rows with total number of zero
####################################
# 6.1.2 Insert length of contigs
####################################
Mastertable$length <- rownames(Mastertable)
Mastertable$length <- gsub("*_cov.*","",Mastertable$length)
Mastertable$length <- gsub(".*length_","", Mastertable$length)
Mastertable$length <- as.numeric(Mastertable$length)
Mastertable$length <- Mastertable$length/1000
ncol(Mastertable)
#View(Mastertable[c(1,7:9,47:56)])
####################################
# 7 viral identification and completeness correction
####################################
Mastertable_viral <- Mastertable[Mastertable$Virus=="phage" | Mastertable$Virus=='eukaryotic_virus',]
Mastertable_HQ_viral <- Mastertable_viral[Mastertable_viral$completeness > 50,]
#View(Mastertable_HQ_viral)
####################################
# 8. Phyloseq
####################################
# 8.1 Metadata
####################################
library(readxl)
setwd("/Users/daan/Desktop/Transfer/metadata")
metadata <- as.data.frame(read_excel("metadata.xlsx", sheet = "metadata"))
rownames(metadata) <- metadata$Sample_ID
metadata[metadata == "NA"] <- NA
#View(metadata)
####################################
# 8.2 Abundance table
####################################
sort(colSums(Mastertable_HQ_viral[1:45]))
Mastertable_HQ_viral$Hae_CF <- NULL
Mastertable_HQ_viral$WL_SW_Small <- NULL
Mastertable_HQ_viral$PMG_SW_6 <- NULL

names <- colnames(Mastertable_HQ_viral[1:42][,colSums(Mastertable_HQ_viral[1:42])!=0])

abundance_table <- Mastertable_HQ_viral
abundance_table$names <- rownames(abundance_table)
abundance_tabe1 <- aggregate(. ~names, FUN = sum, data = abundance_table[,colnames(abundance_table) %in% names | colnames(abundance_table) == 'names'])
rownames(abundance_tabe1) <- abundance_tabe1$names
abundance_tabe1$names <- NULL
####################################
# 8.3 Taxonomy table
####################################
colnames(Mastertable_HQ_viral)
Mastertable_viral_rarefied_2 <- Mastertable_HQ_viral
Mastertable_viral_rarefied_2$Totalnumberofreads <- NULL
Mastertable_viral_rarefied_2$length <- NULL
Mastertable_viral_rarefied_2$temperate <- NULL
Mastertable_viral_rarefied_2$Host_genus <- NULL
Mastertable_viral_rarefied_2$Host_phyla <- NULL
Mastertable_viral_rarefied_2$completeness <- NULL
Mastertable_viral_rarefied_2$completeness_quality <- NULL
colnames(Mastertable_viral_rarefied_2)

vector_1 <- (which(names(Mastertable_viral_rarefied_2)== "Virus"))
vector_2 <- which(names(Mastertable_viral_rarefied_2)== "Genus")
Taxonomy_table <- Mastertable_viral_rarefied_2[,c(vector_1:vector_2)]
colnames(Taxonomy_table)

## Complete taxonomies
names(Taxonomy_table)[names(Taxonomy_table) == "Virus"] <- "Kingdom"
names(Taxonomy_table)[names(Taxonomy_table) == "Order "] <- "Order"
Taxonomy_table$Phylum <- Taxonomy_table$Kingdom
Taxonomy_table$Species <- rownames(Taxonomy_table)
colnames(Taxonomy_table)

# Define the desired column order
desired_order <- c("Genus","Kingdom", "Phylum", "Class", "Order", "Family", "Species")

# Reorder the columns of Taxonomy_table
Taxonomy_table <- Taxonomy_table[, desired_order]
colnames(Taxonomy_table) # Only use Kingdom, Class, order and Family
####################################
# 8.4 Create phyloseq table
####################################
abundance_table_m <- as.matrix(abundance_tabe1)
taxonomy_table_m <- as.matrix(Taxonomy_table)

ABUNDANCE_rarefied_family <- otu_table(abundance_table_m, taxa_are_rows = T)
TAX_rarefied_family <- tax_table(taxonomy_table_m)
samples <- sample_data(metadata)

phyloseq_samples <- phyloseq(ABUNDANCE_rarefied_family,TAX_rarefied_family, samples)
phyloseq_samples
####################################
# 8.5 Visualize virome
####################################
# 8.5.1 Complete dataset
####################################
phyloseq_samples_virome <- phyloseq_samples

phyloseq_virome_transformed <- phyloseq_samples_virome %>% 
  microbiome::transform(transform = "log") %>%
  aggregate_rare(level = "Species", detection = 0.0005/100, prevalence = 0.00001/100)

# Distances
phyloseq_virome_transformed_bray <- phyloseq::distance(phyloseq_virome_transformed, method="bray")
phyloseq_virome_transformed_pcOA <- ordinate(phyloseq_virome_transformed, method="PCoA", distance=phyloseq_virome_transformed_bray)

# PcoA 1: Substrate_type shaped by Sample.type
plot_ordination(phyloseq_virome_transformed, phyloseq_virome_transformed_pcOA, color = "Substrate_type") +
  geom_point(alpha = 1, size = 4.5) +
  theme_bw() +
  scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (24.8%)") +
  ylab("PC1 (25.8%)") +
  labs(color = "Substrate") +
  theme(axis.text.x = element_text(size = 11, colour = "black"),  # Set x-axis text size
      axis.text.y = element_text(size = 11, colour = "black"),  # Set y-axis text size
      panel.border = element_rect(colour = "black", size = 0.7))  # Set border color and thickness

# PcoA 2:
fillable_shapes <- c(21, 22, 23, 24, 25)  # Reusing shapes as needed

plot_ordination(phyloseq_virome_transformed, phyloseq_virome_transformed_pcOA, shape = "Substrate_type") +
  geom_point(aes(fill = Substrate_type), size = 4.5, color = "black") +  # Use fill aesthetic and set point outline color
  theme_bw() +
  scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +  # Define fill colors
  scale_shape_manual(values = fillable_shapes) +  # Use shapes that support filling
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (24.8%)") +
  ylab("PC1 (25.8%)") +
  labs(color = "Substrate") +
  theme(axis.text.x = element_text(size = 11, colour = "black"),  # Set x-axis text size
        axis.text.y = element_text(size = 11, colour = "black"),  # Set y-axis text size
        panel.border = element_rect(colour = "black", size = 0.7))  # Set border color and thicknes

# PcoA 3: Substrate_type - Categories
plot_ordination(phyloseq_virome_transformed, phyloseq_virome_transformed_pcOA, color = "Substrate_type", shape = "Category") +
  geom_point(aes(fill = Substrate_type), size = 4.5, color = "black") +  # Use fill aesthetic and set point outline color
    theme_bw() +
    scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +  # Define fill colors
  scale_shape_manual(values = fillable_shapes) +  # Use shapes that support filling
    geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
    geom_vline(xintercept = 0, linetype = 3, colour = "black") +
    xlab("PC2 (24.8%)") +
    ylab("PC1 (25.8%)") +
    labs(color = "substrate type", shape = "category") +
  theme(axis.text.x = element_text(size = 11, colour = "black"),  # Set x-axis text size
        axis.text.y = element_text(size = 11, colour = "black"),  # Set y-axis text size
        panel.border = element_rect(colour = "black", size = 0.7))  # Set border color and thickness
####################################
# 8.5.2 Exclude diet (selection Gut, Heamolymph and whole BSF)
####################################
# Subset the phyloseq object to include only samples with Sample.type == 'gut'
head(sample_data(phyloseq_virome_transformed))
phyloseq_samples_gut <- subset_samples(phyloseq_samples_virome, !Category == 'diet')

phyloseq_virome_gut <- phyloseq_samples_gut %>% 
  microbiome::transform(transform = "log") %>%
  aggregate_rare(level = "Species", detection = 0.0005/100, prevalence = 0.00001/100)

# Distances
phyloseq_samples_gut_bray <- phyloseq::distance(phyloseq_virome_gut, method="bray")
phyloseq_samples_gut_bray_pcOA <- ordinate(phyloseq_virome_gut, method="PCoA", distance=phyloseq_samples_gut_bray)

# PcoA 1: 
fillable_shapes <- c(21, 22, 23, 24, 25)  # Reusing shapes as needed

plot_ordination(phyloseq_virome_gut, phyloseq_samples_gut_bray_pcOA, color = "Substrate_type", shape = "Anatomy") +
  geom_point(aes(shape = Anatomy, fill = Substrate_type), size = 4.5, color = "black") +  # Use fill aesthetic and set point outline color
  theme_bw() +
  scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +  # Define fill colors
  scale_shape_manual(values = fillable_shapes) +  # Use shapes that support filling
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (20.5%)") +
  ylab("PC1 (37.6%)") +
  labs(color = "Substrate type", shape = "Anatomy") +
  theme(axis.text.x = element_text(size = 11, colour = "black"),  # Set x-axis text size
        axis.text.y = element_text(size = 11, colour = "black"),  # Set y-axis text size
        panel.border = element_rect(colour = "black", size = 0.7))  # Set border color and thickness

# PcoA 2: 
plot_ordination(phyloseq_virome_gut, phyloseq_samples_gut_bray_pcOA, color = "Substrate_type") +
  geom_point(size = 4.5) +  # Set stroke to 1 to add border around points
  theme_bw() +
  scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +  # Define fill colors
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (20.5%)") +
  ylab("PC1 (37.6%)") +
  labs(color = "Substrate") +
  theme(axis.text.x = element_text(size = 11, colour = "black"),  # Set x-axis text size
        axis.text.y = element_text(size = 11, colour = "black"),  # Set y-axis text size
        panel.border = element_rect(colour = "black", size = 0.7))  # Set border color and thickness

plot_ordination(phyloseq_virome_gut, phyloseq_samples_gut_bray_pcOA, shape = "Substrate_type") +
  geom_point(aes(fill = Substrate_type), size = 4.5, color = "black") +  # Use fill aesthetic and set point outline color
  theme_bw() +
  scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +  # Define fill colors
  scale_shape_manual(values = fillable_shapes) +  # Use shapes that support filling
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (20.5%)") +
  ylab("PC1 (37.6%)") +
  labs(color = "Substrate type") +
  theme(axis.text.x = element_text(size = 11, colour = "black"),  # Set x-axis text size
        axis.text.y = element_text(size = 11, colour = "black"),  # Set y-axis text size
        panel.border = element_rect(colour = "black", size = 0.7))  # Set border color and thickness
####################################
# 9. Linear regression analysis
####################################
# 9.1 Regression analysis
####################################
# 9.1.1 Code variables
####################################
distance = "bray"
fdr = 0.05
####################################
# 9.1.2 Univariate anlysis
####################################
# A. Input
####################################
phyloseq_rarefied_phages_transformed1 <- phyloseq_virome_transformed
phyloseq_rarefied_phages_transformed1 <- phyloseq_virome_gut # This is by excluding diet
sample_data(phyloseq_rarefied_phages_transformed1)$BSF_ID <- NULL
sample_data(phyloseq_rarefied_phages_transformed1)$Sample_ID <- NULL
head(sample_data(phyloseq_rarefied_phages_transformed1)) 

GENUS <-as.data.frame(t(otu_table(phyloseq_rarefied_phages_transformed1))) #genus table
m <- data.frame(sample_data(phyloseq_rarefied_phages_transformed1))

# All
head(m)
m$Sample_ID <- NULL
####################################
# B. Rows should be ordered the same
####################################
meta

metadata_rows <- rownames(m)
genus_rows <- rownames(GENUS)

# Check if the row names are ordered the same
if (!identical(metadata_rows, genus_rows)) {
  # Reorder the metadata data frame based on the row order of the genus data frame
  m <- m[genus_rows, ]
  
  # Update the row names of the metadata data frame
  rownames(m) <- genus_rows
  
  print("The rows of metadata have been reordered to match the genus table.")
} else {
  print("The rows of metadata are already ordered the same as the genus table.")
}
####################################
# C. Univariate regression analyses
####################################
all <- c()

for (ii in 1:ncol(m)) { 
  capsc <- capscale(GENUS ~ m[,ii], distance = distance) 
  an <- anova.cca(capsc) 
  pval <- an["Pr(>F)"][[1]][[1]] 
  Fa <- an["F"][[1]][[1]] 
  r2 <- RsquareAdj(capsc)[[1]] 
  r2adj <- RsquareAdj(capsc)[[2]] 
  all <- rbind(all,cbind(Fa,r2,r2adj,pval))
}

FDR = p.adjust(all[,"pval"],method="BH")
all = cbind(all,FDR)
colnames(all) <- c("F","r2","r2adj","p-value","FDR") #generate table
row.names(all) <- colnames(m)
all <- as.data.frame(all)
all_DF <- all[rev(order(all$r2adj)),]
all_DF

# 23 of the 24 variables are significant
# View(metadata)
####################################
# D. Multivariate analysis
####################################
head(m)

sig.vars = row.names(all[all[,"FDR"] < fdr & all[,"p-value"] < 0.05,])
# sig.vars <- sig.vars[2:5]

m0 <- get_variable(phyloseq_rarefied_phages_transformed1,sig.vars) 
m <- na.exclude(m0)
m <- m0
m

GENUS0 <- t(otu_table(phyloseq_rarefied_phages_transformed1))
GENUS=GENUS0[row.names(GENUS0) %in% row.names(m),] 

attach(m)
mod0=capscale(GENUS ~ 1, distance=distance, na.action=na.exclude) # Model with intercept only
mod1=capscale(GENUS ~ . , data=m, distance=distance, na.action=na.exclude) # Model with all expl. variables
set.seed(1)
step.res<-ordiR2step(mod0, scope=formula(mod1), data=m, direction="forward", Pin = 0.05, R2scope = FALSE, pstep = 1000, trace = F)
step.res$anova
anova(step.res) # Summary table of stepwise ordination

ordiR2step.tab = data.matrix(step.res$anova[,1:5])
row.names(ordiR2step.tab) = gsub("\\+ ","",row.names(ordiR2step.tab))
ordiR2step.tab
all_DF
detach(m)
####################################
# E. Envfit (excluding diet)
####################################
# Select metadata
metadata <- data.frame(sample_data(phyloseq_virome_gut))
colnames(metadata)
metadata_patient_1 <- metadata[6]
head(metadata_patient_1)

# Create DF 
scrs <- as.data.frame(phyloseq_samples_gut_bray_pcOA$vectors[,1:2])

# Combine with metadata
scrs <- merge(scrs,metadata_patient_1,by=0,all=F)
rownames(scrs) <- scrs$Row.names
scrs$Row.names <- NULL

# Envfit 
set.seed(123)
phyloseq_rarefied_phages_transformed_1 <- data.frame(sample_data(phyloseq_virome_gut))
phyloseq_rarefied_phages_transformed_2 <- phyloseq_rarefied_phages_transformed_1[,5:6]
vf <- envfit(phyloseq_samples_gut_bray_pcOA$vectors, phyloseq_rarefied_phages_transformed_2, perm = 999)
vf

# Select factors you want to plot
spp.scrs <- as.data.frame(scores(vf, display = "factors"))
spp.scrs$Group <- rownames(spp.scrs)

spp.scrs$Group
spp.scrs <- spp.scrs[spp.scrs$Group=="Substrate_typechicken feed" |
                       spp.scrs$Group=="Substrate_typeswill",]

spp.scrs

# Figure1 with arrows
ggplot(scrs) +
  geom_point(mapping = aes(x = Axis.1, y = Axis.2, colour = Substrate_type), alpha=1, size=4.5) +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (20.5%)") +
  ylab("PC1 (37.6%)") +
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = Axis.1, y = 0, yend = Axis.2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + 
  geom_text(data = spp.scrs, aes(x = Axis.1, y = Axis.2, label = Group, colour="black"),
            size = 3)  +
  theme(axis.text.x = element_text(size = 11),  # Set x-axis text size
        axis.text.y = element_text(size = 11),  # Set y-axis text size
        panel.border = element_rect(colour = "black", size = 0.7))  # Set border color and thickness

# Figure2 with arrows
fillable_shapes <- c(21, 22, 23, 24, 25)  # Reusing shapes as needed

PcoA1 <- ggplot(scrs) +
  geom_point(mapping = aes(x = Axis.1, y = Axis.2, colour = Substrate_type, shape = Substrate_type, fill = Substrate_type), alpha=1, size=4.5, color = "black") +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +  
  scale_shape_manual(values = fillable_shapes) +  
  geom_hline(yintercept = 0, linetype = 3, colour = "black") + 
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab("PC2 (20.5%)") +
  ylab("PC1 (37.6%)") +
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = Axis.1, y = 0, yend = Axis.2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + 
  geom_text(data = spp.scrs, aes(x = Axis.1, y = Axis.2, label = Group, colour="black"),
            size = 3)  +
  theme(axis.text.x = element_text(size = 11),  # Set x-axis text size
        axis.text.y = element_text(size = 11),  # Set y-axis text size
        panel.border = element_rect(colour = "black", size = 0.7),  #,  # Set border color and thickness
      legend.position = "bottom")  # Move legend to bottom

# Boxplot adding to PcoA plot
library(ggsignif)

wilcox_result <- wilcox.test(Axis.1 ~ Substrate_type, data = scrs)
wilcox_result

ggplot(scrs, aes(x = Substrate_type, y = Axis.1, color = Substrate_type)) +
  geom_boxplot(lwd=0.8) +
  scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  theme_minimal() +
  geom_signif(comparisons = list(c("chicken feed", "swill")), annotations = c(paste("p-value =", signif(wilcox_result$p.value, digits = 2))), textsize = 5, colour = "black")

# Make arrow Better in Illustrator & Add boxplot underneath figure (let it look good)
# Add insert in PocA with F-vlaue, adjR, variable name)
# Add table in Excel with P-value
# Repeat similar for whole dataset --> Skew by diet (thus exclude --> Suppl. Figure)
####################################
# 10. Diversity analysis between substrate types
####################################
# 10.1 Observed richness 
####################################
phyloseq_samples_gut <- subset_samples(phyloseq_samples_virome, !Category == 'diet')

# Rarefy to a depth of 100,000 reads per sample
phyloseq_samples_gut_rarefied <- rarefy_even_depth(phyloseq_samples_gut, sample.size = 10000)

phyloseq_samples_gut_rarefied_genus <- phyloseq_samples_gut_rarefied %>%
  aggregate_rare(level = "Species", detection = 0.0005/100, prevalence = 0.00001/100)

# Extract the metadata
metadata <- sample_data(phyloseq_samples_gut_rarefied)

# Convert 'Substrate_type' to a factor if it's not already
metadata$Substrate_type <- as.factor(metadata$Substrate_type)

# Calculate richness for each category of 'Substrate_type'
richness <- estimate_richness(phyloseq_samples_gut_rarefied_genus, 
                              measures = c("Observed", "Simpson", "Shannon"), 
                              split = TRUE)


richness_met <- merge(metadata,richness,by=0,all=F)
median(richness_met$Observed[richness_met$Substrate_type == "chicken feed"]) # median of 6
nrow(richness_met[richness_met$Substrate_type == "chicken feed",]) # 21 samples 
median(richness_met$Observed[richness_met$Substrate_type == "swill"]) #. median of 3
nrow(richness_met[richness_met$Substrate_type == "swill",]) # 15 samples
# Total: 36 samples in total (+ 6 samples dieet revoved and thus 42 samples total, of which 3 removed with little virus, thus 45 original)
# Figure: richness
library(dplyr)
library(ggplot2)
library(ggpubr)

set.seed(123)  # Set a seed for reproducibility of jitter plot!

# Perform Wilcoxon test
wilcox_result <- wilcox.test(Observed ~ Substrate_type, data = filter(richness_met, Observed >= 0 & Observed <= 11L))
wilcox_effsize(Observed ~ Substrate_type, data = filter(richness_met, Observed >= 0 & Observed <= 11L))

print(wilcox_result)

# Plot with p-values
Rochness <- richness_met %>%
  filter(Observed >= 0 & Observed <= 11L) %>%
  ggplot() +
  aes(x = Substrate_type, y = Observed, fill = Substrate_type) +
  geom_boxplot() +
 # geom_point(position = position_jitter(width = 0.1), aes(color = Substrate_type), alpha = 1, size = 2.5, colour = "black") +  # Add jittered data points with geom_point
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +  # Match colors for jittered points
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14),  # Increase x-axis title size
        axis.title.y = element_text(size = 14),  # Increase y-axis title size
        panel.border = element_rect(colour = "black", size = 0.7)) +
  ylab("Observed richness") +
  xlab("") +
  geom_signif(comparisons = list(c("chicken feed", "swill")), 
              annotations = c(paste("p-value =", signif(wilcox_result$p.value, digits = 2))), 
              textsize = 4.5, colour = "black") +
  ylim(0,10)
####################################
# 10.2 Shannon diversity
####################################
phyloseq_samples_gut <- subset_samples(phyloseq_samples_virome, !Category == 'diet')

# Calculate richness for each category of 'Substrate_type'
richness <- estimate_richness(phyloseq_samples_gut, 
                              measures = c( "Simpson", "Shannon"), 
                              split = TRUE)


richness_met <- merge(metadata,richness,by=0,all=F)

# Figure: richness
library(dplyr)
library(ggplot2)
library(ggpubr)

# Perform Wilcoxon test for Shannon diversity
wilcox_result_shannon <- wilcox.test(Shannon ~ Substrate_type, data = filter(richness_met, !is.na(Shannon)))
wilcox_effsize(Shannon ~ Substrate_type, data = richness_met)

# Perform Wilcoxon test for Simpson diversity
wilcox_result_simpson <- wilcox.test(Simpson ~ Substrate_type, data = filter(richness_met, !is.na(Simpson)))

# Print Wilcoxon test results
print(wilcox_result_shannon)
print(wilcox_result_simpson)

# Plot Shannon diversity with p-values
set.seed(123)  # Set a seed for reproducibility of jitter plot!

Shannon <- ggplot(richness_met, aes(x = Substrate_type, y = Shannon, fill = Substrate_type)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.1), aes(color = Substrate_type), alpha = 0, size = 0, colour = "black") +  # Add jittered data points with geom_point
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +  # Match colors for jittered points
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14),  # Increase x-axis title size
        axis.title.y = element_text(size = 14),  # Increase y-axis title size
        panel.border = element_rect(colour = "black", size = 0.7)) +
  ylab("Shannon Diversity") +
  xlab("") +
  geom_signif(comparisons = list(c("chicken feed", "swill")), 
              annotations = c(paste("p-value =", signif(wilcox_result_shannon$p.value, digits = 2))), 
              textsize = 5, colour = "black")

Shannon

# Plot Simpson diversity with p-values
set.seed(123)  # Set a seed for reproducibility of jitter plot!
ggplot(richness_met, aes(x = Substrate_type, y = Simpson, fill = Substrate_type)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.1), aes(color = Substrate_type), alpha = 1, size = 2.5, colour = "black") +  # Add jittered data points with geom_point
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  scale_color_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +  # Match colors for jittered points
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14),  # Increase x-axis title size
        axis.title.y = element_text(size = 14),  # Increase y-axis title size
        panel.border = element_rect(colour = "black", size = 0.7)) +
  ylab("Simpson Diversity") +
  xlab("") +
  geom_signif(comparisons = list(c("chicken feed", "swill")), 
              annotations = c(paste("p-value =", signif(wilcox_result_simpson$p.value, digits = 2))), 
              textsize = 5, colour = "black")
####################################
# 11. Find sharing of viruses
####################################
# 11.1 Exclude Diet samples
####################################
phyloseq_samples_gut <- subset_samples(phyloseq_samples_virome, !Category == 'diet')

# Get OTU table
otu_table <- otu_table(phyloseq_samples_gut)

# Identify taxa that are present in each sample
taxa_present <- rownames(otu_table)[rowSums(otu_table) > 0]

# Identify taxa that are present in at least one sample after pruning
taxa_leftover <- taxa_present[!taxa_present %in% otu_table]

# Subset the phyloseq object to retain only the remaining taxa
phyloseq_samples_gut_leftover <- prune_taxa(taxa_leftover, phyloseq_samples_gut)
phyloseq_samples_gut_leftover
####################################
# 11.2 Sharing of viruses
####################################
phyloseq_samples_gut_leftover

# Subset the phyloseq object based on the 'Substrate_type'.
chicken_feed_samples <- subset_samples(phyloseq_samples_gut_leftover, Substrate_type == "chicken feed")
swill_samples <- subset_samples(phyloseq_samples_gut_leftover, Substrate_type == "swill")

# Extract OTU tables for each subset.
chicken_feed_otu <- otu_table(chicken_feed_samples)
swill_otu <- otu_table(swill_samples)

# Convert OTU Tables to Presence/Absence:
chicken_feed_otu_pa <- chicken_feed_otu
chicken_feed_otu_pa[chicken_feed_otu > 0] <- 1

swill_otu_pa <- swill_otu
swill_otu_pa[swill_otu > 0] <- 1

# Sum across samples
chicken_feed_counts <- rowSums(chicken_feed_otu_pa)
swill_counts <- rowSums(swill_otu_pa)

# Remove contigs not shared in any sample
chicken_feed_counts <- chicken_feed_counts[chicken_feed_counts > 0]
swill_counts <- swill_counts[swill_counts > 0]

# Tabulate the counts
chicken_feed_summary <- as.data.frame(table(chicken_feed_counts))
colnames(chicken_feed_summary) <- c("Counts", "Frequency")
chicken_feed_summary$variable <- "chicken feed"
swill_summary <- as.data.frame(table(swill_counts))
colnames(swill_summary) <- c("Counts", "Frequency")
swill_summary$variable <- "swill"
chicken_feed_swill_summary <- rbind(chicken_feed_summary,swill_summary)
#View(chicken_feed_swill_summary)

# Define the desired order of Counts
count_order <- as.character(1:16)

# Convert Counts column to factor with ordered levels
chicken_feed_swill_summary$Counts <- factor(chicken_feed_swill_summary$Counts, levels = count_order, ordered = TRUE)

# Sort the data frame by increasing values of the "Counts" column
chicken_feed_swill_summary <- chicken_feed_swill_summary %>%
  arrange(Counts)

# Figure
Sharing <- ggplot(chicken_feed_swill_summary) +
  aes(x = Counts, fill = variable, weight = Frequency) +
  geom_bar() +
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  theme_bw() +
  facet_wrap(vars(variable)) +
  ylab("Number of samples") +
  xlab("Number of viral contigs shared") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 12),  # Increase x-axis title size
        axis.title.y = element_text(size = 12),  # Increase y-axis ttle size
        panel.border = element_rect(colour = "black", size = 0.7))


chicken_feed_swill_summary
sum(chicken_feed_swill_summary$Frequency[chicken_feed_swill_summary$variable=="chicken feed"])
(6/18)*100 # 33.3% Unique viral contigs Swill 
100-((6/18)*100) # 66.7% Unique viral contigs Swill 

sum(chicken_feed_swill_summary$Frequency[chicken_feed_swill_summary$variable=="swill"])
(16/25)*100 # 64% Unique viral contigs Swill 
100-((16/25)*100) # 36% Unique viral contigs Swill 

#View(Mastertable_HQ_viral)
#view(chicken_feed_counts)
#view(swill_counts)
####################################
# 11.3  Differential abundance analyses
####################################
phyloseq_samples_gut <- subset_samples(phyloseq_samples_virome, !Category == 'diet')

# Rarefy to a depth of 100,000 reads per sample
phyloseq_samples_gut_rarefied <- rarefy_even_depth(phyloseq_samples_gut, sample.size = 10000)

phyloseq_samples_gut_rarefied_genus <- phyloseq_samples_gut_rarefied %>%
  aggregate_rare(level = "Species", detection = 0.0005/100, prevalence = 0.00001/100)

# Lefse: Disease
mm_lefse <- run_lefse(
  phyloseq_samples_gut_rarefied_genus,
  group = "Substrate_type",
  transform = c("identity"),
  taxa_rank = "Species",
  kw_cutoff = 0.05,
  wilcoxon_cutoff = 0.05,
  strict = "1",
  lda_cutoff = 2)

mm_lefse_table_Disease <- marker_table(mm_lefse)
plot_ef_bar(mm_lefse_table_Disease, label_level = 0) +
  scale_fill_manual(values = c("swill" = "#E6372C", "chicken feed" = "#499B78"))
####################################

######### !!!!!!!!!!!!!!! #########
    # SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########
