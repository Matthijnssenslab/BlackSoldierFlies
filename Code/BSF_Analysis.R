####################################
# SCRIPT: BSF ANALYSES
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
setwd("/Users/daan/Desktop/Transfer/Input/Abundances")
getwd()
dir()

## B) read the abundance text files and store in a list
list <- list()
list_txt <- dir(pattern = "*.magnitudes", full.names = FALSE)
for (k in 1: length(list_txt)){
  list[[k]] <- read.delim(list_txt[k])
}

## C) create names file
names <- sub(pattern = "*.magnitudes", "", list_txt)

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
setwd("/Users/daan/Desktop/Transfer/Input/Taxonomy")
getwd()
dir()

AllSamples_classificationtaxonomy_genomead <- as.data.frame(read_delim("Genomead.txt", "\t", escape_double = FALSE, col_names = TRUE, trim_ws = FALSE))
AllSamples_classificationtaxonomy_genomead[is.na(AllSamples_classificationtaxonomy_genomead)] <- "Unclassified"
AllSamples_classificationtaxonomy_genomead <- AllSamples_classificationtaxonomy_genomead[AllSamples_classificationtaxonomy_genomead$NODE %in% contigs_project,]
rownames(AllSamples_classificationtaxonomy_genomead) <- AllSamples_classificationtaxonomy_genomead$NODE
AllSamples_classificationtaxonomy_genomead$NODE <- NULL
rownames(AllSamples_classificationtaxonomy_genomead) <- gsub("\\.{1}", "_", rownames(AllSamples_classificationtaxonomy_genomead))
#View(AllSamples_classificationtaxonomy_genomead)
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
setwd("/Users/daan/Desktop/Transfer/Input/Host")
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
setwd("/Users/daan/Desktop/Transfer/Input/Function")
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
setwd("/Users/daan/Desktop/Transfer/Input/Taxonomy")
getwd() 
dir()

AllsamplesGenus <- as.data.frame(read_delim("genus_clusters2.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
colnames(AllsamplesGenus) <- c('contig', 'Genus')
rownames(AllsamplesGenus) <- AllsamplesGenus$contig
AllsamplesGenus$contig <- NULL
rownames(AllsamplesGenus) <- gsub("\\.{1}", "_", rownames(AllsamplesGenus))
#View(AllsamplesGenus)
####################################
# 7. Create mastertable: merge scaffolds, abundances and annotations
####################################
COMBO_1 <- merge(table, AllSamples_classificationtaxonomy_genomead, by = 0, all = F)
rownames(COMBO_1) <- COMBO_1$Row.names
COMBO_1$Row.names <- NULL

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
# 7.1 Additional information to Mastertable
####################################
# 7.1.1 Insert total number of trimmed reads
####################################
Mastertable[is.na(Mastertable)] <- 0
Mastertable$Totalnumberofreads <- as.numeric(rowSums(Mastertable[,colnames(Mastertable) %in% names]))
Mastertable <- Mastertable[!(Mastertable$Totalnumberofreads == "0"),] ## remove rows with total number of zero
####################################
# 7.1.2 Insert length of contigs
####################################
Mastertable$length <- rownames(Mastertable)
Mastertable$length <- gsub("*_cov.*","",Mastertable$length)
Mastertable$length <- gsub(".*length_","", Mastertable$length)
Mastertable$length <- as.numeric(Mastertable$length)
Mastertable$length <- Mastertable$length/1000
ncol(Mastertable)
####################################
# 7.1.2 viral identification and completeness correction
####################################
Mastertable_viral <- Mastertable[Mastertable$Virus=="phage" | Mastertable$Virus=='eukaryotic_virus',]
Mastertable_HQ_viral <- Mastertable_viral[Mastertable_viral$completeness > 50,]
####################################
# 8. Phyloseq
####################################
# 8.1 Metadata
####################################
library(readxl)
setwd("/Users/daan/Desktop/Transfer/Input/metadata")
dir()
metadata <- as.data.frame(read_excel("Metadata2.xlsx", sheet = "metadata"))
rownames(metadata) <- metadata$Sample_identifier
metadata[metadata == "NA"] <- NA
####################################
# 8.2 Abundance table
####################################
# Identify sample columns (assuming they're the first 67 columns)
sample_cols <- colnames(Mastertable_HQ_viral)[1:67]

# Calculate read sums for each sample
read_sums <- colSums(Mastertable_HQ_viral[, sample_cols])
read_sums

# Identify outliers (samples with less than 100,000 reads)
outliers <- names(read_sums[read_sums < 100000])

print("Samples removed as outliers:")
print(outliers)

# Remove outlier samples
Mastertable_HQ_viral <- Mastertable_HQ_viral[, !(colnames(Mastertable_HQ_viral) %in% outliers)]

# Recalculate sample columns after removal of outliers
sample_cols <- colnames(Mastertable_HQ_viral)[1:(67 - length(outliers))]

# Create abundance table
abundance_table <- Mastertable_HQ_viral[, sample_cols]
abundance_table$names <- rownames(abundance_table)

# Aggregate data
abundance_table1 <- aggregate(. ~ names, FUN = sum, data = abundance_table)
rownames(abundance_table1) <- abundance_table1$names
abundance_table1$names <- NULL

# Print summary
print(paste("Number of samples retained:", ncol(abundance_table1)))
print("Sample names retained:")
print(colnames(abundance_table1))
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
abundance_table_m <- as.matrix(abundance_table1)
taxonomy_table_m <- as.matrix(Taxonomy_table)

ABUNDANCE_rarefied_family <- otu_table(abundance_table_m, taxa_are_rows = T)
TAX_rarefied_family <- tax_table(taxonomy_table_m)
samples <- sample_data(metadata)

phyloseq_samples <- phyloseq(ABUNDANCE_rarefied_family,TAX_rarefied_family, samples)
phyloseq_samples
colnames(otu_table(phyloseq_samples))
####################################
# 8.5 Subset phyloseq table
####################################
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define the names you want to exclude
names_to_exclude <- c("Hae_CF", "Hae_SW", "Hae6", "Swill6_4", "Swill6_5", "Swill6_6",
                      "WL_CF_5", "WL_CF_Small", "WL_SW_5", "WL_SW_Small")

# Create the final subset
phyloseq_virome <- subset_samples(phyloseq_samples, 
                                  !(Sample_identifier %in% names_to_exclude))

# Print summary of the subsetted dataset
print(phyloseq_virome)
print(paste("Number of samples in phyloseq_virome:", nsamples(phyloseq_virome)))
print(paste("Number of taxa in phyloseq_virome:", ntaxa(phyloseq_virome)))

# Optional: View the sample data of the subsetted dataset
View(sample_data(phyloseq_virome))
####################################
# 8.6 Visualize virome
####################################
# 8.6.1 Complete dataset
####################################
phyloseq_samples_virome <- phyloseq_virome

# Hellinger transformation and rare taxa aggregation
phyloseq_virome_transformed <- phyloseq_samples_virome %>%
  microbiome::transform(transform = "hellinger") %>%
  aggregate_rare(level = "Species", detection = 0.05/100, prevalence = 0.01/100)

# Calculate distances and perform PCoA
phyloseq_virome_transformed_bray <- phyloseq::distance(phyloseq_virome_transformed, method="bray")
phyloseq_virome_transformed_pcOA <- ordinate(phyloseq_virome_transformed, method="PCoA", distance=phyloseq_virome_transformed_bray)

# Extract variance explained by each axis
variance_explained <- phyloseq_virome_transformed_pcOA$values$Relative_eig * 100
pc1_var <- round(variance_explained[1], 1)
pc2_var <- round(variance_explained[2], 1)

# Define color palette and theme
color_palette <- c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")
plot_theme <- theme_bw() +
  theme(axis.text = element_text(size = 11, colour = "black"),
        panel.border = element_rect(colour = "black", size = 0.7),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# PCoA 1: Substrate_type
pcoa_plot1 <- plot_ordination(phyloseq_virome_transformed, phyloseq_virome_transformed_pcOA, color = "Substrate_type") +
  geom_point(alpha = 1, size = 4.5) +
  scale_color_manual(values = color_palette) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") +
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  labs(x = paste0("PC2 (", pc2_var, "%)"),
       y = paste0("PC1 (", pc1_var, "%)"),
       color = "Substrate") +
  plot_theme

# PCoA 2: Substrate_type and Location
fillable_shapes <- c(21, 22, 23, 24, 25)
pcoa_plot2 <- plot_ordination(phyloseq_virome_transformed, phyloseq_virome_transformed_pcOA, shape = "Location") +
  geom_point(aes(fill = Substrate_type), size = 4.5, color = "black") +
  scale_fill_manual(values = color_palette) +
  scale_shape_manual(values = fillable_shapes) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") +
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  labs(x = paste0("PC2 (", pc2_var, "%)"),
       y = paste0("PC1 (", pc1_var, "%)"),
       fill = "Substrate",
       shape = "Location") +
  plot_theme

# Display plots
print(pcoa_plot1)
print(pcoa_plot2)
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
library(vegan)

phyloseq_rarefied_phages_transformed1 <- phyloseq_virome_transformed
head(sample_data(phyloseq_rarefied_phages_transformed1)) 

GENUS <-as.data.frame(t(otu_table(phyloseq_rarefied_phages_transformed1))) #genus table
m <- data.frame(sample_data(phyloseq_rarefied_phages_transformed1))
m$Sample_type <- NULL
m$Black_soldier_fly_identifier <- NULL
m$Category <- NULL
head(m) # when documenting variation, exclude from metadata larvae age

# Function to perform univariate regression analyses
perform_univariate_analysis <- function(GENUS, m) {
  all <- c()
  for (ii in 1:ncol(m)) {
    # Check if the column has more than one unique value
    if (length(unique(m[,ii])) > 1) {
      tryCatch({
        capsc <- capscale(GENUS ~ m[,ii], distance = "bray")
        an <- anova.cca(capsc)
        pval <- an["Pr(>F)"][[1]][[1]]
        Fa <- an["F"][[1]][[1]]
        r2 <- RsquareAdj(capsc)[[1]]
        r2adj <- RsquareAdj(capsc)[[2]]
        all <- rbind(all, cbind(Fa, r2, r2adj, pval))
      }, error = function(e) {
        # If there's an error, add a row of NAs
        all <- rbind(all, cbind(NA, NA, NA, NA))
      })
    } else {
      # If the column has only one unique value, add a row of NAs
      all <- rbind(all, cbind(NA, NA, NA, NA))
    }
  }
  if (nrow(all) > 0) {
    FDR <- p.adjust(all[,"pval"], method="BH")
    all <- cbind(all, FDR)
    colnames(all) <- c("F", "r2", "r2adj", "p-value", "FDR")
    rownames(all) <- colnames(m)
    all <- as.data.frame(all)
    all_DF <- all[rev(order(all$r2adj)),]
    return(all_DF)
  } else {
    return(NULL)
  }
}

# Perform analysis for all data
all_DF <- perform_univariate_analysis(GENUS, m)
print("All data:")
print(all_DF)

# Filter data for Belgium
belgium_samples <- m$Location == "Belgium"
GENUS_belgium <- GENUS[belgium_samples, ]
m_belgium <- m[belgium_samples, ]
m_belgium$Location <- NULL  # Remove Location column as it's constant

# Perform analysis for Belgium
belgium_DF <- perform_univariate_analysis(GENUS_belgium, m_belgium)
print("Belgium data:")
print(belgium_DF)

# Filter data for Italy
italy_samples <- m$Location == "Italy"
GENUS_italy <- GENUS[italy_samples, ]
m_italy <- m[italy_samples, ]
m_italy$Location <- NULL  # Remove Location column as it's constant

# Perform analysis for Italy
italy_DF <- perform_univariate_analysis(GENUS_italy, m_italy)
print("Italy data:")
print(italy_DF)

# Filter data for excluding whole gut
excluding_whole_gut <- !m$Anatomy == "whole gut"
GENUS_excluded_whole_gut <- GENUS[excluding_whole_gut, ]
m_with_whole_gut <- m[excluding_whole_gut, ]
m_with_whole_gut$Location <- NULL
m_with_whole_gut$Larval_age_DAH <- NULL

# Perform analysis for BSF-samples without 'whole gut' samples (difference between migdut & hindgut)
with_whole_gut_DF <- perform_univariate_analysis(GENUS_excluded_whole_gut, m_with_whole_gut)
print("withut whole gut data:")
print(with_whole_gut_DF)
# Univariate variation increases from 26% to 35%

####################################
# D. Multivariate analysis
####################################
head(m)

sig.vars = row.names(all_DF[all_DF[,"FDR"] < fdr & all_DF[,"p-value"] < 0.05,])
sig.vars <- sig.vars[6:9] # Exclude BSF identfier, since it absorbs all variation (and is known in virome analyses)
sig.vars

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
# E. Create barplot of univariate and multivariate R2
####################################
# Prepare data for plotting
all_DF # Univariate 
ordiR2step.tab # Multivariate

univariate_r2 <- all_DF[all_DF$FDR < fdr & all_DF$`p-value` < 0.05, ]
univariate_r2$Variable <- rownames(univariate_r2)
univariate_r2$Analysis <- "Univariate"
univariate_r2$r2 <- univariate_r2$r2adj

# Make sure to use the correct column name for R2 in ordiR2step.tab
multivariate_r2 <- data.frame(
  Variable = rownames(ordiR2step.tab),
  r2 = ordiR2step.tab[, "R2.adj"],  # Make sure this column exists
  Analysis = "Multivariate"
)

# Combine and filter data
combined_r2 <- rbind(
  univariate_r2[, c("Variable", "r2", "Analysis")],
  multivariate_r2
) %>%
  group_by(Variable) %>%
  filter(n() == 2) %>%
  ungroup()

# Order variables based on multivariate R2
var_order <- combined_r2 %>%
  filter(Analysis == "Multivariate") %>%
  arrange(desc(r2)) %>%
  pull(Variable)

combined_r2$Variable <- factor(combined_r2$Variable, levels = var_order)
combined_r2
# Convert R2 to percentage
combined_r2$r2_percent <- combined_r2$r2 * 100

# Create the barplot
r2_barplot <- ggplot(combined_r2, aes(x = Variable, y = r2_percent, fill = Analysis)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Multivariate" = "grey20", "Univariate" = "darkgrey")) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = expression(R^2 * " / Effect size (%)"),
    title = expression("Univariate and Multivariate " * R^2 * " Values")
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(color = "black", size = 10),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.spacing = unit(0.5, "lines")
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(combined_r2$r2_percent) * 1.1))

print(r2_barplot)

# Save the plot using ggsave
ggsave("~/Desktop/r2_barplot.pdf", plot = r2_barplot, width = 10, height = 8, device = "pdf")
####################################
# F. Envfit
####################################
# Select metadata
metadata <- data.frame(sample_data(phyloseq_virome_transformed))
colnames(metadata)
metadata_patient_1 <- metadata[4]
head(metadata_patient_1)

# Create DF 
scrs <- as.data.frame(phyloseq_virome_transformed_pcOA$vectors[,1:2])

# Combine with metadata
scrs <- merge(scrs,metadata_patient_1,by=0,all=F)
rownames(scrs) <- scrs$Row.names
scrs$Row.names <- NULL

# Envfit 
set.seed(123)
phyloseq_rarefied_phages_transformed_1 <- data.frame(sample_data(phyloseq_virome_transformed))
phyloseq_rarefied_phages_transformed_2 <- phyloseq_rarefied_phages_transformed_1[,4:5]
vf <- envfit(phyloseq_virome_transformed_pcOA$vectors, phyloseq_rarefied_phages_transformed_2, perm = 999)
vf

# Select factors you want to plot
spp.scrs <- as.data.frame(scores(vf, display = "factors"))
spp.scrs$Group <- rownames(spp.scrs)

spp.scrs$Group
spp.scrs <- spp.scrs[spp.scrs$Group=="Substrate_typechicken feed" |
                       spp.scrs$Group=="Substrate_typeswill",]

# Assuming 'metadata' contains the full metadata information
# and 'scrs' contains the PCoA coordinates

# Add Anatomy2 and Location information to scrs
scrs$Anatomy2 <- metadata$Anatomy2[match(rownames(scrs), rownames(metadata))]
scrs$Location <- metadata$Location[match(rownames(scrs), rownames(metadata))]

# Get unique Anatomy2 values
unique_anatomies <- unique(scrs$Anatomy2)

# Create a named vector for shape mapping (using filled shapes)
anatomy_shapes <- setNames(21:25, unique_anatomies)  # Adjust range if more shapes are needed

# Calculate the proportion of variance explained by each axis
variance_explained <- phyloseq_virome_transformed_pcOA$values$Relative_eig * 100

# Create custom polygons for each location
create_polygon <- function(data, location) {
  subset_data <- data[data$Location == location, ]
  hull <- chull(subset_data$Axis.1, subset_data$Axis.2)
  return(subset_data[hull, ])
}

belgium_polygon <- create_polygon(scrs, "Belgium")
italy_polygon <- create_polygon(scrs, "Italy")


# Create the plot
pcoa_plot <- ggplot(scrs, aes(x = Axis.1, y = Axis.2)) +
  # Add custom polygons for Location
  geom_polygon(data = belgium_polygon, aes(fill = Location), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = italy_polygon, aes(fill = Location), alpha = 0.2, show.legend = FALSE) +
  geom_point(aes(fill = Substrate_type, shape = Anatomy2),
             color = "black",  # Set border color to black
             alpha = 1, size = 6) +
  coord_fixed() +
  theme_bw() +
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92",
                               "Belgium" = "#e41a1c", "Italy" = "#4daf4a")) +
  scale_shape_manual(values = anatomy_shapes) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") +
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab(paste0("PC1 (", round(variance_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  geom_text(data = spp.scrs,
            aes(x = Axis.1, y = Axis.2, label = Group),
            colour = "black", size = 3, vjust = -0.5) +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    panel.border = element_rect(colour = "black", size = 0.7),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(
    fill = guide_legend(title = "Substrate type", override.aes = list(shape = 21)),
    shape = guide_legend(title = "Anatomy2")
  )

# Display the plot
print(pcoa_plot)

# Create the plot
pcoa_plot <- ggplot(scrs, aes(x = Axis.1, y = Axis.2)) +
  # Add custom polygons for Location
  geom_polygon(data = belgium_polygon, aes(fill = Location), alpha = 0.2, show.legend = FALSE) +
  geom_polygon(data = italy_polygon, aes(fill = Location), alpha = 0.2, show.legend = FALSE) +
  geom_point(aes(fill = Substrate_type, shape = Anatomy2),
             color = "black",  # Set border color to black
             alpha = 1, size = 4.5) +
  coord_fixed() +
  theme_bw() +
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92",
                               "Belgium" = "#e41a1c", "Italy" = "#4daf4a")) +
  scale_shape_manual(values = anatomy_shapes) +
  geom_hline(yintercept = 0, linetype = 3, colour = "black") +
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  xlab(paste0("PC1 (", round(variance_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    panel.border = element_rect(colour = "black", size = 0.7),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(
    fill = guide_legend(title = "Substrate type", override.aes = list(shape = 21)),
    shape = guide_legend(title = "Anatomy2")
  )

# Display the plot
print(pcoa_plot)

# Save the plot
ggsave("~/Desktop/pcoa_plot_with_anatomy_location_custom.pdf", pcoa_plot, width = 12, height = 10, device = "pdf")
####################################
# 10. Diversity analysis between substrate types
####################################
# 10.1 Observed richness 
####################################
library(phyloseq)
library(ggplot2)
library(dplyr)
library(rstatix)

# Subset the phyloseq object
phyloseq_samples_gut <- subset_samples(phyloseq_samples_virome, Category != 'diet')

# Rarefy to a depth of 100,000 reads per sample
phyloseq_samples_gut_rarefied <- rarefy_even_depth(phyloseq_samples_gut, sample.size = 100000, rngseed = 123)

# Aggregate rare taxa
phyloseq_samples_gut_rarefied_genus <- phyloseq_samples_gut_rarefied %>%
  aggregate_rare(level = "Species", detection = 0.0005/100, prevalence = 0.00001/100)

# Calculate richness
richness <- estimate_richness(phyloseq_samples_gut_rarefied_genus,
                              measures = c("Observed", "Simpson", "Shannon"))

# Combine richness with metadata
richness_met <- data.frame(sample_data(phyloseq_samples_gut_rarefied_genus), richness)

# Perform statistical tests (only within-location)
within_location_tests <- richness_met %>%
  group_by(Location) %>%
  wilcox_test(Observed ~ Substrate_type) %>%
  add_significance()

# Apply BH correction
all_tests <- within_location_tests %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Create the plot
richness_plot <- ggplot(richness_met, aes(x = Location, y = Observed, fill = Substrate_type)) +
  geom_boxplot(width = 0.7, position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14),
        panel.border = element_rect(colour = "black", size = 0.7),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(y = "Observed richness", x = "Location", fill = "Substrate type")

# Add statistical annotations
max_y <- max(richness_met$Observed)
y_offset <- max_y * 0.05

for (loc in c("Belgium", "Italy")) {
  test_result <- all_tests %>% filter(Location == loc)
  p_label <- paste("p =", round(test_result$p.adj, 3))
  
  # Calculate positions for the horizontal line
  x_pos <- ifelse(loc == "Belgium", 0.7, 1.7)
  y_pos <- max_y + y_offset
  
  richness_plot <- richness_plot +
    annotate("text", x = x_pos + 0.3, y = y_pos + y_offset * 0.5, label = p_label, size = 3) +
    annotate("segment", x = x_pos, xend = x_pos + 0.6, y = y_pos, yend = y_pos)
}

# Adjust the y-axis limit to accommodate the annotations
richness_plot <- richness_plot + 
  coord_cartesian(ylim = c(0, max_y * 1.2))

# Display the plot
print(richness_plot)

# Print all statistical test results
print(all_tests)

# Calculate median richness and sample size for each group
summary_stats <- aggregate(Observed ~ Location + Substrate_type, data = richness_met,
                           FUN = function(x) c(Median_richness = median(x), Sample_size = length(x)))

# Convert the result to a data frame with separate columns
summary_stats_df <- data.frame(
  Location = summary_stats$Location,
  Substrate_type = summary_stats$Substrate_type,
  Median_richness = summary_stats$Observed[,1],
  Sample_size = summary_stats$Observed[,2]
)

# Print the summary statistics
print(summary_stats_df)

# Save the plot using ggsave
ggsave("~/Desktop/richness_plot.pdf", plot = richness_plot, width = 10, height = 8, device = "pdf")
####################################
# 10.2 Shannon diversity
####################################
library(phyloseq)
library(ggplot2)
library(dplyr)
library(rstatix)

# Subset the phyloseq object
phyloseq_samples_gut <- subset_samples(phyloseq_samples_virome, Category != 'diet')

# Calculate richness for each category of 'Substrate_type'
richness <- estimate_richness(phyloseq_samples_gut,
                              measures = c("Simpson", "Shannon"),
                              split = TRUE)

shannon_met <- merge(sample_data(phyloseq_samples_gut), richness, by = 0, all = F)

# Perform statistical tests
# Within location comparisons only
within_location_tests <- shannon_met %>%
  group_by(Location) %>%
  wilcox_test(Shannon ~ Substrate_type) %>%
  add_significance()

# Apply BH correction
all_tests <- within_location_tests %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Create the base plot
shannon_plot <- shannon_met %>%
  ggplot(aes(x = Location, y = Shannon, fill = Substrate_type)) +
  geom_boxplot(width = 0.7, position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(colour = "black", size = 0.7),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  ylab("Shannon diversity") +
  xlab("Location") +
  labs(fill = "Substrate type")

median(shannon_met$Shannon[shannon_met$Location=="Belgium"])
median(shannon_met$Shannon[shannon_met$Location=="Belgium" & shannon_met$Substrate_type=="chicken feed"])
median(shannon_met$Shannon[shannon_met$Location=="Belgium" & shannon_met$Substrate_type=="swill"])

median(shannon_met$Shannon[shannon_met$Location=="Italy"])
median(shannon_met$Shannon[shannon_met$Location=="Italy" & shannon_met$Substrate_type=="chicken feed"])
median(shannon_met$Shannon[shannon_met$Location=="Italy" & shannon_met$Substrate_type=="swill"])

# Add statistical annotations
max_y <- max(shannon_met$Shannon)
y_offset <- max_y * 0.05

# Within location comparisons
for (loc in c("Belgium", "Italy")) {
  test_result <- all_tests %>% filter(Location == loc)
  p_label <- paste("p =", round(test_result$p.adj, 3))
  
  # Calculate positions for the horizontal line
  x_pos <- ifelse(loc == "Belgium", 0.7, 1.7)
  y_pos <- max_y + y_offset
  
  shannon_plot <- shannon_plot +
    annotate("text", x = x_pos + 0.3, y = y_pos + y_offset * 0.5, label = p_label, size = 3) +
    annotate("segment", x = x_pos, xend = x_pos + 0.6, y = y_pos, yend = y_pos)
}

# Adjust the y-axis limit to accommodate the annotations
shannon_plot <- shannon_plot + 
  coord_cartesian(ylim = c(0, max_y * 1.2))

# Display the plot
print(shannon_plot)

# Print all statistical test results
print(all_tests)

# Save the plot using ggsave
ggsave("~/Desktop/Shannon_diversity.pdf", plot = shannon_plot, width = 10, height = 8, device = "pdf")
####################################
# 10.3 Create Richness and Diversity plot
####################################
library(phyloseq)
library(ggplot2)
library(dplyr)
library(rstatix)
library(patchwork)
library(ggtext)

# Subset the phyloseq object
phyloseq_samples_gut <- subset_samples(phyloseq_samples_virome, Category != 'diet')

# Rarefy to a depth of 100,000 reads per sample
phyloseq_samples_gut_rarefied <- rarefy_even_depth(phyloseq_samples_gut, sample.size = 100000, rngseed = 123)

# Aggregate rare taxa
phyloseq_samples_gut_rarefied_genus <- phyloseq_samples_gut_rarefied %>%
  aggregate_rare(level = "Species", detection = 0.0005/100, prevalence = 0.00001/100)

# Calculate richness
richness <- estimate_richness(phyloseq_samples_gut_rarefied_genus,
                              measures = c("Observed", "Simpson", "Shannon"))

# Combine richness with metadata
richness_met <- data.frame(sample_data(phyloseq_samples_gut_rarefied_genus), richness)

# Function to create plot
create_diversity_plot <- function(data, measure, y_label) {
  plot <- ggplot(data, aes(x = Location, y = !!sym(measure), fill = Substrate_type)) +
    geom_boxplot(width = 0.7, position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")) +
    theme_bw() +
    theme(legend.position = "right",
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14),
          panel.border = element_rect(colour = "black", size = 0.7),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.title.x = element_blank()) +  # Remove x-axis title
    labs(y = y_label, fill = "Substrate type")
  
  # Perform statistical tests
  within_location_tests <- data %>%
    group_by(Location) %>%
    do(wilcox_test(data = ., as.formula(paste(measure, "~ Substrate_type")))) %>%
    add_significance() %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  
  # Between-location test
  between_location_test <- wilcox_test(data = data, as.formula(paste(measure, "~ Location")))
  
  # Add statistical annotations
  max_y <- max(data[[measure]])
  y_offset <- max_y * 0.05
  
  for (loc in c("Belgium", "Italy")) {
    test_result <- within_location_tests %>% filter(Location == loc)
    p_label <- paste("p =", format.pval(test_result$p.adj, digits = 3))
    
    x_pos <- ifelse(loc == "Belgium", 0.7, 1.7)
    y_pos <- max_y + y_offset
    
    plot <- plot +
      annotate("text", x = x_pos + 0.3, y = y_pos + y_offset * 0.5, label = p_label, size = 3) +
      annotate("segment", x = x_pos, xend = x_pos + 0.6, y = y_pos, yend = y_pos)
  }
  
  # Add between-location p-value
  between_p_label <- paste("p =", format.pval(between_location_test$p, digits = 3))
  plot <- plot +
    annotate("text", x = 1.5, y = max_y + 3*y_offset, label = between_p_label, size = 3) +
    annotate("segment", x = 1, xend = 2, y = max_y + 2.5*y_offset, yend = max_y + 2.5*y_offset)
  
  # Adjust the y-axis limit to accommodate the annotations
  plot <- plot + coord_cartesian(ylim = c(0, max_y * 1.3))
  
  return(plot)
}

# Create individual plots
richness_plot <- create_diversity_plot(richness_met, "Observed", "Observed richness") +
  ggtitle("Observed Richness") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

shannon_plot <- create_diversity_plot(richness_met, "Shannon", "Shannon diversity") +
  ggtitle("Shannon Diversity") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Remove legend from richness plot (we'll use the legend from shannon plot)
richness_plot <- richness_plot + theme(legend.position = "none")

# Combine plots using patchwork
combined_plot <- richness_plot + shannon_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Diversity Measures by Location and Substrate Type",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16)
    )
  )

# Display the combined plot
print(combined_plot)

# Print p-values
# Within Location
print("P-values for Observed Richness:")
print(wilcox_test(data = richness_met, Observed ~ Location))
median(richness_met$Observed[richness_met$Location=="Italy"])
median(richness_met$Observed[richness_met$Location=="Belgium"])


print("P-values for Shannon Diversity:")
print(wilcox_test(data = richness_met, Shannon ~ Location))

# Between location
within_location_tests

# Save the plot if needed
ggsave("~/Desktop/combined_diversity_plot.pdf", plot = combined_plot, width = 12, height = 6, device = "pdf")
####################################
# 10.4 Create Richness and Diversity plot (without whole gut) - variable Anatomy2
####################################
# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(rstatix)
library(patchwork)
library(ggtext)

# Subset the phyloseq object
phyloseq_samples_gut <- subset_samples(phyloseq_samples_virome, Anatomy != 'whole gut')

# Rarefy to a depth of 100,000 reads per sample
phyloseq_samples_gut_rarefied <- rarefy_even_depth(phyloseq_samples_gut, sample.size = 100000, rngseed = 123)

# Aggregate rare taxa
phyloseq_samples_gut_rarefied_genus <- phyloseq_samples_gut_rarefied %>%
  aggregate_rare(level = "Species", detection = 0.0005/100, prevalence = 0.00001/100)

# Calculate richness
richness <- estimate_richness(phyloseq_samples_gut_rarefied_genus,
                              measures = c("Observed", "Simpson", "Shannon"))

# Combine richness with metadata
richness_met <- data.frame(sample_data(phyloseq_samples_gut_rarefied_genus), richness)

# Reorder Anatomy2 levels
richness_met$Anatomy2 <- factor(richness_met$Anatomy2, 
                                levels = c("anterior midgut", 
                                           "middle midgut", 
                                           "posterior midgut",
                                           "hindgut"))

# Function to create plot
create_diversity_plot <- function(data, measure, y_label) {
  plot <- ggplot(data, aes(x = Anatomy2, y = !!sym(measure), fill = Anatomy2)) +
    geom_boxplot(width = 0.7, alpha = 0.7) +  # Added alpha for transparency
    scale_fill_manual(values = c("anterior midgut" = "#4D4D4D",    # Dark grey
                                 "middle midgut" = "#737373",        # Medium dark grey
                                 "posterior midgut" = "#A6A6A6",     # Medium light grey
                                 "hindgut" = "#CCCCCC")) +          # Light grey
    theme_bw() +
    theme(legend.position = "none",  # Remove legend since colors are just for visualization
          axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14),
          panel.border = element_rect(colour = "black", size = 0.7),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_line(color = "grey95"),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.title.x = element_blank()) +  # Remove x-axis title
    labs(y = y_label)
  
  # Perform Kruskal-Wallis test
  kw_test <- kruskal.test(as.formula(paste(measure, "~ Anatomy2")), data = data)
  
  # Add statistical annotation
  max_y <- max(data[[measure]])
  y_offset <- max_y * 0.05
  
  # Add p-value
  p_label <- paste("p =", format.pval(kw_test$p.value, digits = 3))
  plot <- plot +
    annotate("text", x = length(unique(data$Anatomy2))/2, 
             y = max_y + 2*y_offset, 
             label = p_label, size = 3)
  
  # Adjust the y-axis limit to accommodate the annotation
  plot <- plot + coord_cartesian(ylim = c(0, max_y * 1.2))
  
  return(plot)
}

# Create individual plots
richness_plot <- create_diversity_plot(richness_met, "Observed", "Observed richness") +
  ggtitle("Observed Richness") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

shannon_plot <- create_diversity_plot(richness_met, "Shannon", "Shannon diversity") +
  ggtitle("Shannon Diversity") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Combine plots using patchwork
combined_plot <- richness_plot + shannon_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Diversity Measures by Gut Region",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16)
    )
  )

# Display the combined plot
print(combined_plot)

# Print p-values
print("Kruskal-Wallis test results for Observed Richness:")
print(kruskal.test(Observed ~ Anatomy2, data = richness_met))

print("Kruskal-Wallis test results for Shannon Diversity:")
print(kruskal.test(Shannon ~ Anatomy2, data = richness_met))
####################################
# 10.5 Create Richness and Diversity plot (without whole gut) - variable Anatomy
####################################
# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(rstatix)
library(patchwork)
library(ggtext)

# Subset the phyloseq object
phyloseq_samples_gut <- subset_samples(phyloseq_samples_virome, Anatomy != 'whole gut')

# Rarefy to a depth of 100,000 reads per sample
phyloseq_samples_gut_rarefied <- rarefy_even_depth(phyloseq_samples_gut, sample.size = 100000, rngseed = 123)

# Aggregate rare taxa
phyloseq_samples_gut_rarefied_genus <- phyloseq_samples_gut_rarefied %>%
  aggregate_rare(level = "Species", detection = 0.0005/100, prevalence = 0.00001/100)

# Calculate richness
richness <- estimate_richness(phyloseq_samples_gut_rarefied_genus,
                              measures = c("Observed", "Simpson", "Shannon"))

# Combine richness with metadata
richness_met <- data.frame(sample_data(phyloseq_samples_gut_rarefied_genus), richness)

# Reorder Anatomy levels
richness_met$Anatomy <- factor(richness_met$Anatomy, 
                               levels = c("midgut", "hindgut"))

# Function to create plot
create_diversity_plot <- function(data, measure, y_label) {
  plot <- ggplot(data, aes(x = Anatomy, y = !!sym(measure), fill = Anatomy)) +
    geom_boxplot(width = 0.7, alpha = 0.7) +  # Added alpha for transparency
    scale_fill_manual(values = c("midgut" = "#737373",     # Medium grey
                                 "hindgut" = "#CCCCCC")) +   # Light grey
    theme_bw() +
    theme(legend.position = "none",  # Remove legend since colors are just for visualization
          axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14),
          panel.border = element_rect(colour = "black", size = 0.7),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_line(color = "grey95"),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.title.x = element_blank()) +  # Remove x-axis title
    labs(y = y_label)
  
  # Perform Wilcoxon test (since we now have only 2 groups)
  wilcox_test_result <- wilcox.test(as.formula(paste(measure, "~ Anatomy")), data = data)
  
  # Add statistical annotation
  max_y <- max(data[[measure]])
  y_offset <- max_y * 0.05
  
  # Add p-value
  p_label <- paste("p =", format.pval(wilcox_test_result$p.value, digits = 3))
  plot <- plot +
    annotate("text", x = 1.5, 
             y = max_y + 2*y_offset, 
             label = p_label, size = 3)
  
  # Adjust the y-axis limit to accommodate the annotation
  plot <- plot + coord_cartesian(ylim = c(0, max_y * 1.2))
  
  return(plot)
}

# Create individual plots
richness_plot <- create_diversity_plot(richness_met, "Observed", "Observed richness") +
  ggtitle("Observed Richness") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

shannon_plot <- create_diversity_plot(richness_met, "Shannon", "Shannon diversity") +
  ggtitle("Shannon Diversity") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Combine plots using patchwork
combined_plot <- richness_plot + shannon_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Diversity Measures by Gut Region",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16)
    )
  )

# Display the combined plot
print(combined_plot)
View(richness_met)

# Print p-values
print("Wilcoxon test results for Observed Richness:")
print(wilcox.test(Observed ~ Anatomy, data = richness_met))

print("Wilcoxon test results for Shannon Diversity:")
print(wilcox.test(Shannon ~ Anatomy, data = richness_met))
####################################
# 11.3  Differential abundance analyses
####################################
library(phyloseq)
library(microbiomeMarker)
library(ggplot2)

phyloseq_samples_virome

perform_lefse_analysis <- function(phyloseq_obj, location) {
  print(paste("Analyzing location:", location))
  
  # Subset the data for the given location
  sample_names <- sample_names(phyloseq_obj)[sample_data(phyloseq_obj)$Location == location & 
                                               sample_data(phyloseq_obj)$Category != "diet"]
  phyloseq_location <- prune_samples(sample_names, phyloseq_obj)
  
  print(paste("Number of samples after subsetting:", nsamples(phyloseq_location)))
  
  # Rarefy to a depth of 10,000 reads per sample
  phyloseq_location_rarefied <- rarefy_even_depth(phyloseq_location, sample.size = 100000, rngseed = 123)
  
  print(paste("Number of samples after rarefaction:", nsamples(phyloseq_location_rarefied)))
  
  # Aggregate rare taxa at the Genus level
  phyloseq_location_rarefied_genus <- phyloseq_location_rarefied %>%
    aggregate_rare(level = "Species", detection = 0.0005/100, prevalence = 0.00001/100)
  
  print(paste("Number of samples after aggregation:", nsamples(phyloseq_location_rarefied_genus)))
  
  # Run LEfSe analysis
  lefse_result <- run_lefse(
    phyloseq_location_rarefied_genus,
    group = "Substrate_type",
    transform = "identity",
    taxa_rank = "Species",
    kw_cutoff = 0.05,
    wilcoxon_cutoff = 0.05,
    lda_cutoff = 3
  )
  
  # Get the marker table
  lefse_table <- marker_table(lefse_result)
  
  # Create the plot
  lefse_plot <- plot_ef_bar(lefse_table, label_level = 0) +
    scale_fill_manual(values = c("swill" = "#2a8a92", "chicken feed" = "#bd7f1c")) +
    ggtitle(paste("LEfSe Analysis -", location))
  
  # Return results
  return(list(table = lefse_table, plot = lefse_plot))
}

# Perform analysis for Italy and Belgium
set.seed(123)  # For reproducibility

italy_results <- perform_lefse_analysis(phyloseq_samples_virome, "Italy")
belgium_results <- perform_lefse_analysis(phyloseq_samples_virome, "Belgium")

# Print results
print("LEfSe results for Italy:")
print(italy_results$table)
print(italy_results$plot)

print("LEfSe results for Belgium:")
print(belgium_results$table)
print(belgium_results$plot)

Mastertable$names <- rownames(Mastertable)
Mastertable$Class[Mastertable$names =="NODE_B1_length_58442_cov_3283_912036_L_CF_2"]
Mastertable$Class[Mastertable$names =="NODE_A2_length_65724_cov_688_901580_G_SW_2A"]
Mastertable$Class[Mastertable$names =="NODE_B2_length_37010_cov_100_664338_L_CF_1"]
Mastertable$Class[Mastertable$names =="NODE_A2_length_34355_cov_1074_783301_G_CF_1C"]

Mastertable$Family[Mastertable$names=="NODE_B10_length_7424_cov_15461_806043_G_SW_3A"]
Mastertable$Class[Mastertable$names=="NODE_A3_length_50881_cov_2074_477994_G_SW_3C"]
Mastertable$Class[Mastertable$names=="NODE_B1_length_58581_cov_63_473609_G_SW_3A"]
Mastertable$Family[Mastertable$names=="NODE_A8_length_6848_cov_14972_125388_Swill6_6"]
Mastertable$Class[Mastertable$names=="NODE_A5_length_45221_cov_336_145512_G_SW_2B"]
Mastertable$Class[Mastertable$names=="NODE_A6_length_39704_cov_139_412724_G_SW_3A"]

Mastertable$names <- rownames(Mastertable)
####################################
# 12. Dietary Viruses and Their Impact on BSF Haemolymph Virome Composition
####################################
library(phyloseq)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# Add the new sample names to the subset
names1 <- c("Hae_CF", "Hae_SW", "Hae6","Swill6_4", "Swill6_5", "Swill6_6")
additional_samples <- c("CF6_6", "CF6_4", "CF6_5", "Swill6_4", "Swill6_6", "Swill6_5")
names1 <- c(names1, "Hae_CF", additional_samples)

# Ensure we have all samples
subset1 <- subset_samples(phyloseq_samples, Sample_identifier %in% names1)

# Extract the OTU table and convert to data frame
otu_table <- as.data.frame(as(otu_table(subset1), "matrix"))

print("Samples in otu_table:")
print(colnames(otu_table))

# If Hae_CF is not in the otu_table, add it with zero counts
if(!"Hae_CF" %in% colnames(otu_table)) {
  otu_table$Hae_CF <- 0
  print("Added Hae_CF column with zero counts")
}

# Check which Swill samples are present
swill_samples <- c("Swill6_4", "Swill6_5", "Swill6_6")
present_swill_samples <- intersect(swill_samples, colnames(otu_table))
print("Present Swill samples:")
print(present_swill_samples)

# Rarefy Swill samples to the same number of reads (e.g., 100,000 reads)
if(length(present_swill_samples) > 0) {
  swill_subset <- prune_samples(present_swill_samples, subset1)
  swill_subset_rarefied <- rarefy_even_depth(swill_subset, sample.size = 100000, rngseed = 123)
  swill_otu_table <- as.data.frame(as(otu_table(swill_subset_rarefied), "matrix"))
  
  # Merge rarefied OTU table back into the original OTU table
  for (sample in present_swill_samples) {
    otu_table[, sample] <- 0  # Initialize with zeros
    otu_table[rownames(swill_otu_table), sample] <- swill_otu_table[, sample]
  }
  print("Rarefied Swill samples to 100,000 reads")
}

# Combine available Swill samples into Diet (Swill)
if(length(present_swill_samples) > 0) {
  otu_table$Diet_Swill <- rowSums(otu_table[, present_swill_samples, drop = FALSE])
  print(paste("Combined", length(present_swill_samples), "Swill samples into Diet (Swill)"))
} else {
  otu_table$Diet_Swill <- 0
  print("No Swill samples found. Added empty Diet (Swill) column.")
}

# Add an empty Diet (CF) column
otu_table$Diet_CF <- 0
print("Added empty Diet (CF) column")

# Remove individual Swill6_* columns if they exist
keep_columns <- setdiff(colnames(otu_table), swill_samples)
otu_table <- otu_table[, keep_columns]

# Remove OTUs with zero counts in all samples
otu_table <- otu_table[rowSums(otu_table) > 0, ]

# Calculate relative abundance within each sample
otu_rel <- apply(otu_table, 2, function(x) x / sum(x))
# Replace NaN with 0 for the empty samples
otu_rel[is.nan(otu_rel)] <- 0

# Filter out OTUs with relative abundance below 1% in all samples
otu_rel <- otu_rel[rowSums(otu_rel >= 0.01) > 0, ]

# Get taxonomy table
tax_table <- as.data.frame(tax_table(subset1))

# Create a data frame with OTU abundances and their corresponding taxonomic information
otu_with_tax <- data.frame(
  OTU = rownames(otu_rel),
  Class = tax_table[rownames(otu_rel), "Class"],
  Family = tax_table[rownames(otu_rel), "Family"],
  otu_rel
)

# Order the data frame by Class
otu_with_tax <- otu_with_tax[order(otu_with_tax$Class), ]

# Separate the abundance data and taxonomic information
otu_matrix <- as.matrix(otu_with_tax[, -(1:3)])
rownames(otu_matrix) <- otu_with_tax$OTU

# Reorder columns in otu_matrix according to the specified order
desired_order <- c("Diet_Swill", "Hae_SW", "Diet_CF", "Hae_CF", "Hae6")
# Ensure only existing columns are included in the desired order
desired_order <- intersect(desired_order, colnames(otu_matrix))
otu_matrix <- otu_matrix[, desired_order]

# Create annotation for rows (OTUs)
row_annotation <- data.frame(
  Class = otu_with_tax$Class,
  Family = otu_with_tax$Family
)
rownames(row_annotation) <- rownames(otu_matrix)

# Update Substrate_type information
substrate_type <- sample_data(subset1)$Substrate_type
names(substrate_type) <- sample_names(subset1)

# Ensure Hae_CF has the correct Substrate_type if it wasn't in the original data
if("Hae_CF" %in% colnames(otu_matrix) && !("Hae_CF" %in% names(substrate_type))) {
  substrate_type["Hae_CF"] <- "chicken feed"
}

# Add substrate types for Diet columns
substrate_type["Diet_Swill"] <- "swill"
substrate_type["Diet_CF"] <- "chicken feed"

# Ensure substrate_type matches the order of columns in otu_matrix
substrate_type <- substrate_type[colnames(otu_matrix)]

# Create annotation for columns (samples)
col_annotation <- data.frame(Substrate_type = substrate_type)
rownames(col_annotation) <- colnames(otu_matrix)

# Get unique classes and families
unique_classes <- unique(otu_with_tax$Class)
unique_families <- unique(otu_with_tax$Family)

# Create color palettes for classes and families
n_classes <- length(unique_classes)
n_families <- length(unique_families)
class_colors <- colorRampPalette(c("#7AAB90", "grey88"))(n_classes)
family_colors <- colorRampPalette(c("#7AAB90", "grey88"))(n_families)
names(class_colors) <- unique_classes
names(family_colors) <- unique_families

# Create color palettes
color_palette <- c("white", colorRampPalette(brewer.pal(9, "YlOrRd"))(99))
substrate_colors <- c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")

# Create the heatmap
heatmap_plot <- pheatmap(otu_matrix,
         annotation_row = row_annotation,
         annotation_col = col_annotation,
         annotation_colors = list(
           Class = class_colors,
           Family = family_colors,
           Substrate_type = substrate_colors
         ),
         show_colnames = TRUE,
         show_rownames = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = color_palette,
         main = "Dietary Viruses and Their Impact on BSF Haemolymph Virome Composition",
         fontsize_row = 8,
         fontsize_col = 10,
         annotation_names_row = TRUE,
         annotation_legend = TRUE,
         scale = "none",
         angle_col = 45,  # Rotate column labels for better readability
         cellwidth = 60,
         cellheight = 50,
         annotation_legend_side = "right",
         legend_labels = c("0", "0.25", "0.5", "0.75", "1"),
         legend_breaks = c(0, 0.25, 0.5, 0.75, 1))

heatmap_plot

# Save the plot using ggsave
library(ggplot2)
ggsave("~/Desktop/viral_diet_haemolymp_heatmap.pdf", plot = heatmap_plot$gtable, width = 12, height = 10, device = "pdf")
####################################
# 13. Dietary Viruses and Their Impact on Small Larvae Virome Composition
####################################
# Add the new sample names to the subset
additional_samples <- c("CF6_6", "CF6_4", "CF6_5", "Swill6_4", "Swill6_6", "Swill6_5", "WL_CF_5", "WL_CF_Small", "WL_SW_5", "WL_SW_Small")
names1 <- c(names2, additional_samples)

# Ensure we have all samples
subset1 <- subset_samples(phyloseq_samples, Sample_identifier %in% names1)

# Extract the OTU table and convert to data frame
otu_table <- as.data.frame(as(otu_table(subset1), "matrix"))
print("Samples in otu_table:")
print(colnames(otu_table))

# Check which Swill samples are present
swill_samples <- c("Swill6_4", "Swill6_5", "Swill6_6")
present_swill_samples <- intersect(swill_samples, colnames(otu_table))
print("Present Swill samples:")
print(present_swill_samples)

# Rarefy Swill samples to the same number of reads (e.g., 100,000 reads)
if(length(present_swill_samples) > 0) {
  swill_subset <- prune_samples(present_swill_samples, subset1)
  swill_subset_rarefied <- rarefy_even_depth(swill_subset, sample.size = 100000, rngseed = 123)
  swill_otu_table <- as.data.frame(as(otu_table(swill_subset_rarefied), "matrix"))
  
  # Merge rarefied OTU table back into the original OTU table
  for (sample in present_swill_samples) {
    otu_table[, sample] <- 0  # Initialize with zeros
    otu_table[rownames(swill_otu_table), sample] <- swill_otu_table[, sample]
  }
  print("Rarefied Swill samples to 100,000 reads")
}

# Combine available Swill samples into Diet (Swill)
if(length(present_swill_samples) > 0) {
  otu_table$Diet_Swill <- rowSums(otu_table[, present_swill_samples, drop = FALSE])
  print(paste("Combined", length(present_swill_samples), "Swill samples into Diet (Swill)"))
} else {
  otu_table$Diet_Swill <- 0
  print("No Swill samples found. Added empty Diet (Swill) column.")
}

# Add an empty Diet (CF) column
otu_table$Diet_CF <- 0
print("Added empty Diet (CF) column")

# Remove individual Swill6_* columns if they exist
keep_columns <- setdiff(colnames(otu_table), swill_samples)
otu_table <- otu_table[, keep_columns]

# Remove OTUs with zero counts in all samples
otu_table <- otu_table[rowSums(otu_table) > 0, ]

# Calculate relative abundance within each sample
otu_rel <- apply(otu_table, 2, function(x) x / sum(x))
# Replace NaN with 0 for the empty samples
otu_rel[is.nan(otu_rel)] <- 0

# Filter out OTUs with relative abundance below 1% in all samples
otu_rel <- otu_rel[rowSums(otu_rel >= 0.01) > 0, ]

# Get taxonomy table
tax_table <- as.data.frame(tax_table(subset1))

# Create a data frame with OTU abundances and their corresponding taxonomic information
otu_with_tax <- data.frame(
  OTU = rownames(otu_rel),
  Class = tax_table[rownames(otu_rel), "Class"],
  Family = tax_table[rownames(otu_rel), "Family"],
  otu_rel
)

# Order the data frame by Class
otu_with_tax <- otu_with_tax[order(otu_with_tax$Class), ]

# Separate the abundance data and taxonomic information
otu_matrix <- as.matrix(otu_with_tax[, -(1:3)])
rownames(otu_matrix) <- otu_with_tax$OTU

# Reorder columns in otu_matrix according to the specified order
colnames(otu_matrix)
desired_order <- c("Diet_Swill", "WL_SW_5","WL_SW_Small", "Diet_CF", "WL_CF_5", "WL_CF_Small")

# Ensure only existing columns are included in the desired order
desired_order <- intersect(desired_order, colnames(otu_matrix))
otu_matrix <- otu_matrix[, desired_order]

# Create annotation for rows (OTUs)
row_annotation <- data.frame(
  Class = otu_with_tax$Class,
  Family = otu_with_tax$Family
)
rownames(row_annotation) <- rownames(otu_matrix)

# Update Substrate_type information
substrate_type <- sample_data(subset1)$Substrate_type
names(substrate_type) <- sample_names(subset1)

# Add substrate types for Diet columns
substrate_type["Diet_Swill"] <- "swill"
substrate_type["Diet_CF"] <- "chicken feed"

# Ensure substrate_type matches the order of columns in otu_matrix
substrate_type <- substrate_type[colnames(otu_matrix)]

# Create annotation for columns (samples)
col_annotation <- data.frame(Substrate_type = substrate_type)
rownames(col_annotation) <- colnames(otu_matrix)

# Get unique classes and families
unique_classes <- unique(otu_with_tax$Class)
unique_families <- unique(otu_with_tax$Family)

# Create color palettes for classes and families
n_classes <- length(unique_classes)
n_families <- length(unique_families)
class_colors <- colorRampPalette(c("#7AAB90", "grey88"))(n_classes)
family_colors <- colorRampPalette(c("#7AAB90", "grey88"))(n_families)
names(class_colors) <- unique_classes
names(family_colors) <- unique_families

# Create color palettes
color_palette <- c("white", colorRampPalette(brewer.pal(9, "YlOrRd"))(99))
substrate_colors <- c("chicken feed" = "#bd7f1c", "swill" = "#2a8a92")

# Create the heatmap
heatmap_plot <- pheatmap(otu_matrix,
                         annotation_row = row_annotation,
                         annotation_col = col_annotation,
                         annotation_colors = list(
                           Class = class_colors,
                           Family = family_colors,
                           Substrate_type = substrate_colors
                         ),
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         color = color_palette,
                         main = "Diet and its Impact on BSF Small Larvae Virome Composition",
                         fontsize_row = 8,
                         fontsize_col = 10,
                         annotation_names_row = TRUE,
                         annotation_legend = TRUE,
                         scale = "none",
                         angle_col = 45,  # Rotate column labels for better readability
                         cellwidth = 60,
                         cellheight = 50,
                         annotation_legend_side = "right",
                         legend_labels = c("0", "0.25", "0.5", "0.75", "1"),
                         legend_breaks = c(0, 0.25, 0.5, 0.75, 1))

heatmap_plot

# Save the plot using ggsave
library(ggplot2)
ggsave("~/Desktop/viral_diet_small_larvae_heatmap.pdf", plot = heatmap_plot$gtable, width = 12, height = 10, device = "pdf")

####################################
# 14 Complete virome BSFL figures
####################################
# 14.1 Barplot Life cycle
####################################
# Define columns to exclude
names_to_exclude <- c("Hae_CF", "Hae_SW", "Hae6", "Swill6_4", "Swill6_5", "Swill6_6",
                      "WL_CF_5", "WL_CF_Small", "WL_SW_5", "WL_SW_Small")

# Filter for high quality viral contigs
Mastertable_HQ_viral <- Mastertable_viral[Mastertable_viral$completeness > 50,]

# Remove excluded columns
Mastertable_HQ_viral <- Mastertable_HQ_viral[, !colnames(Mastertable_HQ_viral) %in% names_to_exclude]

# Identify sample columns (assuming they're the first 67 columns)
sample_cols <- colnames(Mastertable_HQ_viral)[1:57]

# Calculate read sums for each sample
read_sums <- colSums(Mastertable_HQ_viral[, sample_cols])

# Identify outliers (samples with less than 100,000 reads)
outliers <- names(read_sums[read_sums < 100000])

print("Samples removed as outliers:")
print(outliers)

# Remove outlier samples
Mastertable_HQ_viral <- Mastertable_HQ_viral[, !(colnames(Mastertable_HQ_viral) %in% outliers)]
# Recalculate sample columns after removal of outliers
sample_cols <- colnames(Mastertable_HQ_viral)[1:(57 - length(outliers))]

# Create abundance table with samples AND temperate column
abundance_table <- Mastertable_HQ_viral[, c(sample_cols, "temperate")]

# First aggregate the data by temperate/virulent status
lifecycle_abundances <- aggregate(abundance_table[, sample_cols], 
                                by = list(Lifecycle = abundance_table$temperate), 
                                FUN = sum)

rownames(lifecycle_abundances) <- lifecycle_abundances$Lifecycle
lifecycle_abundances$Lifecycle <- NULL
#View(lifecycle_abundances)

# RA calculation
lifecycle_abundances1 <- sweep(lifecycle_abundances, 2, colSums(lifecycle_abundances), '/')
lifecycle_abundances1[is.na(lifecycle_abundances1)] <- 0
lifecycle_abundances1

# Merge with metadata
lifecycle_abundances_meta <- merge(t(lifecycle_abundances1), samples, by = 0, all = F)
rownames(lifecycle_abundances_meta) <- lifecycle_abundances_meta$Row.names
lifecycle_abundances_meta$Row.names <- NULL
colnames(lifecycle_abundances_meta)

# Melt
lifecycle_abundances_m <- melt(lifecycle_abundances_meta, id.vars = c("Sample_identifier","Black_soldier_fly_identifier","Sample_type","Substrate_type", "Category", "Anatomy", "Anatomy2", "Location", "Larval_age_DAH"))

# Calculate sample sizes per group
n_per_group <- table(lifecycle_abundances_m$Substrate_type[lifecycle_abundances_m$variable == "Temperate"])
print(n_per_group)

median(lifecycle_abundances_m$value[lifecycle_abundances_m$variable=="Temperate" & lifecycle_abundances_m$Substrate_type=="chicken feed"])
median(lifecycle_abundances_m$value[lifecycle_abundances_m$variable=="Temperate" & lifecycle_abundances_m$Substrate_type=="swill"])

lifecycle_abundances_m$value
lifecycle_abundances_m$variable

# Perform Wilcoxon test
wilcox_result <- wilcox.test(value ~ Substrate_type, 
                             data = subset(lifecycle_abundances_m, variable == "Temperate"), 
                             exact = FALSE)
print(wilcox_result)

# Create plot with sample sizes in the title
ggplot(subset(lifecycle_abundances_m, variable == "Temperate"), 
       aes(x = Substrate_type, y = value, fill = Substrate_type)) + 
  geom_boxplot() + 
  scale_fill_hue(direction = 1) + 
  theme_bw() +
  labs(title = paste("Temperate Phage Abundance by Substrate Type\n",
                     "n =", paste(names(n_per_group), n_per_group, sep = ": ", collapse = ", ")),
       y = "Relative Abundance",
       x = "Substrate Type")
####################################
# 14.2 Host prediction phages
####################################
# Define columns to exclude
names_to_exclude <- c("Hae_CF", "Hae_SW", "Hae6", "Swill6_4", "Swill6_5", "Swill6_6",
                      "WL_CF_5", "WL_CF_Small", "WL_SW_5", "WL_SW_Small")

# Filter for high quality viral contigs
Mastertable_HQ_viral <- Mastertable_viral[Mastertable_viral$completeness > 50,]

# Remove excluded columns
Mastertable_HQ_viral <- Mastertable_HQ_viral[, !colnames(Mastertable_HQ_viral) %in% names_to_exclude]

# Identify sample columns (assuming they're the first 67 columns)
sample_cols <- colnames(Mastertable_HQ_viral)[1:57]

# Calculate read sums for each sample
read_sums <- colSums(Mastertable_HQ_viral[, sample_cols])

# Identify outliers (samples with less than 100,000 reads)
outliers <- names(read_sums[read_sums < 100000])

print("Samples removed as outliers:")
print(outliers)

# Remove outlier samples
Mastertable_HQ_viral <- Mastertable_HQ_viral[, !(colnames(Mastertable_HQ_viral) %in% outliers)]
# Recalculate sample columns after removal of outliers
sample_cols <- colnames(Mastertable_HQ_viral)[1:(57 - length(outliers))]

# Create abundance table with samples AND Host_phyla column
abundance_table <- Mastertable_HQ_viral[, c(sample_cols, "Host_phyla")]

# First aggregate the data by temperate/virulent status
lifecycle_abundances <- aggregate(abundance_table[, sample_cols], 
                                  by = list(Host = abundance_table$Host_phyla), 
                                  FUN = sum)

lifecycle_abundances$Host[lifecycle_abundances$Host=="0"] <- "Unclassified"
rownames(lifecycle_abundances) <- lifecycle_abundances$Host
lifecycle_abundances$Host <- NULL
lifecycle_abundances

# RA calculation
lifecycle_abundances1 <- sweep(lifecycle_abundances, 2, colSums(lifecycle_abundances), '/')
lifecycle_abundances1[is.na(lifecycle_abundances1)] <- 0
lifecycle_abundances1

# Merge with metadata
lifecycle_abundances_meta <- merge(t(lifecycle_abundances1), samples, by = 0, all = F)
rownames(lifecycle_abundances_meta) <- lifecycle_abundances_meta$Row.names
lifecycle_abundances_meta$Row.names <- NULL
colnames(lifecycle_abundances_meta)

# Melt
lifecycle_abundances_m <- melt(lifecycle_abundances_meta, id.vars = c("Sample_identifier","Black_soldier_fly_identifier","Sample_type","Substrate_type", "Category", "Anatomy", "Anatomy2", "Location", "Larval_age_DAH"))

# Function to run Wilcoxon test for each host phylum
run_stats <- function(data, phylum) {
  subset_data <- subset(data, variable == phylum)
  test <- wilcox.test(value ~ Substrate_type, data = subset_data, exact = FALSE)
  
  # Calculate means for each group
  means <- tapply(subset_data$value, subset_data$Substrate_type, mean)
  
  return(c(
    p.value = test$p.value,
    mean_chicken = means["chicken feed"],
    mean_swill = means["swill"],
    n_chicken = sum(subset_data$Substrate_type == "chicken feed"),
    n_swill = sum(subset_data$Substrate_type == "swill")
  ))
}

# Run tests for each phylum
phyla <- unique(lifecycle_abundances_m$variable)
stats_results <- lapply(phyla, run_stats, data = lifecycle_abundances_m)
names(stats_results) <- phyla

# Extract p-values and correct for multiple testing
p_values <- sapply(stats_results, function(x) x["p.value"])
p_adjusted <- p.adjust(p_values, method = "BH")  # Benjamini-Hochberg correction
p_values
p_adjusted

# Create plot
ggplot(lifecycle_abundances_m, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = Substrate_type)) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(
    y = "Relative Abundance",
    x = "Host Phylum",
    fill = "Substrate Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 12)
  ) +
  annotate(
    "text",
    x = seq_along(phyla),
    y = max(lifecycle_abundances_m$value) * 1.1,
    label = sprintf("p = %.3f\nadj.p = %.3f", p_values, p_adjusted)
  )
####################################
# 14.3 Sharing of viruses between swill-fed BSF
####################################
# Dataset (n=19 swill, n=45 total)
phyloseq_samples_virome <- phyloseq_virome
colnames(sample_data(phyloseq_samples_virome))

# DF
phyloseq_samples_virome_swill <- subset_samples(phyloseq_samples_virome, Substrate_type == "swill")
phyloseq_samples_virome_swill <- prune_samples(sample_sums(phyloseq_samples_virome_swill) > 0, phyloseq_samples_virome_swill)
RA_matrix_p <- as.data.frame(otu_table(phyloseq_samples_virome_swill))
presenceabsencep<-RA_matrix_p
presenceabsencep[presenceabsencep>0]<-1
presenceabsencep <- presenceabsencep[rowSums(presenceabsencep[, -1])>0, ]
ncol(presenceabsencep)

# Which Viruses are shared the most
presenceabsencep_Sharing <- presenceabsencep
presenceabsencep_Sharing$sharing <- rowSums(presenceabsencep_Sharing) # which viral 
presenceabsencep_Sharing <- presenceabsencep_Sharing[rev(order(presenceabsencep_Sharing$sharing)),]
rownames(presenceabsencep_Sharing)
presenceabsencep_Sharing$sharing
nrow(presenceabsencep_Sharing)

# Merge with metadata
presenceabsence_agg_metp <- merge(t(presenceabsencep), samples, by = 0, all =F)
rownames(presenceabsence_agg_metp) <-  presenceabsence_agg_metp$Row.names
presenceabsence_agg_metp$Row.names <- NULL

# Calculate shared phage families
pre <-aggregate(presenceabsence_agg_metp[, 1:nrow(presenceabsencep)], by=list(presenceabsence_agg_metp$Sample_identifier),FUN=max)
pre[,-1]<-lapply(pre[,-1],as.numeric)
rowSums(unique(pre[,-1]))
pre$totalcontigs <-rowSums(unique(pre[,-1]))
sum(pre$totalcontigs)

#find sharing of contigs between infants
sharedcontigs <-as.data.frame(colSums(pre[,-1]))
sharedcontigs <-  as.data.frame(t(sharedcontigs))
sharedcontigs$totalcontigs <- NULL
sharedcontigs <-  as.data.frame(t(sharedcontigs))
sharedcontigs_table <- as.data.frame(table(sharedcontigs$`colSums(pre[, -1])`))
colnames(sharedcontigs_table) <- c("Individuals", "Shared")
sharedcontigs_table
pre$totalcontig
pre$totalcontigs <- NULL
pre$Group.1 <- NULL
pre_DF <- as.data.frame(colSums(pre))
pre_DF$percentage <- (pre_DF$`colSums(pre)`/19)*100
nrow(pre_DF[pre_DF$`colSums(pre)`=="1",])

# Figure
ggplot(sharedcontigs_table) +
  aes(x = Individuals, y = Shared) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values="cadetblue") +
  xlab("Number of individuals") +
  ylab("Number of shared viruses") +
  ggtitle("Swill") +
  theme_classic()
####################################
# 14.4 Sharing of viruses between chicken-fed BSF
####################################
# Dataset (n=45 all; n=26 chicken feed)
phyloseq_samples_virome <- phyloseq_virome
colnames(sample_data(phyloseq_samples_virome))

# DF
phyloseq_samples_virome_CF <- subset_samples(phyloseq_samples_virome, Substrate_type == "chicken feed")
phyloseq_samples_virome_CF <- prune_samples(sample_sums(phyloseq_samples_virome_CF) > 0, phyloseq_samples_virome_CF)
RA_matrix_p <- as.data.frame(otu_table(phyloseq_samples_virome_CF))
presenceabsencep<-RA_matrix_p
presenceabsencep[presenceabsencep>0]<-1
presenceabsencep <- presenceabsencep[rowSums(presenceabsencep[, -1])>0, ]
nrow(presenceabsencep)

# Which Viruses are shared the most
presenceabsencep_Sharing <- presenceabsencep
presenceabsencep_Sharing$sharing <- rowSums(presenceabsencep_Sharing) # which viral 
presenceabsencep_Sharing <- presenceabsencep_Sharing[rev(order(presenceabsencep_Sharing$sharing)),]

# Merge with metadata
presenceabsence_agg_metp <- merge(t(presenceabsencep), samples, by = 0, all =F)
rownames(presenceabsence_agg_metp) <-  presenceabsence_agg_metp$Row.names
presenceabsence_agg_metp$Row.names <- NULL

# Calculate shared phage families
pre <- aggregate(presenceabsence_agg_metp[, 1:nrow(presenceabsencep)],by=list(presenceabsence_agg_metp$Sample_identifier),FUN=max)
pre[,-1] <- lapply(pre[,-1], as.numeric)
pre$totalcontigs <- rowSums(pre[,-1])
#View(pre)

# Merge with metadata
presenceabsence_agg_metp <- merge(t(presenceabsencep), samples, by = 0, all =F)
rownames(presenceabsence_agg_metp) <-  presenceabsence_agg_metp$Row.names
presenceabsence_agg_metp$Row.names <- NULL

#find sharing of contigs between infants
sharedcontigs <-as.data.frame(colSums(pre[,-1]))
sharedcontigs <-  as.data.frame(t(sharedcontigs))
sharedcontigs$totalcontigs <- NULL
sharedcontigs <-  as.data.frame(t(sharedcontigs))
sharedcontigs_table <- as.data.frame(table(sharedcontigs$`colSums(pre[, -1])`))
colnames(sharedcontigs_table) <- c("Individuals", "Shared")
sharedcontigs_table
pre$totalcontig
pre$totalcontigs <- NULL
pre$Group.1 <- NULL
pre_DF <- as.data.frame(colSums(pre))
pre_DF$percentage <- (pre_DF$`colSums(pre)`/26)*100
nrow(pre_DF[pre_DF$`colSums(pre)`=="1",])

# Figure
ggplot(sharedcontigs_table) +
  aes(x = Individuals, y = Shared) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values="cadetblue") +
  xlab("Number of individuals") +
  ylab("Number of shared viruses") +
  ggtitle("Chicken fed") +
  theme_classic()
####################################

######### !!!!!!!!!!!!!!! #########

    # SAVE Global Environment #
######### !!!!!!!!!!!!!!! #########

