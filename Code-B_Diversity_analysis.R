#####
# Supplementary code B
# Diversity analysis pipeline v1
# by Ostaizka Aizpurua and Antton Alberdi
# ostaizka.aizpurua@sund.ku.dk, antton.alberdi@sund.ku.dk
# January 2021
#####

#Load required libraries
library(reshape2)
library(gplots)
library(ggplot2)
library(ggrepel)
library(purrr)
library(hilldiv)
library(ape)
library(phytools)
library(phylobase)
library(phylosignal)
library(ggplot2)
library(vegan)
library(tidyverse)
library(phytools)
library(plyr)
library(dplyr)
library(dendextend)
library(ggpubr)
library(dmetar)
library(meta)
library(metafor)
library(TreeDist)

#Set working directory
setwd("~/github/Wild_Captive_16S")

###############
# INDEX
###############

# B1) PREPARE WORKING ABUNDANCE-TABLES
# B2) ACROSS-SPECIES DIFFERENCES
# B3) DIVERSITY SUMMARIES
# B4) OVERALL WILD vs CAPTIVE DIVERSITY META-ANALYSES
# B5) TAXON-SPECIFIC DIVERSITY META-ANALYSES
# B6) COMPOSITIONAL DIFFERENCES BETWEEN WILD AND CAPTIVE ANIMALS (DIVERSITY PARTITIONING)
# B7) PHYLOGENETIC SIGNAL
# B8) COMPOSITIONAL DIFFERENCES (ANALYSIS OF VARIANCE)
# B9) NULL MODELS
# B10) APROXIMATE BAYESIAN COMPUTATION
# B11) DISTRIBUTION OF THE ORIGIN OF THE DETECTED GENERA
# B12) PAIRWISE DISSIMILARITY CORRELATION BETWEEN WILD AND CAPTIVE
# B13) NEAREST DATASET
# B14) HIERARCHICAL CLUSTERING AND TOPOLOGICAL DIFFERENCES
# B15) TAXONOMIC ENRICHMENT ANALYSIS
# B16) NMDS plot

###############
# B1) PREPARE WORKING ABUNDANCE-TABLES
###############

#### Declare dataset (species) code names ####
list.dirs <- function(path="Species", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
                    full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}

code.list <- as.character(read.table("Data/sp_code.txt")[,1])

##
## B1.1) Generate genus-level read-abundance tables for each dataset
##

for (code in code.list){
  print(code)
  #Load individual genus-level abundance tables
  files.g = list.files(path=(paste("Species/",code,sep="")), pattern="*.tax_genus.txt")
  names.g = gsub(".tax_genus.txt","",files.g)
  filelist.g = lapply(paste("Species",code,files.g,sep="/"), function(x)read.table(x, header=T,sep="\t"))
  #Transform list to matrix
  filelist.g %>% purrr::reduce(full_join, by = "Taxa") -> count.table.g
  rownames(count.table.g) <- count.table.g[,1]
  count.table.g <- count.table.g[,-1]
  #Name columns
  colnames(count.table.g) <- names.g
  #Replace NAs by 0s
  count.table.g[is.na(count.table.g)] <- 0
  write.table(count.table.g, paste("Tables/count_genus_",code,".tsv",sep=""))
}

##
## B1.2) Filter genus-level read-abundance tables
##

# Load all the tables created
for (code in code.list){
  print(code)
  count.table <- read.table(paste("Tables/count_genus_",code,".tsv",sep=""))
  assign(paste("count.table.",code,sep=""),get("count.table"))
}

# Filter out samples with low sequencing depth (1000)
for (code in code.list){
  count.filtdepth <- depth_filt(get(paste("count.table.",code,sep="")),1000)
  count.filtdepth <- count.filtdepth[rowSums(count.filtdepth)>0, ]
  assign(paste("count.filtdepth.",code,sep=""),get("count.filtdepth"))
  write.table(get(paste("count.filtdepth.",code,sep="")),paste("Tables/countfiltdepth_",code,".tsv",sep=""))
}

# Filter out samples with low diversity coverage
for (code in code.list){
  filtcov <- as.data.frame(depth_cov(get(paste("count.filtdepth.",code,sep="")),qvalue=0))
  count.filtcov <- get(paste("count.filtdepth.",code,sep=""))
  count.filtcov <- count.filtcov[,rownames(filtcov[filtcov$Coverage > 99,])]
  assign(paste("count.filtcov.",code,sep=""),get("count.filtcov"))
  write.table(get(paste("count.filtcov.",code,sep="")),paste("Tables/countfiltcov_",code,".tsv",sep=""))
}

# Remove non-bacterial or taxonomy-lacking genera from the dataset and simplify names
taxonlist_NA <- read.table("Data/NA_taxa.txt")
taxonlist_NA <-  as.character(taxonlist_NA[,1])

for (code in code.list){
  count.filtdepth <- read.table(paste("Tables/countfiltcov_",code,".tsv",sep=""))
  count.filtered <- count.filtdepth[! rownames(count.filtdepth) %in% taxonlist_NA,]
  rownames(count.filtered) <- gsub(" ","_",rownames(count.filtered))
  rownames(count.filtered) <- gsub("[","",rownames(count.filtered),fixed = TRUE)
  rownames(count.filtered) <- gsub("]","",rownames(count.filtered),fixed = TRUE)
  rownames(count.filtered) <- gsub("(","-",rownames(count.filtered),fixed = TRUE)
  rownames(count.filtered) <- gsub(")","-",rownames(count.filtered),fixed = TRUE)
  assign(paste("count.filtered.",code,sep=""),get("count.filtered"))
  write.table(get(paste("count.filtered.",code,sep="")),paste("Tables/countfiltered_",code,".tsv",sep=""))
}

# Calculate read percentage of unclassified

for (code in code.list){
  count.filtdepth <- read.table(paste("Tables/countfiltcov_",code,".tsv",sep=""))
  count.filtdepth.sums <- colSums(count.filtdepth)
  count.filtered <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  count.filtered.sums <- colSums(count.filtered)
  unclassified <- mean(100 - (count.filtered.sums / count.filtdepth.sums * 100))
  print(code)
  print(unclassified)
}

##
## B1.3) Merge all species-specific tables into a single table
##

#Load filtered genus tables
files.g = list.files(path="Tables",pattern="countfiltered_*")
filelist.g = lapply(paste("Tables",files.g,sep="/"), function(x)read.table(x, header=T))
rownameslist = lapply(filelist.g, function(x) rownames(x))
for (l in c(1:length(filelist.g))){
  filelist.g[[l]]$taxa <- rownames(filelist.g[[l]])
}

#Transform list to matrix
filelist.g %>% purrr::reduce(full_join, by = "taxa") -> count.table.all.g
rownames(count.table.all.g) <- count.table.all.g[,"taxa"]
count.table.all.g <- count.table.all.g[,colnames(count.table.all.g) != "taxa"]

#Replace NAs by 0s
count.table.all.g[is.na(count.table.all.g)] <- 0
write.table(count.table.all.g, "Tables/count_Genus_all.tsv")

#Filter metadata
metadata <- read.table("Data/metadata.csv", sep =";",row.names=1)
metadata.filtered <- metadata[which(metadata[,"Sample"] %in% colnames(count.table.all.g)),]
write.table(metadata.filtered, "Data/metadata.filtered.csv",sep=";")

###############
# B2) ACROSS-SPECIES DIFFERENCES
###############

# Load required data
count.table.all.g <- read.table("Tables/count_Genus_all.tsv")
metadata.filtered <- read.table("Data/metadata.filtered.csv", sep =";",row.names=1)
code.list <- as.character(read.table("Data/sp_code.txt")[,1])

##
## B2.1) Alpha diversity differences across species
##

#dR
div_dR.all <- div_test(count.table.all.g,qvalue=0,hierarchy=metadata.filtered[,c("Sample","Species")])
saveRDS(div_dR.all,"Results/RDS/div_dR.all.RData")

#dRE
div_dRE.all <- div_test(count.table.all.g,qvalue=1,hierarchy=metadata.filtered[,c("Sample","Species")])
saveRDS(div_dRE.all,"Results/RDS/div_dRE.all.RData")

#dRER
capwild.tree <- read.tree("Data/genustree.tre")
tree_filtered <- match_data(count.table.all.g,capwild.tree,output="tree")
div_dRER.all <- div_test(count.table.all.g,qvalue=1,hierarchy=metadata.filtered[,c("Sample","Species")],tree = tree_filtered)
saveRDS(div_dRER.all,"Results/RDS/div_dRER.all.RData")

##
## B2.2) Compositional differences across species
##

# dR
pairdis_dR.all <- pair_dis(count.table.all.g,qvalue=0,hierarchy=metadata.filtered[,c("Sample","Species")])
saveRDS(pairdis_dR.all,"Results/RDS/pairdis_dR.all.RData")
pairdis_dR.all <- readRDS("Results/RDS/pairdis_dR.all.RData")
dR.data <- pairdis_dR.all$L1_UqN
dR.dist <- as.dist(dR.data)
ps.disper.dR.species <- betadisper(dR.dist, metadata.filtered$Species)
permutest(ps.disper.dR.species, pairwise = TRUE)
adonis(dR.dist ~ Species, data =metadata.filtered, permutations = 999)

# dRE
pairdis_dRE.all <- pair_dis(count.table.all.g,qvalue=1,hierarchy=metadata.filtered[,c("Sample","Species")])
saveRDS(pairdis_dRE.all,"Results/RDS/pairdis_dRE.all.RData")
pairdis_dRE.all <- readRDS("Results/RDS/pairdis_dRE.all.RData")
dRE.data <- pairdis_dRE.all$L1_UqN
dRE.dist <- as.dist(dRE.data)
ps.disper.dRE.species <- betadisper(dRE.dist, metadata.filtered$Species)
permutest(ps.disper.dRE.species, pairwise = TRUE)
adonis(dRE.dist ~ Species, data = metadata.filtered, permutations = 999)

# dRER
capwild.tree <- read.tree("Data/genustree.tre")
tree_filtered <- match_data(count.table.all.g,capwild.tree,output="tree")
 #Note the next step may take multiple days to be completed
pairdis_dRER.all <- pair_dis(count.table.all.g,qvalue=1,hierarchy=metadata.filtered[,c("Sample","Species")], tree = tree_filtered)
saveRDS(pairdis_dRER.all,"Results/RDS/pairdis_dRER.all.RData")
pairdis_dRER.all <- readRDS("Results/RDS/pairdis_dRER.all.RData")
dRER.data <- pairdis_dRER.all$L1_UqN
dRER.dist <- as.dist(dRER.data)
ps.disper.dRER.species <- betadisper(dRER.dist, metadata.filtered$Species)
permutest(ps.disper.dRER.species, pairwise = TRUE)
adonis(dRER.dist ~ Species, data = metadata.filtered, permutations = 999)

###############
# B3) DIVERSITY SUMMARIES
###############

# Load required data
count.table.all.g <- read.table("Tables/count_Genus_all.tsv")
metadata.filtered <- read.table("Data/metadata.filtered.csv", sep =";",row.names=1)
code.list <- as.character(read.table("Data/sp_code.txt")[,1])
capwild.tree <- read.tree("Data/genustree.tre")

##
## B3.1) Summary of diversity values based on richness (R)
##

summary_dR <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  metadata.filtered.subset <- metadata.filtered[which(metadata.filtered[,1] %in% colnames(final.table)),]
  divqR <- div_test(final.table,qvalue=0,hierarchy=metadata.filtered.subset[,c("Sample","Origin")])
  divqR_table <- divqR$data
  n <- nrow(divqR_table)
  mean <- mean(divqR_table[,c("Value")])
  sd <- sd(divqR_table[,c("Value")])
  divqR_wild <- divqR_table[divqR_table$Group == "Wild",]
  n_wild <- nrow(divqR_wild)
  mean_wild <- mean(divqR_wild[,c("Value")])
  sd_wild <- sd(divqR_wild[,c("Value")])
  divqR_captive <- divqR_table[divqR_table$Group == "Captivity",]
  n_captive <- nrow(divqR_captive)
  mean_captive <- mean(divqR_captive[,c("Value")])
  sd_captive <- sd(divqR_captive[,c("Value")])
  summary_dR <- rbind(summary_dR,c(n,round(mean,3),round(sd,3),round(n_wild,3),round(mean_wild,3),round(sd_wild,3),round(n_captive,3),round(mean_captive,3),round(sd_captive,3)))
}

colnames(summary_dR) <- c("n","mean","sd","n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_dR) <- code.list
write.table(summary_dR, "Results/summary_diversity_dR.tsv")

##
## B3.1) Summary of diversity values based on richness (R)
##

summary_dRE <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  metadata.filtered.subset <- metadata.filtered[which(metadata.filtered[,1] %in% colnames(final.table)),]
  divqRE <- div_test(final.table,qvalue=1,hierarchy=metadata.filtered.subset[,c("Sample","Origin")])
  divqRE_table <- divqRE$data
  n <- nrow(divqRE_table)
  mean <- mean(divqRE_table[,c("Value")])
  sd <- sd(divqRE_table[,c("Value")])
  divqRE_wild <- divqRE_table[divqRE_table$Group == "Wild",]
  n_wild <- nrow(divqRE_wild)
  mean_wild <- mean(divqRE_wild[,c("Value")])
  sd_wild <- sd(divqRE_wild[,c("Value")])
  divqRE_captive <- divqRE_table[divqRE_table$Group == "Captivity",]
  n_captive <- nrow(divqRE_captive)
  mean_captive <- mean(divqRE_captive[,c("Value")])
  sd_captive <- sd(divqRE_captive[,c("Value")])
  summary_dRE <- rbind(summary_dRE,c(n,round(mean,3),round(sd,3),round(n_wild,3),round(mean_wild,3),round(sd_wild,3),round(n_captive,3),round(mean_captive,3),round(sd_captive,3)))
}

colnames(summary_dRE) <- c("n","mean","sd","n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_dRE) <- code.list
write.table(summary_dRE, "Results/summary_diversity_dRE.tsv")

##
## B3.3) Summary of diversity values based on richness+eveness+regularity (REH)
##

summary_dRER <- c()
for (code in code.list){
  print(code)
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  metadata.filtered.subset <- metadata.filtered[which(metadata.filtered[,1] %in% colnames(final.table)),]
  tree_filtered <- match_data(final.table,capwild.tree,output="tree")
  divqRER <- div_test(final.table,qvalue=1,hierarchy=metadata.filtered.subset[,c("Sample","Origin")], tree=tree_filtered)
  divqRER_table <- divqRER$data
  n <- nrow(divqRER_table)
  mean <- mean(divqRER_table[,c("Value")])
  sd <- sd(divqRER_table[,c("Value")])
  divqRER_wild <- divqRER_table[divqRER_table$Group == "Wild",]
  n_wild <- nrow(divqRER_wild)
  mean_wild <- mean(divqRER_wild[,c("Value")])
  sd_wild <- sd(divqRER_wild[,c("Value")])
  divqRER_captive <- divqRER_table[divqRER_table$Group == "Captivity",]
  n_captive <- nrow(divqRER_captive)
  mean_captive <- mean(divqRER_captive[,c("Value")])
  sd_captive <- sd(divqRER_captive[,c("Value")])
  summary_dRER <- rbind(summary_dRER,c(n,round(mean,3),round(sd,3),round(n_wild,3),round(mean_wild,3),round(sd_wild,3),round(n_captive,3),round(mean_captive,3),round(sd_captive,3)))
}
colnames(summary_dRER) <- c("n","mean","sd","n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_dRER) <- code.list
write.table(summary_dRER, "Results/summary_diversity_dRER.tsv")

###############
# B4) OVERALL WILD vs CAPTIVE DIVERSITY META-ANALYSES
###############

##
## B4.1) Meta-analysis based on richness (dR)
##

summary_dR <- read.table("Results/summary_diversity_dR.tsv")
summary_dR <- as.data.frame(summary_dR)
meta_dR_ready <- tibble::rownames_to_column(summary_dR,"Author")
rownames(meta_dR_ready) <- meta_dR_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_dR_ready <- meta_dR_ready[sp_sorted,]
meta_dR.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                        data = meta_dR_ready,
                        studlab = paste(Author),
                        comb.fixed = FALSE,
                        comb.random = TRUE,
                        method.tau = "SJ",
                        hakn = TRUE,
                        prediction = TRUE,
                        sm = "SMD")
saveRDS(meta_dR.raw, "Results/RDS/meta_dR.all.RData")

#Forest plot dR
meta_dR.raw <- readRDS("Results/RDS/meta_dR.all.RData")
pdf("Results/Plots/forest_dR.pdf", width=13, height=8)
forest(meta_dR.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect",
       leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),lab.e = "Captivity",lab.c="Wild")
dev.off()

#Finding outliers
meta_dR.raw <- readRDS("Results/RDS/meta_dR.all.RData")
find.outliers(meta_dR.raw)

#Sensitivity analysis
summary_dR <- read.table("Results/summary_diversity_dR.tsv")
summary_dR_outlier <- summary_dR[-c(9,17,24),]
summary_dR_outlier_ready <- tibble::rownames_to_column(summary_dR_outlier,"Author")
rownames(summary_dR_outlier_ready) <- summary_dR_outlier_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","PELE","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
summary_dR_outlier_ready <- summary_dR_outlier_ready[sp_sorted,]
summary_dR_outlier_ready.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                                      data = summary_dR_outlier_ready,
                                      studlab = paste(Author),
                                      comb.fixed = FALSE,
                                      comb.random = TRUE,
                                      method.tau = "SJ",
                                      hakn = TRUE,
                                      prediction = TRUE,
                                      sm = "SMD")

#GOSH Plot Analysis
meta_dR.raw <- readRDS("Results/RDS/meta_dR.all.RData")
m.rma.dR <- rma(yi = meta_dR.raw$TE,
             sei = meta_dR.raw$seTE,
             method = meta_dR.raw$method.tau,
             test = "knha")
#GOSH plot dR
dat.gosh.dR <- gosh(m.rma.dR)
plot(dat.gosh.dR, alpha= 0.1, col = "blue")

#Clusters detection in the GOSH plot data
gosh.diagnostics(dat.gosh.dR)
#Re-run the meta-analysis without the species identified in the GLOSH plot analysis
summary_dR <- read.table("Results/summary_diversity_dR.tsv")
summary_dR <- as.data.frame(summary_dR)
meta_dR_ready <- tibble::rownames_to_column(summary_dR,"Author")
rownames(meta_dR_ready) <- meta_dR_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_dR_ready <- meta_dR_ready[sp_sorted,]
meta_dR.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                       data = meta_dR_ready,
                       studlab = paste(Author),
                       comb.fixed = FALSE,
                       comb.random = TRUE,
                       method.tau = "SJ",
                       hakn = TRUE,
                       prediction = TRUE,
                       sm = "SMD",
                       exclude = c(3,11,12,15,22))#Species to be excluded from the meta-analysis

##
## B4.2) Meta-analysis based on richness+eveness (dRE)
##

summary_dRE <- read.table("Results/summary_diversity_dRE.tsv")
summary_dRE <- as.data.frame(summary_dRE)
meta_dRE_ready <- tibble::rownames_to_column(summary_dRE,"Author")
rownames(meta_dRE_ready) <- meta_dRE_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_dRE_ready <- meta_dRE_ready[sp_sorted,]
meta_dRE.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                      data = meta_dRE_ready,
                      studlab = paste(Author),
                      comb.fixed = FALSE,
                      comb.random = TRUE,
                      method.tau = "SJ",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD")
saveRDS(meta_dRE.raw, "Results/RDS/meta_dRE.all.RData")

#Forest plot dR
meta_dRE.raw <- readRDS("Results/RDS/meta_dRE.all.RData")
pdf("Results/Plots/forest_dRE.pdf", width=13, height=8)
forest(meta_dRE.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect", leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),lab.e = "Captivity",lab.c="Wild")
dev.off()

#Finding outliers
meta_dRE.raw <- readRDS("Results/RDS/meta_dRE.all.RData")
find.outliers(meta_dRE.raw)

#Sensitivity analysis
summary_dRE <- read.table("Results/summary_diversity_dRE.tsv")
summary_dRE_outlier <- summary_dRE[-c(9,17,24),]
summary_dRE_outlier_ready <- tibble::rownames_to_column(summary_dRE_outlier,"Author")
rownames(summary_dRE_outlier_ready) <- summary_dRE_outlier_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","PELE","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
summary_dRE_outlier_ready <- summary_dRE_outlier_ready[sp_sorted,]
summary_dRE_outlier_ready.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                                            data = summary_dRE_outlier_ready,
                                            studlab = paste(Author),
                                            comb.fixed = FALSE,
                                            comb.random = TRUE,
                                            method.tau = "SJ",
                                            hakn = TRUE,
                                            prediction = TRUE,
                                            sm = "SMD")

#GOSH Plot Analysis
meta_dRE.raw <- readRDS("Results/RDS/meta_dRE.all.RData")
m.rma.dRE <- rma(yi = meta_dRE.raw$TE,
             sei = meta_dRE.raw$seTE,
             method = meta_dRE.raw$method.tau,
             test = "knha")

#GOSH plot dR
dat.gosh.dRE <- gosh(m.rma.dRE)
plot(dat.gosh.dRE, alpha= 0.1, col = "blue")

#Clusters detection in the GOSH plot data
gosh.diagnostics(dat.gosh.dRE)
#Re-run the meta-analysis without the species identified in the GLOSH plot analysis
summary_dRE <- read.table("Results/summary_diversity_dRE.tsv")
summary_dRE <- as.data.frame(summary_dRE)
meta_dRE_ready <- tibble::rownames_to_column(summary_dRE,"Author")
rownames(meta_dRE_ready) <- meta_dRE_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_dRE_ready <- meta_dRE_ready[sp_sorted,]
meta_dRE.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                         data = meta_dRE_ready,
                         studlab = paste(Author),
                         comb.fixed = FALSE,
                          comb.random = TRUE,
                          method.tau = "SJ",
                          hakn = TRUE,
                          prediction = TRUE,
                          sm = "SMD",
                          exclude = c(3,11,12,15,22))#Species to be excluded from the meta-analysis


##
## B4.3) Meta-analysis based on richness+eveness+regularity (dRER)
##

summary_dRER <- read.table("Results/summary_diversity_dRER.tsv")
summary_dRER <- as.data.frame(summary_dRER)
meta_dRER_ready <- tibble::rownames_to_column(summary_dRER,"Author")
rownames(meta_RER_ready) <- meta_RER_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_dRER_ready <- meta_dRER_ready[sp_sorted,]
meta_dRER.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                             data = meta_dRER_ready,
                             studlab = paste(Author),
                             comb.fixed = FALSE,
                             comb.random = TRUE,
                             method.tau = "SJ",
                             hakn = TRUE,
                             prediction = TRUE,
                             sm = "SMD")
saveRDS(meta_dRER.raw, "Results/RDS/meta_dRER.all.RData")

#Forest plot dRER
meta_dRER.raw <- readRDS("Results/RDS/meta_dRER.all.RData")
pdf("Results/Plots/forest_dRER.pdf",width=13, height=8)
forest(meta_dRER.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect",
       leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),lab.e = "Captivity",lab.c="Wild")
dev.off()

#Finding outliers
meta_dRER.raw <- readRDS("Results/RDS/meta_dRER.all.RData")
find.outliers(meta_dRER.raw)

#Sensitivity analysis
summary_dRER <- read.table("Results/summary_diversity_dRER.tsv")
summary_dRER_outlier <- summary_dRER[-c(1,5,11,17,18,20,24),]
summary_dRER_outlier_ready <- tibble::rownames_to_column(summary_dRER_outlier,"Author")
rownames(summary_dRER_outlier_ready) <- summary_dRER_outlier_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","PAAN","PATR","PELE","BOGA","ELDA","EQKI","PATI","MYTR","SAHA1","SAHA2","LALT")
summary_dRER_outlier_ready.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                                      data = summary_dRER_outlier_ready,
                                      studlab = paste(Author),
                                      comb.fixed = FALSE,
                                      comb.random = TRUE,
                                      method.tau = "SJ",
                                      hakn = TRUE,
                                      prediction = TRUE,
                                      sm = "SMD")

#GOSH Plot Analysis
meta_dRER.raw <- readRDS("Results/RDS/meta_dRER.all.RData")
m.rma.dRER <- rma(yi = meta_dRER.raw$TE,
             sei = meta_dRER.raw$seTE,
             method = meta_dRER.raw$method.tau,
             test = "knha")
#GOSH plot dR
dat.gosh.dRER <- gosh(m.rma.dRER)
plot(dat.gosh.dRER, alpha= 0.1, col = "blue")

#Clusters detection in the GOSH plot data
gosh.diagnostics(dat.gosh.dRER)
#Re-run the meta-analysis without the species identified in the GLOSH plot analysis
summary_dRER <- read.table("Results/summary_diversity_dRER.tsv")
summary_dRER <- as.data.frame(summary_dRER)
meta_dRER_ready <- tibble::rownames_to_column(summary_dRER,"Author")
rownames(meta_dRER_ready) <- meta_dRER_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_dRER.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                       data = meta_dRER_ready,
                       studlab = paste(Author),
                       comb.fixed = FALSE,
                       comb.random = TRUE,
                       method.tau = "SJ",
                       hakn = TRUE,
                       prediction = TRUE,
                       sm = "SMD",
                       exclude = c(4,5,8,9,14,15,22,24))#Species to be excluded from the meta-analysis


###############
# B5) TAXON-SPECIFIC DIVERSITY META-ANALYSES
###############

##
## B5.1) Meta-analysis of primates based on richness+eveness+regularity (dR)
##

#Subset diversity table
summary_dR <- read.table("Results/summary_diversity_dR.tsv")
sublist1 <- c("RHBR","PYNE","PAAN","PATR","GOGO")
summary_dR.subset1 <- summary_dR[rownames(summary_dR) %in% sublist1,]

#Run meta-analysis
meta_dR_sub1_ready <- tibble::rownames_to_column(summary_dR.subset1,"Author")
rownames(meta_dR_sub1_ready) <- meta_dR_sub1_ready$Author
meta_dR_sub1_ready <- meta_dR_sub1_ready[sublist1,]
meta_dR_sub1.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                            data = meta_dR_sub1_ready,
                            studlab = paste(Author),
                            comb.fixed = FALSE,
                            comb.random = TRUE,
                            method.tau = "SJ",
                            hakn = TRUE,
                            prediction = TRUE,
                            sm = "SMD")
saveRDS(meta_dR_sub1.raw,"Results/RDS/meta_dR.primates.RData")

##
## B5.2) Meta-analysis of primates based on richness+eveness (dRE)
##

#Subset diversity table
summary_dRER <- read.table("Results/summary_diversity_dRE.tsv")
sublist1 <- c("RHBR","PYNE","PAAN","PATR","GOGO")
summary_dRER.subset1 <- summary_dRER[rownames(summary_dRER) %in% sublist1,]

#Run meta-analysis
meta_dRER_sub1_ready <- tibble::rownames_to_column(summary_dRER.subset1,"Author")
rownames(meta_dRER_sub1_ready) <- meta_dRER_sub1_ready$Author
meta_dRER_sub1_ready <- meta_dRER_sub1_ready[sublist1,]
meta_dRER_sub1.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                                 data = meta_dRER_sub1_ready,
                                 studlab = paste(Author),
                                 comb.fixed = FALSE,
                                 comb.random = TRUE,
                                 method.tau = "SJ",
                                 hakn = TRUE,
                                 prediction = TRUE,
                                 sm = "SMD")

saveRDS(meta_dRER_sub1.raw,"Results/RDS/meta_dRE.primates.RData")

##
## B5.3) Meta-analysis of primates based on richness+eveness+regularity (dRER)
##

#Subset diversity table
summary_dRER <- read.table("Results/summary_diversity_dRER.tsv")
sublist1 <- c("RHBR","PYNE","PAAN","PATR","GOGO")
summary_dRER.subset1 <- summary_dRER[rownames(summary_dRER) %in% sublist1,]

#Run meta-analysis
meta_dRER_sub1_ready <- tibble::rownames_to_column(summary_dRER.subset1,"Author")
rownames(meta_dRER_sub1_ready) <- meta_dRER_sub1_ready$Author
meta_dRER_sub1_ready <- meta_dRER_sub1_ready[sublist1,]
meta_dRER_sub1.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                                 data = meta_dRER_sub1_ready,
                                 studlab = paste(Author),
                                 comb.fixed = FALSE,
                                 comb.random = TRUE,
                                 method.tau = "SJ",
                                 hakn = TRUE,
                                 prediction = TRUE,
                                 sm = "SMD")

saveRDS(meta_dRER_sub1.raw,"Results/RDS/meta_dRER.primates.RData")

##
## B5.4) Meta-analysis of cetartiodactylans based on richness (dR)
##

#Subset diversity table
summary_dR <- read.table("Results/summary_diversity_dR.tsv")
sublist2 <- c("TUTR","MOCH","BOGA","ELDA","CENI")
summary_dR.subset2 <- summary_dR[rownames(summary_dR) %in% sublist2,]

#Run meta-analysis
meta_dR_sub2_ready <- tibble::rownames_to_column(summary_dR.subset2,"Author")
rownames(meta_dR_sub2_ready) <- meta_dR_sub2_ready$Author
meta_dR_sub2_ready <- meta_dR_sub2_ready[sublist2,]
meta_dR_sub2.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                            data = meta_dR_sub2_ready,
                            studlab = paste(Author),
                            comb.fixed = FALSE,
                            comb.random = TRUE,
                            method.tau = "SJ",
                            hakn = TRUE,
                            prediction = TRUE,
                            sm = "SMD")
saveRDS(meta_dR_sub2.raw,"Results/RDS/meta_dR.cetartiodactyla.RData")

##
## B5.5) Meta-analysis of cetartiodactylans based on richness+eveness (dRE)
##

#Subset diversity table
summary_dRER <- read.table("Results/summary_diversity_dRE.tsv")
sublist2 <- c("TUTR","MOCH","BOGA","ELDA","CENI")
summary_dRER.subset2 <- summary_dRER[rownames(summary_dRER) %in% sublist2,]

#Run meta-analysis
meta_RER_sub2_ready <- tibble::rownames_to_column(summary_dRER.subset2,"Author")
rownames(meta_RER_sub2_ready) <- meta_RER_sub2_ready$Author
meta_RER_sub2_ready <- meta_RER_sub2_ready[sublist2,]
meta_RER_sub2.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                                 data = meta_RER_sub2_ready,
                                 studlab = paste(Author),
                                 comb.fixed = FALSE,
                                 comb.random = TRUE,
                                 method.tau = "SJ",
                                 hakn = TRUE,
                                 prediction = TRUE,
                                 sm = "SMD")

saveRDS(meta_RER_sub2.raw,"Results/RDS/meta_dRE.cetartiodactyla.RData")


##
## B5.6) Meta-analysis of cetartiodactylans based on richness+eveness+regularity (dRER)
##

#Subset diversity table
summary_dRER <- read.table("Results/summary_diversity_dRER.tsv")
sublist2 <- c("TUTR","MOCH","BOGA","ELDA","CENI")
summary_dRER.subset2 <- summary_dRER[rownames(summary_dRER) %in% sublist2,]

#Run meta-analysis
meta_RER_sub2_ready <- tibble::rownames_to_column(summary_dRER.subset2,"Author")
rownames(meta_RER_sub2_ready) <- meta_RER_sub2_ready$Author
meta_RER_sub2_ready <- meta_RER_sub2_ready[sublist2,]
meta_RER_sub2.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                                 data = meta_RER_sub2_ready,
                                 studlab = paste(Author),
                                 comb.fixed = FALSE,
                                 comb.random = TRUE,
                                 method.tau = "SJ",
                                 hakn = TRUE,
                                 prediction = TRUE,
                                 sm = "SMD")

saveRDS(meta_RER_sub2.raw,"Results/RDS/meta_dRER.cetartiodactyla.RData")

###############
# B6) COMPOSITIONAL DIFFERENCES BETWEEN WILD AND CAPTIVE ANIMALS (DIVERSITY PARTITIONING)
###############

# Load required data
metadata.filtered <- read.table("Data/metadata.filtered.csv", sep =";",row.names=1)
code.list <- as.character(read.table("Data/sp_code.txt")[,1])
capwild.tree <- read.tree("Data/genustree.tre")

##
## B6.1) Diversity partitioning based on richness (dR)
##

betadis_dR_results <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  ##Filter the metadata
  metadata.filtered.subset <- metadata.filtered[which(metadata.filtered[,1] %in% colnames(final.table)),]
  #Check whether both lists are identical
  identical(sort(colnames(final.table)),sort(as.character(metadata.filtered.subset[,"Sample"])))
  divpart.dR <- div_part(final.table,qvalue=0,hierarchy=metadata.filtered.subset[,c("Sample","Origin")])
  betadis.dR <- beta_dis(divpart.dR)
  betadis_dR_results <- append(betadis_dR_results,betadis.dR$UqN[2])
}
names(betadis_dR_results) <- code.list

#Statistics
mean(betadis_dR_results)
sd(betadis_dR_results)
max(betadis_dR_results)
min(betadis_dR_results)

##
## B6.2) Diversity partitioning based on richness+eveness (dRE)
##

betadis_dRE_results <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  ##Filter the metadata
  metadata.filtered.subset <- metadata.filtered[which(metadata.filtered[,1] %in% colnames(final.table)),]
  #Check whether both lists are identical
  identical(sort(colnames(final.table)),sort(as.character(metadata.filtered.subset[,"Sample"])))
  divpart.dRE <- div_part(final.table,qvalue=1,hierarchy=metadata.filtered.subset[,c("Sample","Origin")])
  betadis.dRE <- beta_dis(divpart.dRE)
  betadis_dRE_results <- append(betadis_dRE_results,betadis.dRE$UqN[2])
}
names(betadis_dRE_results) <- code.list

#Statistics
mean(betadis_dRE_results)
sd(betadis_dRE_results)
max(betadis_dRE_results)
min(betadis_dRE_results)

##
## B6.3) Diversity partitioning based on richness+eveness+regularity (dRER)
##

betadis_dRER_results <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  ##Filter the metadata
  metadata.filtered.subset <- metadata.filtered[which(metadata.filtered[,1] %in% colnames(final.table)),]
  #Check whether both lists are identical
  identical(sort(colnames(final.table)),sort(as.character(metadata.filtered.subset[,"Sample"])))
  tree_filtered <- match_data(final.table,capwild.tree,output="tree")
  divpart.dRER <- div_part(final.table,qvalue=1,hierarchy=metadata.filtered.subset[,c("Sample","Origin")],tree=tree_filtered)
  betadis.dRER <- beta_dis(divpart.dRER)
  betadis_dRER_results <- append(betadis_dRER_results,betadis.dRER$UqN[2])
}
names(betadis_dRER_results) <- code.list

#Statistics
mean(betadis_dRER_results)
sd(betadis_dRER_results)
max(betadis_dRER_results)
min(betadis_dRER_results)

#Join results
betadis_results <- cbind(betadis_dR_results,betadis_dRE_results,betadis_dRER_results)
write.table(betadis_results,"Results/summary_betadis.tsv",sep="\t")

###############
# B7) PHYLOGENETIC SIGNAL
###############

##
## B7.1) Diversity
##

#Load host phylogeny
host_tree <- read.tree("Data/host_phylogeny.tre")

#Load diversity data
meta_dR.raw <- readRDS("Results/RDS/meta_dR.all.RData")
meta_dRE.raw <- readRDS("Results/RDS/meta_dRE.all.RData")
meta_dRER.raw <- readRDS("Results/RDS/meta_dRER.all.RData")
results <- data.frame(dR=meta_dR.raw$TE,dRE=meta_dRE.raw$TE,dRER=meta_dRER.raw$TE)
rownames(results) <- meta_dR.raw$studlab

#Create phylo4d object
tree4d <-phylo4d(host_tree, tip.data = results)

#Run phylogenetic signal analysis
phyloSignal(tree4d, methods = "Cmean", reps = 9999, W = NULL)

##
## B7.2) Composition
##

#Load host phylogeny
host_tree <- read.tree("Data/host_phylogeny.tre")

#Load compositional data
results <- read.table("Results/summary_betadis.tsv")

#Create phylo4d object
tree4d <-phylo4d(host_tree, tip.data = results)

#Run phylogenetic signal analysis
phyloSignal(tree4d, methods = "Cmean", reps = 9999, W = NULL)

###############
# B8) COMPOSITIONAL DIFFERENCES (ANALYSIS OF VARIANCE)
###############

# Load required data
metadata.filtered <- read.table("Data/metadata.filtered.csv", sep =";",row.names=1)
code.list <- as.character(read.table("Data/sp_code.txt")[,1])
capwild.tree <- read.tree("Data/genustree.tre")

##
## B8.1) Permutational analysis of variance based on richness (dR)
##

permanovaR_results <- c()
permutestR_results <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  ##Filter the metadata
  metadata.filtered.subset <- metadata.filtered[which(metadata.filtered[,1] %in% colnames(final.table)),]
  #Check whether both lists are identical
  identical(sort(colnames(final.table)),sort(as.character(metadata.filtered.subset[,"Sample"])))
  pairdis.R <- pair_dis(final.table,qvalue=0,hierarchy=metadata.filtered.subset[,c("Sample","Origin")])
  u0n <- pairdis.R$L1_UqN
  u0n.dist <- as.dist(u0n)
  ps.disper.u0n.origin <- betadisper(u0n.dist, metadata.filtered.subset$Origin)
  permutestR <- permutest(ps.disper.u0n.origin, pairwise = TRUE)
  permutestR_results <- append(permutestR_results,permutestR[1])
  permanovaR <- adonis(u0n.dist ~ Origin, data = metadata.filtered.subset, permutations = 999)
  print(code)
  permanovaR_results <- append(permanovaR_results,permanovaR[1])
}
names(permanovaR_results) <- code.list
names(permutestR_results) <- code.list
saveRDS(permanovaR_results, "Results/RDS/permanova_dR.RData")
saveRDS(permutestR_results, "Results/RDS/permutest_dR.RData")

##
## B8.2) Permutational analysis of variance based on richness+eveness+regularity (dRER)
##

permanovaRER_results <- c()
permutestRER_results <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  tree_filtered <- match_data(final.table,capwild.tree,output="tree")
  ##Filter the metadata
  metadata.filtered.subset <- metadata.filtered[which(metadata.filtered[,1] %in% colnames(final.table)),]
  #Check whether both lists are identical
  identical(sort(colnames(final.table)),sort(as.character(metadata.filtered.subset[,"Sample"])))

  pairdis.RER <- pair_dis(final.table,qvalue=1,hierarchy=metadata.filtered.subset[,c("Sample","Origin")], tree=tree_filtered)
  u1n <- pairdis.RER$L1_UqN
  u1n.dist <- as.dist(u1n)
  ps.disper.u1n.origin <- betadisper(u1n.dist, metadata.filtered.subset$Origin)
  permutestRER <- permutest(ps.disper.u1n.origin, pairwise = TRUE)
  permutestRER_results <- append(permutestRER_results,permutestRER[1])
  permanovaRER <- adonis(u1n.dist ~ Origin, data =metadata.filtered.subset, permutations = 999)
  permanovaRER_results <- append(permanovaRER_results,permanovaRER[1])
}
names(permanovaRER_results) <- code.list
names(permutestRER_results) <- code.list
saveRDS(permanovaRER_results,paste("Results/RDS/permanova_dRER.RData",sep=""))
saveRDS(permutestRER_results,paste("Results/RDS/permutest_dRER.RData",sep=""))

###############
# B9) NULL MODELS (python)
###############

countdata <- read.csv("genustable.csv",row.names=1)
countdata <- tss(countdata)
metadata <- read.csv("metadata.csv")

#Get lis of studies
study_list <- as.character(unique(metadata$Study))

#Subset per species
for (study in study_list){

#Aggregate wild
wild_samples <- as.character(metadata[metadata$Study == study & metadata$Sample.type == "Wild","Individual"])
wild <- rowMeans(countdata[,wild_samples])

#Aggregate captivity
captivity_samples <- as.character(metadata[metadata$Study == study & metadata$Sample.type == "Captivity","Individual"])
captivity <- rowMeans(countdata[,captivity_samples])

study_table <- data.frame(wild,captivity)
study_table <- round(study_table * 10000)
write.csv(study_table,paste0("data/nullmodels",study,".csv"))
}

#In python

conda activate qdiv_env

import qdiv
from qdiv import stats
from qdiv import null
import pandas as pd
import numpy as np
import matplotlib

study='VAHI' #iterate across studies
tab='data/nullmodels'+study+'.csv'

#load files
obj = qdiv.files.load(path="data/nullmodels", tab=tab, fasta='None', tree='genustree.tre', sep=',', addTaxonPrefix=False, orderSeqs=True)
# Raup-Crick dR
dR=null.rcq(obj, randomization='frequency', iterations=999, divType='Jaccard', q=0)
# Raup-Crick dRE
dRE=null.rcq(obj, randomization='frequency', iterations=999, divType='Jaccard', q=1)
# Raup-Crick dRER
dRER=null.rcq(obj, randomization='frequency', iterations=999, divType='phyl', q=1)

print(study)
print(dR)
print(dRE)
print(dRER)

###############
# B10) APROXIMATE BAYESIAN COMPUTATION
###############

# Specify the mean vector and variances for the multinormal distribution

SVD <- function(mu,sigma) {
  n <- dim(sigma)[1]
  s <- svd((sigma + t(sigma))/2) # Guarantees symmetry
  s$d <- abs(zapsmall(s$d))
  m <- sum(s$d > 0)

  # Generate normals in the lower dimension
  n.sample <- 100   # Number of realizations
  x <- matrix(rnorm(m*n.sample),nrow=m)

  # Embed in the higher dimension and apply square root
  # of sigma obtained from its SVD
  x <- rbind(x, matrix(0, nrow=n-m, ncol=n.sample))
  y <- s$u %*% diag(sqrt(s$d)) %*% x + mu
  t <- t(y)
  colnames(t) <- c("W","S","C")
  return(t)
}


#Scenario A
mu <- c(55,35,10)
sigma <- matrix(c(40,-20,-20,-20,40,-20,-20,-20,40),3,3)
A_scenario <- SVD(mu,sigma)

#Scenario B
mu <- c(10,35,55)
sigma <- matrix(c(40,-20,-20,-20,40,-20,-20,-20,40),3,3)
B_scenario <- SVD(mu,sigma)

#Scenario C
mu <- c(10,80,10)
sigma <- matrix(c(40,-20,-20,-20,40,-20,-20,-20,40),3,3)
C_scenario <- SVD(mu,sigma)

#Scenario D
mu <- c(30,40,30)
sigma <- matrix(c(40,-20,-20,-20,40,-20,-20,-20,40),3,3)
D_scenario <- SVD(mu,sigma)

#Scenario E
mu <- c(45,10,45)
sigma <- matrix(c(40,-20,-20,-20,40,-20,-20,-20,40),3,3)
E_scenario <- SVD(mu,sigma)

#Concatenate all scenarios
values <- rbind(A_scenario,B_scenario,C_scenario,D_scenario,E_scenario)
scenarios <- c(rep("A",100),rep("B",100),rep("C",100),rep("D",100),rep("E",100))
origin <- read.csv("origin.csv",row.names=1)

# Plot density plot
##########
# Plot densities of the three components for a given scenario (example C)
all <- as.data.frame(cbind(c(C_scenario[,1],C_scenario[,2],C_scenario[,3]),c(rep("W",100),rep("S",100),rep("C",100))))
colnames(all) <- c("values","group")
all[,1] <- as.numeric(as.character(all[,1]))
all[,2] <- as.factor(all[,2])
ggplot(all, aes(values, fill = group)) + geom_density(alpha = 0.9)

# Plot simulated and observed data
##########

allvalues <- rbind(A_scenario,B_scenario,C_scenario,D_scenario,E_scenario,round(origin,0))
allgroups <- c(rep("A",100),rep("B",100),rep("C",100),rep("D",100),rep("E",100),rownames(origin))
allvalues$G <- allgroups
allvalues$Size <- c(rep(1,100),rep(1,100),rep(1,100),rep(1,100),rep(1,100),rep(2,25))

col = c("#e5e5e5","#5d4d72","#72845c","#417491","#c5c5c5","#a45d6d","#747474","#a47273","#97cd6b","#4f4f4f","#333333","#847f7c","#a18381","#fec689","#a94422","#b39186","#949599","#f9c160","#7358a5","#b69e2e","#f26549","#f69779","#c48b56","#cdc774","#c6c176","#927452","#938f60","#c5db6e","#6169a7","#71c4ee")
pca_res <- prcomp(allvalues[,c(1:3)], scale. = TRUE)
pca_res2 <- data.frame(pca_res$x,Group=as.character(allvalues$G),Size=as.numeric(allvalues$Size))

pdf("pca.pdf",width=8,height=5)
ggplot(pca_res2,aes(x=PC1,y=PC2,col=Group,size=Size))+
   geom_point(alpha=0.8)+
   scale_color_manual(values = col)+
   theme_classic()
dev.off()

# Perform ABC
##########

sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")

modtable <- c()
for (sp in sp_sorted){
modsel <- postpr(origin[sp,], scenarios, values, tol=.8, method="mnlogistic")
modtable <-rbind(modtable,modsel$pred)
}
rownames(modtable) <- sp_sorted
round(modtable,3)
modtable

# Plot heatmap of posterior probabilities
coul <- colorRampPalette(brewer.pal(8, "YlOrBr"))(25)
pdf("heatmap.pdf",width=8,height=5)
heatmap(round(modtable,3),Colv=NA, col = coul)
dev.off()

#
res.gfit.bott=gfit(target=origin["VAHI",], sumstat=values[scenarios=="C",], statistic=mean, nb.replicate=100)
plot(res.gfit.bott, main="Histogram under H0")

gfitpca(target=origin["VAHI",], sumstat=values, index=scenarios, cprob=.05)

###############
# B11) DISTRIBUTION OF THE ORIGIN OF THE DETECTED GENERA
###############

# Load required data
metadata.filtered <- read.table("Data/metadata.filtered.csv", sep =";",row.names=1)
code.list <- as.character(read.table("Data/sp_code.txt")[,1])

shared.taxa <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  ##Filter the metadata
  metadata.filtered.subset <- metadata.filtered[which(metadata.filtered[,1] %in% colnames(final.table)),c("Sample","Origin")]
  filtered.pecent <- tss(final.table)*100
  filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
  table.W <- tibble::rownames_to_column(filtered.pecent.t, "Sample")
  #Wild taxa
  wild.capt <- merge(table.W,metadata.filtered.subset,by="Sample")
  wild.capt <- wild.capt[,-c(1)]
  wild <- wild.capt[wild.capt$Origin == "Wild",]
  wild <- wild[,-ncol(wild)]
  wild <- wild[,colSums(wild) > 0]
  wild.taxa <- colnames(wild)
  #Captive taxa
  captive <- wild.capt[wild.capt$Origin == "Captivity",]
  captive <- captive[,-ncol(captive)]
  captive <- captive[,colSums(captive) > 0]
  captive.taxa <- colnames(captive)
  #Distribution
  both.taxa <- intersect(wild.taxa,captive.taxa)
  wildonly.taxa <- wild.taxa[!wild.taxa %in% captive.taxa]
  captiveonly.taxa <- captive.taxa[! captive.taxa %in% wild.taxa]
  alltaxa <- c(both.taxa,wildonly.taxa,captiveonly.taxa)
  #Vector
  shared <- cbind(length(wildonly.taxa),length(both.taxa),length(captiveonly.taxa),length(wildonly.taxa)/length(alltaxa)*100,length(both.taxa)/length(alltaxa)*100,length(captiveonly.taxa)/length(alltaxa)*100)
  shared.taxa <- rbind(shared.taxa,shared)
}
colnames(shared.taxa) <- c("wild_only", "both", "captive_only","wild_only_per", "both_per", "captive_only_per")
rownames(shared.taxa) <- code.list
write.table(shared.taxa,"Results/shared_taxa.tsv")

#Mean stats
colMeans(shared.taxa)

###############
# B12) PAIRWISE DISSIMILARITY CORRELATION BETWEEN WILD AND CAPTIVE
###############

##
## B10.1) dR
##

# Load required data
count.table.all.g <- read.table("Tables/count_Genus_all.tsv")
metadata.filtered <- read.table("Data/metadata.filtered.csv", sep =";",row.names=1)
host_tree <- read.tree("Data/host_phylogeny.tre")
sp_code <- read.table("Data/sp_code.txt")

#Pairwise dissimilarity matrix of wild individuals
count.table.wild <- count.table.all.g[,colnames(count.table.all.g) %in% metadata.filtered[metadata.filtered$Origin == "Wild","Sample"]]
pair_dis_dR_wild <- pair_dis(count.table.wild,qvalue=0,hierarchy=metadata.filtered[metadata.filtered$Origin == "Wild",c("Sample","Species")], tree = tree_filtered))
saveRDS(pair_dis_dR_wild,"Results/RDS/pairdis_dR.wild.RData")

#Pairwise dissimilarity matrix of captive individuals
count.table.captive <- count.table.all.g[,colnames(count.table.all.g) %in% metadata.filtered[metadata.filtered$Origin == "Captivity","Sample"]]
pair_dis_dR_captive <- pair_dis(count.table.captive,qvalue=0,hierarchy=metadata.filtered[metadata.filtered$Origin == "Captivity",c("Sample","Species")], tree = tree_filtered))
saveRDS(pair_dis_dR_captive,"Results/RDS/pairdis_dR.captive.RData")

#Load dissimilarity files
pair_dis_dR_wild <- readRDS("Results/RDS/pairdis_dR.wild.RData")
pair_dis_dR_captive <- readRDS("Results/RDS/pairdis_dR.captive.RData")

#Prepare dissimilarity matrices
pair_dis_dR_wild_L2_UqN <- pair_dis_dR_wild$L2_UqN
pair_dis_dR_wild_L2_UqN <- pair_dis_dR_wild_L2_UqN[!is.na(pair_dis_dR_wild_L2_UqN)]
pair_dis_dR_captive_L2_UqN <- pair_dis_dR_captive$L2_UqN
pair_dis_dR_captive_L2_UqN <- pair_dis_dR_captive_L2_UqN[!is.na(pair_dis_dR_captive_L2_UqN)]

#Obtain TMRCA and sort data
myr_split <- cophenetic.phylo(host_tree)/2
myr_split <- myr_split[colnames(pair_dis_dR_wild$L2_UqN),rownames(pair_dis_dR_wild$L2_UqN)]
myr_split_vector <- myr_split[lower.tri(myr_split, diag = FALSE)]

#Correlation table
cortable <- as.data.frame(cbind(pair_dis_dR_wild_L2_UqN,pair_dis_dR_captive_L2_UqN,myr_split_vector))
colnames(cortable) <- c("Y","X","D")

#Correlation test
cor.test(cortable$Y,cortable$X)

#Correlation plot
corplot <- ggplot(cortable, aes(x=X, y=Y, color=D)) +
  geom_point() +
  scale_color_gradient(low="#009bb4", high="#efcf00") +
  geom_smooth(method=lm , color="#cc2b7c", fill="#cc2b7c", se=TRUE) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey'))
ggsave("Results/Plots/dissimilarity_correlation_dR.pdf",corplot)

##
## B12.2) dRE
##

# Load required data
count.table.all.g <- read.table("Tables/count_Genus_all.tsv")
metadata.filtered <- read.table("Data/metadata.filtered.csv", sep =";",row.names=1)
host_tree <- read.tree("Data/host_phylogeny.tre")
sp_code <- read.table("Data/sp_code.txt")

#Pairwise dissimilarity matrix of wild individuals
count.table.wild <- count.table.all.g[,colnames(count.table.all.g) %in% metadata.filtered[metadata.filtered$Origin == "Wild","Sample"]]
pair_dis_dRE_wild <- pair_dis(count.table.wild,qvalue=1,hierarchy=metadata.filtered[metadata.filtered$Origin == "Wild",c("Sample","Species")])
saveRDS(pair_dis_dRE_wild,"Results/RDS/pairdis_dRE.wild.RData")

#Pairwise dissimilarity matrix of captive individuals
count.table.captive <- count.table.all.g[,colnames(count.table.all.g) %in% metadata.filtered[metadata.filtered$Origin == "Captivity","Sample"]]
pair_dis_dRE_captive <- pair_dis(count.table.captive,qvalue=1,hierarchy=metadata.filtered[metadata.filtered$Origin == "Captivity",c("Sample","Species")])
saveRDS(pair_dis_dRE_captive,"Results/RDS/pairdis_dRE.captive.RData")

#Load dissimilarity files
pair_dis_dRE_wild <- readRDS("Results/RDS/pairdis_dRE.wild.RData")
pair_dis_dRE_captive <- readRDS("Results/RDS/pairdis_dRE.captive.RData")

#Prepare dissimilarity matrices
pair_dis_dRE_wild_L2_UqN <- pair_dis_dRE_wild$L2_UqN
pair_dis_dRE_wild_L2_UqN <- pair_dis_dRE_wild_L2_UqN[!is.na(pair_dis_dRE_wild_L2_UqN)]
pair_dis_dRE_captive_L2_UqN <- pair_dis_dRE_captive$L2_UqN
pair_dis_dRE_captive_L2_UqN <- pair_dis_dRE_captive_L2_UqN[!is.na(pair_dis_dRE_captive_L2_UqN)]

#Obtain TMRCA and sort data
myr_split <- cophenetic.phylo(host_tree)/2
myr_split <- myr_split[colnames(pair_dis_dRE_wild$L2_UqN),rownames(pair_dis_dRE_wild$L2_UqN)]
myr_split_vector <- myr_split[lower.tri(myr_split, diag = FALSE)]

#Correlation table
cortable <- as.data.frame(cbind(pair_dis_dRE_wild_L2_UqN,pair_dis_dRE_captive_L2_UqN,myr_split_vector))
colnames(cortable) <- c("Y","X","D")

#Correlation test
cor.test(cortable$Y,cortable$X)

#Correlation plot
corplot <- ggplot(cortable, aes(x=X, y=Y, color=D)) +
  geom_point() +
  scale_color_gradient(low="#009bb4", high="#efcf00") +
  geom_smooth(method=lm , color="#cc2b7c", fill="#cc2b7c", se=TRUE) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey'))
ggsave("Results/Plots/dissimilarity_correlation_dRE.pdf",corplot)

##
## B12.2) dRER
##

# Load required data
count.table.all.g <- read.table("Tables/count_Genus_all.tsv")
metadata.filtered <- read.table("Data/metadata.filtered.csv", sep =";",row.names=1)
host_tree <- read.tree("Data/host_phylogeny.tre")
sp_code <- read.table("Data/sp_code.txt")
capwild.tree <- read.tree("Data/genustree.tre")

#Pairwise dissimilarity matrix of wild individuals
count.table.wild <- count.table.all.g[,colnames(count.table.all.g) %in% metadata.filtered[metadata.filtered$Origin == "Wild","Sample"]]
tree_filtered <- match_data(count.table.wild,capwild.tree,output="tree")
pair_dis_dRER_wild <- pair_dis(count.table.wild,qvalue=1,hierarchy=metadata.filtered[metadata.filtered$Origin == "Wild",c("Sample","Species")],tree=tree_filtered)
saveRDS(pair_dis_dRER_wild,"Results/RDS/pairdis_dRER.wild.RData")

#Pairwise dissimilarity matrix of captive individuals
count.table.captive <- count.table.all.g[,colnames(count.table.all.g) %in% metadata.filtered[metadata.filtered$Origin == "Captivity","Sample"]]
tree_filtered <- match_data(count.table.captive,capwild.tree,output="tree")
pair_dis_dRER_captive <- pair_dis(count.table.captive,qvalue=1,hierarchy=metadata.filtered[metadata.filtered$Origin == "Captivity",c("Sample","Species")],tree=tree_filtered)
saveRDS(pair_dis_dRER_captive,"Results/RDS/pairdis_dRER.captive.RData")

#Load dissimilarity files
pair_dis_dRER_wild <- readRDS("Results/RDS/pairdis_dRER.wild.RData")
pair_dis_dRER_captive <- readRDS("Results/RDS/pairdis_dRER.captive.RData")

#Prepare dissimilarity matrices
pair_dis_dRER_wild_L2_UqN <- pair_dis_dRER_wild$L2_UqN
pair_dis_dRER_wild_L2_UqN <- pair_dis_dRER_wild_L2_UqN[!is.na(pair_dis_dRER_wild_L2_UqN)]
pair_dis_dRER_captive_L2_UqN <- pair_dis_dRER_captive$L2_UqN
pair_dis_dRER_captive_L2_UqN <- pair_dis_dRER_captive_L2_UqN[!is.na(pair_dis_dRER_captive_L2_UqN)]

#Obtain TMRCA and sort data
myr_split <- cophenetic.phylo(host_tree)/2
myr_split <- myr_split[colnames(pair_dis_dRER_wild$L2_UqN),rownames(pair_dis_dRER_wild$L2_UqN)]
myr_split_vector <- myr_split[lower.tri(myr_split, diag = FALSE)]

#Correlation table
cortable <- as.data.frame(cbind(pair_dis_dRER_wild_L2_UqN,pair_dis_dRER_captive_L2_UqN,myr_split_vector))
colnames(cortable) <- c("Y","X","D")

#Correlation test
cor.test(cortable$Y,cortable$X)

#Correlation plot
corplot <- ggplot(cortable, aes(x=X, y=Y, color=D)) +
  geom_point() +
  scale_color_gradient(low="#009bb4", high="#efcf00") +
  geom_smooth(method=lm , color="#cc2b7c", fill="#cc2b7c", se=TRUE) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey'))
ggsave("Results/Plots/dissimilarity_correlation_dRER.pdf",corplot)


###############
# B13) NEAREST DATASET (only dR)
###############

# Load required data
count.table.all.g <- read.table("Tables/count_Genus_all.tsv")
metadata.filtered <- read.table("Data/metadata.filtered.csv", sep =";",row.names=1)

#Add group variable to the metadata matrix by combining species and origin
metadata.filtered$Group <- paste(metadata.filtered$Species,metadata.filtered$Origin,sep="_")

#Compute pairwise dissimilarities
pair_dis_dR.groups <- pair_dis(count.table.all.g,qvalue=0,hierarchy=metadata.filtered[,c("Sample","Group")])
saveRDS(pair_dis_dR.groups,"Results/RDS/pairdis_dR.groups.RData")
pairdis_dR.groups <- readRDS("Results/RDS/pairdis_dR.groups.RData")

#Obtain distance matrix
distmatrix <- pairdis_dR.groups$L2_UqN

#Calculate position of the captive or wild counterpart
vector <- c()
for (g in as.character(unique(metadata.filtered$Group))){
sp <- gsub("_.*","",g)
origin <- gsub(".*_","",g)

if (origin == "Wild") {
ref=paste(sp,"Captivity",sep="_")
}

if (origin == "Captivity") {
ref=paste(sp,"Wild",sep="_")
}

distmatrix_sub <- c(distmatrix[,g],distmatrix[g,])
distmatrix_sub <- distmatrix_sub[!is.na(distmatrix_sub)]
distmatrix_sub <- sort(distmatrix_sub)
position <- which(names(distmatrix_sub) == ref)
vector <- c(vector,position)
}
names(vector) <- as.character(unique(metadata.filtered$Group))

#Percentage of groups in which the wild or captive counterpart is the most similar group
length(vector[vector == 1]) / length(vector) * 100

#Percentage of groups in which the wild or captive counterpart is not the most similar group
length(vector[vector > 1]) / length(vector) * 100

###############
# B14) HIERARCHICAL CLUSTERING AND TOPOLOGICAL DIFFERENCES
###############

##
## B14.1) dR
##

#Load dissimilarity files
pair_dis_dR_wild <- readRDS("Results/RDS/pairdis_dR.wild.RData")
pair_dis_dR_captive <- readRDS("Results/RDS/pairdis_dR.captive.RData")

#Load host phylogeny
host_tree <- read.tree("Data/host_phylogeny.tre")

#Hierarchical clustering
hclust_wild <- as.phylo(hclust(as.dist(pair_dis_dR_wild$L2_UqN), method="average"))
hclust_captive <- as.phylo(hclust(as.dist(pair_dis_dR_captive$L2_UqN),method="average"))

#Calculate tree dissimilarities
JaccardRobinsonFoulds(hclust_wild, hclust_captive)
JaccardRobinsonFoulds(host_tree, hclust_wild)
JaccardRobinsonFoulds(host_tree, hclust_captive)

# Phylosymbiosis effect size comparison
study_list <- host_tree$tip.label

valuetable <- c()
for (i in c(1:100)){
selected <- study_list[sample(c(1:25),10)]
host_tree_selected <- keep.tip(host_tree,selected)
hclust_wild_selected <- keep.tip(hclust_wild,selected)
hclust_captive_selected <- keep.tip(hclust_captive,selected)
row <- c(WC=JaccardRobinsonFoulds(hclust_wild_selected, hclust_captive_selected),HW=JaccardRobinsonFoulds(host_tree_selected, hclust_wild_selected),HC=JaccardRobinsonFoulds(host_tree_selected, hclust_captive_selected))
valuetable <- rbind(valuetable,row)
}

boxplot(valuetable)
t.test(valuetable[,2],valuetable[,3])

#Plot tanglegram
pdf("Results/Plots/clustering_tanglegram_dR.pdf",width=8,height=6)
tanglegram(untangle_labels(hclust_wild, hclust_captive,method="random"), highlight_distinct_edges = FALSE,highlight_branches_lwd = FALSE)
dev.off()

##
## B14.2) dRE
##

# Load dissimilarity files
pair_dis_dRE_wild <- readRDS("Results/RDS/pairdis_dRE.wild.RData")
pair_dis_dRE_captive <- readRDS("Results/RDS/pairdis_dRE.captive.RData")

# Load host phylogeny
host_tree <- read.tree("Data/host_phylogeny.tre")

# Hierarchical clustering
hclust_wild <- as.phylo(hclust(as.dist(pair_dis_dRE_wild$L2_UqN), method="average"))
hclust_captive <- as.phylo(hclust(as.dist(pair_dis_dRE_captive$L2_UqN),method="average"))

# Calculate tree dissimilarities
JaccardRobinsonFoulds(hclust_wild, hclust_captive)
JaccardRobinsonFoulds(host_tree, hclust_wild)
JaccardRobinsonFoulds(host_tree, hclust_captive)

# Phylosymbiosis effect size comparison
study_list <- host_tree$tip.label

valuetable <- c()
for (i in c(1:100)){
selected <- study_list[sample(c(1:25),10)]
host_tree_selected <- keep.tip(host_tree,selected)
hclust_wild_selected <- keep.tip(hclust_wild,selected)
hclust_captive_selected <- keep.tip(hclust_captive,selected)
row <- c(WC=JaccardRobinsonFoulds(hclust_wild_selected, hclust_captive_selected),HW=JaccardRobinsonFoulds(host_tree_selected, hclust_wild_selected),HC=JaccardRobinsonFoulds(host_tree_selected, hclust_captive_selected))
valuetable <- rbind(valuetable,row)
}

boxplot(valuetable)
t.test(valuetable[,2],valuetable[,3])

#Plot tanglegram
pdf("Results/Plots/clustering_tanglegram_dRE.pdf",width=8,height=6)
tanglegram(untangle_labels(hclust_wild, hclust_captive,method="random"), highlight_distinct_edges = FALSE,highlight_branches_lwd = FALSE)
dev.off()

##
## B14.3) dRER
##

#Load dissimilarity files
pair_dis_dRER_wild <- readRDS("Results/RDS/pairdis_dRER.wild.RData")
pair_dis_dRER_captive <- readRDS("Results/RDS/pairdis_dRER.captive.RData")

#Load host phylogeny
host_tree <- read.tree("Data/host_phylogeny.tre")

# Hierarchical clustering
hclust_wild <- as.phylo(hclust(as.dist(pair_dis_dRER_wild$L2_UqN), method="average"))
hclust_captive <- as.phylo(hclust(as.dist(pair_dis_dRER_captive$L2_UqN),method="average"))

# Calculate tree dissimilarities
JaccardRobinsonFoulds(hclust_wild, hclust_captive)
JaccardRobinsonFoulds(host_tree, hclust_wild)
JaccardRobinsonFoulds(host_tree, hclust_captive)

# Phylosymbiosis effect size comparison
study_list <- host_tree$tip.label

valuetable <- c()
for (i in c(1:100)){
selected <- study_list[sample(c(1:25),10)]
host_tree_selected <- keep.tip(host_tree,selected)
hclust_wild_selected <- keep.tip(hclust_wild,selected)
hclust_captive_selected <- keep.tip(hclust_captive,selected)
row <- c(WC=JaccardRobinsonFoulds(hclust_wild_selected, hclust_captive_selected),HW=JaccardRobinsonFoulds(host_tree_selected, hclust_wild_selected),HC=JaccardRobinsonFoulds(host_tree_selected, hclust_captive_selected))
valuetable <- rbind(valuetable,row)
}

boxplot(valuetable)
t.test(valuetable[,2],valuetable[,3])

#Plot tanglegram
pdf("Results/Plots/clustering_tanglegram_dRE.pdf",width=8,height=6)
tanglegram(untangle_labels(hclust_wild, hclust_captive,method="random"), highlight_distinct_edges = FALSE,highlight_branches_lwd = FALSE)
dev.off()

###############
# B15) TAXONOMIC ENRICHMENT ANALYSIS
###############

countdata <- read.csv("genustable.csv",row.names=1)
countdata <- tss(countdata)
metadata <- read.csv("metadata.csv")
colourcodes <- read.table("sp_code.txt", header=TRUE)

#Modify taxonomy names to input to metamicrobiomeR
rownames(countdata) <- gsub("^","k__",rownames(countdata))

#Merge count data and metadata
alldata <- merge(metadata[,c("Sample","Species","Sample.type","Origin")],t(countdata),by.x=1,by.y=0)

colnames(alldata) <- gsub(" ","-",colnames(alldata))
alldata$Origin <- as.factor(alldata$Origin)

#Run enrichment analysis in each study and add to column
studylist <- unique(metadata$Species)

taxacompare_all <- c()
for (Study in studylist){
  studydata <- alldata[alldata$Species == Study,]
  taxacompare <- taxa.compare(taxtab=studydata, propmed.rel="gamlss", comvar="Origin", adjustvar="Origin", longitudinal="no")
  taxacompare$study <- Study
  taxacompare_all <- rbind.fill(taxacompare_all,taxacompare)
}

#Run differetial abundance meta-analysis
metataxa <- meta.taxa(taxcomdat = taxacompare_all, summary.measure = "RR", pool.var = "id", studylab = "study", backtransform = FALSE, percent.meta = 0.33, p.adjust.method = "fdr")

#Color test
oribacterium <- countdata["k__Roseburia",]
oribacterium2 <- merge(t(oribacterium),metadata[,c("Sample","Species","Sample.type","Origin")],by.x=0,by.y=1)
aggregate(oribacterium2[,2],by=list(oribacterium2[,5]),FUN=mean)

# Histogram
#####
taxonsubset <- metataxa$random$OriginWild$id
studylist <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")

allestimates <- c()
for(taxon in taxonsubset){
taxoncompare <- taxacompare_all[taxacompare_all$id == taxon,]
  for(study in studylist){
    estimate <- taxoncompare[taxoncompare$study == study,"Estimate.OriginWild"]
    if(length(estimate)==0){estimate=0}
    row <- c(taxon,study,estimate)
    allestimates <- rbind(allestimates,row)
  }
}

allestimates <- as.data.frame(allestimates)
colnames(allestimates) <- c("Taxon","Study","Estimate")

#Order factors
taxonorder <- gsub("k__","",taxonsubset)

allestimates$Taxon <- gsub("k__","",allestimates$Taxon)
allestimates$Taxon <- factor(allestimates$Taxon, levels = rev(taxonorder))
allestimates$Study <- factor(allestimates$Study, levels = c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT"))
allestimates$Estimate <- as.numeric(allestimates$Estimate)

#Plot
pdf("heatmap.pdf",width=8,height=7)
my_colours <- colorRampPalette(c("#854292","#ffffff","#cc2b7c"))(11)
ggplot(allestimates, aes(Study, Taxon, fill= Estimate)) +
  geom_tile() +
  scale_fill_gradientn(colours = my_colours, limits=c(-6, 6)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

# Forest plot
#####

metaestimates <- metataxa$random$OriginWild
#Color code
metaestimates$enrichment <- metaestimates$estimate
metaestimates$enrichment[metaestimates$enrichment < 0] <- "Captivity"
metaestimates$enrichment[metaestimates$enrichment != "Captivity"] <- "Wild"
metaestimates$enrichment <- as.factor(metaestimates$enrichment)
#Alpha code
metaestimates$significance <- metaestimates$p.adjust
metaestimates$significance[metaestimates$significance < 0.05] <- 2
metaestimates$significance[metaestimates$significance < 0.1] <- 1
metaestimates$significance[metaestimates$significance < 1] <- 0
metaestimates$significance[metaestimates$significance == 0] <- 0.5
metaestimates$significance[metaestimates$significance == 1] <- 0.7
metaestimates$significance[metaestimates$significance == 2] <- 1
#Other adjustments
metaestimates$id <- gsub("k__","",metaestimates$id)
#Set order
rownames(metaestimates) <- metaestimates$id
metaestimates <- metaestimates[rev(taxonorder),]
metaestimates$id <- factor(metaestimates$id, levels = rev(taxonorder))

pdf("forest.pdf",width=4,height=7)
ggplot(metaestimates) +
  geom_errorbar(data=metaestimates, mapping=aes(y=id, xmin=ll, xmax=ul, col=enrichment, alpha=significance), width=0.2, size=1) +
  scale_alpha_identity() +
  scale_color_manual("enrichment", values=c("#854292", "#cc2b7c")) +
  geom_point(data=metaestimates, mapping=aes(x=estimate, y=id, fill=enrichment), size=2, shape=21) +
  scale_fill_manual("enrichment", values=c("#854292", "#cc2b7c")) +
  theme_minimal()
dev.off()

# Relative abundances
#####

taxonsubset <- metataxa$random$OriginWild$id
studylist <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")

allabundances <- c()
for(taxon in taxonsubset){
  for(study in studylist){
    value <- sum(alldata[alldata$Species == study,taxon])
    row <- c(taxon,study,value)
    allabundances <- rbind(allabundances,row)
  }
}

allabundances <- as.data.frame(allabundances)
colnames(allabundances) <- c("Taxon","Study","Abundance")
allabundances$Taxon <- gsub("k__","",allabundances$Taxon)
allabundances$Abundance <- as.numeric(allabundances$Abundance)
allabundances$Abundance <- allabundances$Abundance / 25
allabundances$Taxon <- factor(allabundances$Taxon, levels = rev(taxonorder))

pdf("barplot.pdf",width=7,height=7)
ggplot(allabundances, aes(fill=Study, y=Abundance, x=Taxon)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual("Study", values=c(colourcodes[,3]))+
    coord_flip() +
    theme_minimal()
dev.off()

#Total values
aggregate(allabundances$Abundance,by=list(allabundances$Taxon),FUN=sum)
aggregate(allabundances$Abundance,by=list(allabundances$Taxon),function(x) sum(x != 0))

###############
# B16) NMDS plot
###############

# Load required data
metadata.filtered <- read.table("Data/metadata.filtered.csv", sep =";",row.names=1)

# B16.1) NMDS plot dR

#Load pairwise dissimilarity file
pairdis_dR.all <- readRDS("Results/RDS/pairdis_dR.all.RData")

#Extract dissimilarity matrix
dismatrix <- pairdis_dR.all$L1_CqN

#Generate NMDS object
values.NMDS <- metaMDS(as.dist(dismatrix), k = 2, trymax = 200)
NMDS=data.frame(x=values.NMDS$point[,1],y=values.NMDS$point[,2],Species=as.factor(metadata.filtered$Species),Origin=as.factor(metadata.filtered$Origin),Group=as.factor(paste(metadata.filtered$Species,metadata.filtered$Origin,sep="_")))
NMDS$Sample <- rownames(NMDS)

#Add colour
sp_code <- read.table("Data/sp_code.txt")
colnames(sp_code) <- c("Species","Description","Color")
NMDS=merge(NMDS,sp_code[,c(1,3)],by="Species")

#Find group centroids
NMDS.centroids=aggregate(NMDS[,c("x","y")],by=list(NMDS$Group),FUN=mean)
colnames(NMDS.centroids) <- c("Group","x_cen","y_cen")
NMDS=merge(NMDS,NMDS.centroids,by="Group")

#Find species centroids
NMDS.centroids2=aggregate(NMDS[,c("x_cen","y_cen")],by=list(NMDS$Species),FUN=mean)
colnames(NMDS.centroids2) <- c("Species","x_cen2","y_cen2")
NMDS=merge(NMDS,NMDS.centroids2,by="Species")

#Plot
nmds.plot <- ggplot(NMDS, aes(x,y,colour=Species,shape=Origin)) +
  #Centroids
  geom_point(aes(x=x_cen,y=y_cen),size=3) +
  #Lines between centroids and individuals
  geom_segment(aes(x=x_cen, y=y_cen, xend=x, yend=y), alpha=0.2) +
  #Centroid names
  geom_text_repel(aes(x=x_cen,y=y_cen,label=Group),size=3,vjust=0,force=0) +
  #Lines between centroids
  geom_segment(aes(x=x_cen, y=y_cen, xend=x_cen2, yend=y_cen2), alpha=0.8) +
  scale_colour_manual(values = as.character(sp_code[unique(NMDS$Species),3])) +
  #Theme
  theme(panel.background = element_rect(fill = 'white', colour = 'grey'))

ggsave("Results/Plots/NMDS_dR.pdf",nmds.plot, width = 10, height = 6)

# B16.2) NMDS plot dRER

#Load pairwise dissimilarity file
pairdis_dR.all <- readRDS("Results/RDS/pairdis_dRER.all.RData")

#Extract dissimilarity matrix
dismatrix <- pairdis_dR.all$L1_CqN

#Generate NMDS object
values.NMDS <- metaMDS(as.dist(dismatrix), k = 2, trymax = 200)
NMDS=data.frame(x=values.NMDS$point[,1],y=values.NMDS$point[,2],Species=as.factor(metadata.filtered$Species),Origin=as.factor(metadata.filtered$Origin),Group=as.factor(paste(metadata.filtered$Species,metadata.filtered$Origin,sep="_")))
NMDS$Sample <- rownames(NMDS)

#Add colour
sp_code <- read.table("Data/sp_code.txt")
colnames(sp_code) <- c("Species","Description","Color")
NMDS=merge(NMDS,sp_code[,c(1,3)],by="Species")

#Find group centroids
NMDS.centroids=aggregate(NMDS[,c("x","y")],by=list(NMDS$Group),FUN=mean)
colnames(NMDS.centroids) <- c("Group","x_cen","y_cen")
NMDS=merge(NMDS,NMDS.centroids,by="Group")

#Find species centroids
NMDS.centroids2=aggregate(NMDS[,c("x_cen","y_cen")],by=list(NMDS$Species),FUN=mean)
colnames(NMDS.centroids2) <- c("Species","x_cen2","y_cen2")
NMDS=merge(NMDS,NMDS.centroids2,by="Species")

#Plot
nmds.plot <- ggplot(NMDS, aes(x,y,colour=Species,shape=Origin)) +
  #Centroids
  geom_point(aes(x=x_cen,y=y_cen),size=3) +
  #Lines between centroids and individuals
  geom_segment(aes(x=x_cen, y=y_cen, xend=x, yend=y), alpha=0.2) +
  #Centroid names
  geom_text_repel(aes(x=x_cen,y=y_cen,label=Group),size=3,vjust=0,force=0) +
  #Lines between centroids
  geom_segment(aes(x=x_cen, y=y_cen, xend=x_cen2, yend=y_cen2), alpha=0.8) +
  scale_colour_manual(values = as.character(sp_code[unique(NMDS$Species),3])) +
  #Theme
  theme(panel.background = element_rect(fill = 'white', colour = 'grey'))

ggsave("Results/Plots/NMDS_dRER.pdf",nmds.plot, width = 10, height = 6)
