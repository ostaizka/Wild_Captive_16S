#Load the following packages
library(reshape2)
library(gplots)
library(purrr)
library(hilldiv)
library(ape)
library("ggplot2")
library(vegan)
library("RColorBrewer")
library(tidyverse)
library(metagenomeSeq)
library(metagMisc)
library(phytools)
library(dplyr)
library("ggpubr")
library(dmetar)
library(meta)

#Directory
setwd("~/github/Wild_Captive_16S")

#### Load the metadata and create hierarchy table for diversity analysis ####
metadata <- read.table("Data/metadata.tsv", sep =";")
hierarchy_all <- tibble::rownames_to_column(metadata, "Samples")

#### Get code names ####
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
code.list <- list.dirs()

#### Get genus table ####
for (code in code.list){
  print(code)
  #Load files at genus level
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

#### Filtering the tables ####
#Load all the tables created
for (code in code.list){
  print(code)
  count.table <- read.table(paste("Tables/count_genus_",code,".tsv",sep=""))
  assign(paste("count.table.",code,sep=""),get("count.table"))
}

#Filter out samples with low sequencing depth (1000)
for (code in code.list){
  count.filtdepth <- depth_filt(get(paste("count.table.",code,sep="")),1000)
  count.filtdepth <- count.filtdepth[rowSums(count.filtdepth)>0, ]
  assign(paste("count.filtdepth.",code,sep=""),get("count.filtdepth"))
  write.table(get(paste("count.filtdepth.",code,sep="")),paste("Tables/countfiltdepth_",code,".tsv",sep=""))
}

## Remove NA taxa from Genus table
taxonlist_NA <- read.table("Data/NA_taxa.txt")
taxonlist_NA <-  as.character(taxonlist_NA[,1])

for (code in code.list){
  count.filtdepth <- read.table(paste("Tables/countfiltdepth_",code,".tsv",sep=""))
  count.filtered <- count.filtdepth[! rownames(count.filtdepth) %in% taxonlist_NA,]
  assign(paste("count.filtered.",code,sep=""),get("count.filtered"))
  write.table(get(paste("count.filtered.",code,sep="")),paste("Tables/countfiltered_",code,".tsv",sep=""))
}

####Diversity analysis based on abundace: using hilldiv####
##R
summary_R <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  divqR <- div_test(final.table,qvalue=0,hierarchy=hierarchy[,c(1,5)])
  divqR_table <- divqR$data
  divqR_wild <- divqR_table[divqR_table$Group == "Wild",]
  n_wild <- nrow(divqR_wild)
  mean_wild <- mean(divqR_wild[,c("Value")])
  sd_wild <- sd(divqR_wild[,c("Value")])
  divqR_captive <- divqR_table[divqR_table$Group == "Captivity",]
  n_captive <- nrow(divqR_captive)
  mean_captive <- mean(divqR_captive[,c("Value")])
  sd_captive <- sd(divqR_captive[,c("Value")])
  summary_R <- rbind(summary_R,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
}
colnames(summary_R) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_R) <- code.list
write.table(summary_R, "Results/summary_diversity_R.tsv")
# Meta-analysis
summary_R <- read.table("Results/summary_diversity_R.tsv")
ssummary_R <- as.data.frame(summary_R)
meta_R_ready <- tibble::rownames_to_column(summary_R,"Author")
rownames(meta_R_ready) <- meta_R_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_R_ready <- meta_R_ready[sp_sorted,]
meta_R.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                        data = meta_R_ready,
                        studlab = paste(Author),
                        comb.fixed = FALSE,
                        comb.random = TRUE,
                        method.tau = "SJ",
                        hakn = TRUE,
                        prediction = TRUE,
                        sm = "SMD")
saveRDS(meta_R.raw, "Results/RDS/metanalysis_R.RData")

pdf("Results/Plots/forest_R.pdf", width=13, height=8)
forest(meta_R.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect",
       rightlabs = c("g","95% CI","Weight"),leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),
       lab.e = "Captivity",lab.c="Wild")
dev.off()

##REH
capwild.tree <- read.tree("Data/genustree.tre")

summary_REH <- c()
for (code in code.list){
  print(code)
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  rownames(final.table) <- gsub(" ","_",rownames(final.table))
  rownames(final.table) <- gsub("[","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("]","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("(","-",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub(")","-",rownames(final.table),fixed = TRUE)
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  tree_filtered <- match_data(final.table,capwild.tree,output="tree")
  divqREH <- div_test(final.table,qvalue=1,hierarchy=hierarchy[,c(1,5)], tree=tree_filtered)
  divqREH_table <- divqREH$data
  divqREH_wild <- divqREH_table[divqREH_table$Group == "Wild",]
  n_wild <- nrow(divqREH_wild)
  mean_wild <- mean(divqREH_wild[,c("Value")])
  sd_wild <- sd(divqREH_wild[,c("Value")])
  divqREH_captive <- divqREH_table[divqREH_table$Group == "Captivity",]
  n_captive <- nrow(divqREH_captive)
  mean_captive <- mean(divqREH_captive[,c("Value")])
  sd_captive <- sd(divqREH_captive[,c("Value")])
  summary_REH <- rbind(summary_REH,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
}
colnames(summary_REH) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_REH) <- code.list
write.table(summary_REH, "Results/summary_diversity_REH.tsv")
# Meta-analysis
summary_REH <- read.table("Results/summary_diversity_REH.tsv")
summary_REH <- as.data.frame(summary_REH)
meta_REH_ready <- tibble::rownames_to_column(summary_REH,"Author")
rownames(meta_REH_ready) <- meta_REH_ready$Author
sp_sorted <- c("VAHI","APIB","RADY","CHMY","ALGI","SHCR","RHBR","PYNE","PAAN","PATR","GOGO","PEMA","PELE","TUTR","MOCH","BOGA","ELDA","CENI","EQKI","AIME","PATI","MYTR","SAHA1","SAHA2","LALT")
meta_REH_ready <- meta_REH_ready[sp_sorted,]
meta_REH.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                             data = meta_REH_ready,
                             studlab = paste(Author),
                             comb.fixed = FALSE,
                             comb.random = TRUE,
                             method.tau = "SJ",
                             hakn = TRUE,
                             prediction = TRUE,
                             sm = "SMD")
saveRDS(meta_REH.raw, "Results/RDS/metanalysis_REH.RData")

pdf("Results/Plots/forest_REH.pdf",width=13, height=8)
forest(meta_REH.raw,col.diamond = "blue",col.diamond.lines = "black",text.random = "Overall effect",
       rightlabs = c("g","95% CI","Weight"),leftlabs = c("Species", "N","Mean","SD","N","Mean","SD"),
       lab.e = "Captivity",lab.c="Wild")
dev.off()

#Subset of the data set (Primates & Cetartiodactylans)
#Primates
#R
sublist <- c("RHBR","PYNE","PAAN","PATR","GOGO")
summary_R.subset <- c()
for (code in sublist){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  divqR <- div_test(final.table,qvalue=0,hierarchy=hierarchy[,c(1,5)])
  divqR_table <- divqR$data
  divqR_wild <- divqR_table[divqR_table$Group == "Wild",]
  n_wild <- nrow(divqR_wild)
  mean_wild <- mean(divqR_wild[,c("Value")])
  sd_wild <- sd(divqR_wild[,c("Value")])
  divqR_captive <- divqR_table[divqR_table$Group == "Captivity",]
  n_captive <- nrow(divqR_captive)
  mean_captive <- mean(divqR_captive[,c("Value")])
  sd_captive <- sd(divqR_captive[,c("Value")])
  summary_R.subset <- rbind(summary_R.subset,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
}
colnames(summary_R.subset) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_R.subset) <- sublist
summary_R.subset <- as.data.frame(summary_R.subset)
meta_R_sub_ready <- tibble::rownames_to_column(summary_R.subset,"Author")
rownames(meta_R_sub_ready) <- meta_R_sub_ready$Author
sp_sorted <- c("RHBR","PYNE","PAAN","PATR","GOGO")
meta_R_sub_ready <- meta_R_sub_ready[sp_sorted,]
meta_R_sub.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                            data = meta_R_sub_ready,
                            studlab = paste(Author),
                            comb.fixed = FALSE,
                            comb.random = TRUE,
                            method.tau = "SJ",
                            hakn = TRUE,
                            prediction = TRUE,
                            sm = "SMD")

#REH
#sublist <- c("RHBR","PYNE","PAAN","PATR","GOGO")
summary_REH_sub <- c()
for (code in sublist){
  print(code)
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  rownames(final.table) <- gsub(" ","_",rownames(final.table))
  rownames(final.table) <- gsub("[","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("]","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("(","-",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub(")","-",rownames(final.table),fixed = TRUE)
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  tree_filtered <- match_data(final.table,capwild.tree,output="tree")
  divqREH <- div_test(final.table,qvalue=1,hierarchy=hierarchy[,c(1,5)], tree=tree_filtered)
  divqREH_table <- divqREH$data
  divqREH_wild <- divqREH_table[divqREH_table$Group == "Wild",]
  n_wild <- nrow(divqREH_wild)
  mean_wild <- mean(divqREH_wild[,c("Value")])
  sd_wild <- sd(divqREH_wild[,c("Value")])
  divqREH_captive <- divqREH_table[divqREH_table$Group == "Captivity",]
  n_captive <- nrow(divqREH_captive)
  mean_captive <- mean(divqREH_captive[,c("Value")])
  sd_captive <- sd(divqREH_captive[,c("Value")])
  summary_REH_sub <- rbind(summary_REH_sub,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
}
colnames(summary_REH_sub) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_REH_sub) <- sublist
summary_REH_sub <- as.data.frame(summary_REH_sub)
meta_REH_sub_ready <- tibble::rownames_to_column(summary_REH_sub,"Author")
rownames(meta_REH_sub_ready) <- meta_REH_sub_ready$Author
sp_sorted <- c("RHBR","PYNE","PAAN","PATR","GOGO")
meta_REH_sub_ready <- meta_REH_sub_ready[sp_sorted,]
meta_REH_sub_tree.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                                 data = meta_REH_sub_ready,
                                 studlab = paste(Author),
                                 comb.fixed = FALSE,
                                 comb.random = TRUE,
                                 method.tau = "SJ",
                                 hakn = TRUE,
                                 prediction = TRUE,
                                 sm = "SMD")

#Cetartiodactylans
#R
sublist <- c("TUTR","MOCH","BOGA","ELDA","CENI")
summary_R.subset <- c()
for (code in sublist){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  divR <- div_test(final.table,qvalue=0,hierarchy=hierarchy[,c(1,5)])
  divR_table <- divR$data
  divR_wild <- divR_table[divR_table$Group == "Wild",]
  n_wild <- nrow(divR_wild)
  mean_wild <- mean(divR_wild[,c("Value")])
  sd_wild <- sd(divR_wild[,c("Value")])
  divR_captive <- divR_table[divR_table$Group == "Captivity",]
  n_captive <- nrow(divR_captive)
  mean_captive <- mean(divR_captive[,c("Value")])
  sd_captive <- sd(divR_captive[,c("Value")])
  summary_R.subset <- rbind(summary_R.subset,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
}
colnames(summary_R.subset) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_R.subset) <- sublist
summary_R.subset <- as.data.frame(summary_R.subset)
meta_R_sub_ready <- tibble::rownames_to_column(summary_R.subset,"Author")
rownames(meta_R_sub_ready) <- meta_R_sub_ready$Author
sp_sorted <- c("TUTR","MOCH","BOGA","ELDA","CENI")
meta_R_sub_ready <- meta_R_sub_ready[sp_sorted,]
meta_R_sub.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                            data = meta_R_sub_ready,
                            studlab = paste(Author),
                            comb.fixed = FALSE,
                            comb.random = TRUE,
                            method.tau = "SJ",
                            hakn = TRUE,
                            prediction = TRUE,
                            sm = "SMD")

#REH
#sublist <- c("TUTR","MOCH","BOGA","ELDA","CENI")
summary_REH_sub <- c()
for (code in sublist){
  print(code)
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  rownames(final.table) <- gsub(" ","_",rownames(final.table))
  rownames(final.table) <- gsub("[","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("]","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("(","-",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub(")","-",rownames(final.table),fixed = TRUE)
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  tree_filtered <- match_data(final.table,capwild.tree,output="tree")
  divq1 <- div_test(final.table,qvalue=1,hierarchy=hierarchy[,c(1,5)], tree=tree_filtered)
  divq1_table <- divq1$data
  divq1_wild <- divq1_table[divq1_table$Group == "Wild",]
  n_wild <- nrow(divq1_wild)
  mean_wild <- mean(divq1_wild[,c("Value")])
  sd_wild <- sd(divq1_wild[,c("Value")])
  divq1_captive <- divq1_table[divq1_table$Group == "Captivity",]
  n_captive <- nrow(divq1_captive)
  mean_captive <- mean(divq1_captive[,c("Value")])
  sd_captive <- sd(divq1_captive[,c("Value")])
  summary_REH_sub <- rbind(summary_REH_sub,c(n_wild,mean_wild,sd_wild,n_captive,mean_captive,sd_captive))
}
colnames(summary_REH_sub) <- c("n_wild", "mean_wild", "sd_wild","n_captive", "mean_captive", "sd_captive")
rownames(summary_REH_sub) <- sublist
summary_REH_sub <- as.data.frame(summary_REH_sub)
meta_REH_sub_ready <- tibble::rownames_to_column(summary_REH_sub,"Author")
rownames(meta_REH_sub_ready) <- meta_REH_sub_ready$Author
sp_sorted <- c("TUTR","MOCH","BOGA","ELDA","CENI")
meta_REH_sub_ready <- meta_REH_sub_ready[sp_sorted,]
meta_q1_sub_tree.raw <- metacont(n_captive,mean_captive,sd_captive,n_wild,mean_wild,sd_wild,
                                 data = meta_REH_sub_ready,
                                 studlab = paste(Author),
                                 comb.fixed = FALSE,
                                 comb.random = TRUE,
                                 method.tau = "SJ",
                                 hakn = TRUE,
                                 prediction = TRUE,
                                 sm = "SMD")

#### Compositional differences: PERMANOVA #####
#Permutest and Adonis
#Pairwise (dis)similarity computation based on beta diversities
#R
permanovaR_results <- c()
permutestR_results <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  ##Filter the hierarchy
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  #Check whether both lists are identical
  identical(sort(colnames(final.table)),sort(as.character(hierarchy[,1])))
  ##Filter the metadata
  samples.kept <- colnames(final.table)
  metadata.filtered <- metadata[rownames(metadata) %in% samples.kept,]
  metadata.filtered1 <- tibble::rownames_to_column(metadata.filtered, "Code")
  row.names(metadata.filtered1) <- metadata.filtered1$Datafiles
  pairdis.R <- pair_dis(final.table,qvalue=0,hierarchy=hierarchy[,c(1,5)])
  u0n <- pairdis.R$L1_UqN
  u0n.dist <- as.dist(u0n)
  ps.disper.u0n.origin <- betadisper(u0n.dist, metadata.filtered1$Origin)
  permutestR <- permutest(ps.disper.u0n.origin, pairwise = TRUE)
  permutestR_results <- append(permutestR_results,permutestR[1])
  permanovaR <- adonis(u0n.dist ~ Origin, data =metadata.filtered1, permutations = 999)
  print(code)
  permanovaR_results <- append(permanovaR_results,permanovaR[1])
}
names(permanovaR_results) <- code.list
names(permutestR_results) <- code.list
saveRDS(permanovaR_results, "Results/RDS/permanova_R.RData")
saveRDS(permutestR_results, "Results/RDS/permutest_R.RData")

#REH
capwild.tree <- read.tree("Data/genustree.tre")
permanovaREH_results <- c()
permutestREH_results <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  rownames(final.table) <- gsub(" ","_",rownames(final.table))
  rownames(final.table) <- gsub("[","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("]","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("(","-",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub(")","-",rownames(final.table),fixed = TRUE)
  tree_filtered <- match_data(final.table,capwild.tree,output="tree")
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  samples.kept <- colnames(final.table)
  metadata.filtered <- metadata[rownames(metadata) %in% samples.kept,]
  metadata.filtered1 <- tibble::rownames_to_column(metadata.filtered, "Code")
  row.names(metadata.filtered1) <- metadata.filtered1$Datafiles
  pairdis.REH <- pair_dis(final.table,qvalue=1,hierarchy=hierarchy[,c(1,5)], tree=tree_filtered)
  u1n <- pairdis.REH$L1_UqN
  u1n.dist <- as.dist(u1n)
  ps.disper.u1n.origin <- betadisper(u1n.dist, metadata.filtered1$Origin)
  permutestREH <- permutest(ps.disper.u1n.origin, pairwise = TRUE)
  permutestREH_results <- append(permutestREH_results,permutestREH[1])
  permanovaREH <- adonis(u1n.dist ~ Origin, data =metadata.filtered1, permutations = 999)
  permanovaREH_results <- append(permanovaREH_results,permanovaREH[1])
}
names(permanovaREH_results) <- code.list
names(permutestREH_results) <- code.list
saveRDS(permanovaREH_results,paste("Results/permanova_REH.RData",sep=""))
saveRDS(permutestREH_results,paste("Results/permutest_REH.RData",sep=""))

#####Shapiro-Wilk’s method for normality test#####
#Preliminary test to check independent t-test assumptions: normality(Shapiro) & homogeneity in variances(F-test)
#Wild individuals
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  filtered.pecent <- tss(final.table)*100
  filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
  table.W <- tibble::rownames_to_column(filtered.pecent.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild.capt <- wild.capt[,-c(1)]
  #all.columns <- colnames(wild.capt)
  #taxa.all <- all.columns[-length(all.columns)]
  wild <- wild.capt[wild.capt$Origin == "Wild",]
  wild <- wild[,-ncol(wild)]
  wild <- wild[,colSums(wild) > 0]
  wild.taxa <- colnames(wild)
  shapiro_wild <- c()
  for (y in wild.taxa){
    res.shapiro.wild <- shapiro.test(wild[,y])
    shapiro_wild <- rbind(shapiro_wild,c(res.shapiro.wild[[1]],Pvalue=res.shapiro.wild[[2]]))
  }
  rownames(shapiro_wild) <- wild.taxa
  assign(paste("shapiro_wild_",code,sep=""),get("shapiro_wild"))
  write.table(get(paste("shapiro_wild_",code,sep="")),paste("Results/Shapiro/shapiro_wild_",code,".tsv",sep=""))
}
#Captive individuals
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  filtered.pecent <- tss(final.table)*100
  filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
  table.W <- tibble::rownames_to_column(filtered.pecent.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild.capt <- wild.capt[,-c(1)]
  captive <- wild.capt[wild.capt$Origin == "Captivity",]
  captive <- captive[,-ncol(captive)]
  captive <- captive[,colSums(captive) > 0]
  captive.taxa <- colnames(captive)
  shapiro_captive <- c()
  for (y in captive.taxa){
    res.shapiro.captive <- shapiro.test(captive[,y])
    shapiro_captive <- rbind(shapiro_captive,c(res.shapiro.captive[[1]],Pvalue=res.shapiro.captive[[2]]))
  }
  rownames(shapiro_captive) <- captive.taxa
  assign(paste("shapiro_captive_",code,sep=""),get("shapiro_captive"))
  write.table(get(paste("shapiro_captive_",code,sep="")),paste("Results/shapiro_captive_",code,".tsv",sep=""))
}

#### Unpaired Two-Samples Wilcoxon Test ####
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  filtered.pecent <- tss(final.table)*100
  filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
  table.W <- tibble::rownames_to_column(filtered.pecent.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild.capt <- wild.capt[,-c(1)]
  wild.capt.no <- wild.capt[,-ncol(wild.capt)]
  all.taxa <- colnames(wild.capt.no)
  Wilcox_result <- c()
  for (y in all.taxa){
    res.wilcox <- wilcox.test(wild.capt[,y] ~ Origin, data = wild.capt,
                              exact = FALSE, alternative = "less")
    Wilcox_result <- rbind(Wilcox_result,c(res.wilcox[[1]],Pvalue=res.wilcox[[3]]))
  }
  rownames(Wilcox_result) <- all.taxa
  Wilcox_result <- as.data.frame(Wilcox_result)
  Wilcox_result <- subset(Wilcox_result, Pvalue <= 0.01)
  assign(paste("Wilcox_result_",code,sep=""),get("Wilcox_result"))
  write.table(get(paste("Wilcox_result_",code,sep="")),paste("Results/WilcoxonTest/Wilcox_result_",code,".tsv",sep=""))
}

#Summary of relative values of each taxa in each group to interpret Wilcoxon test results
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  filtered.pecent <- tss(final.table)*100
  filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
  table.W <- tibble::rownames_to_column(filtered.pecent.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild.capt <- wild.capt[,-c(1)]
  wild <- wild.capt[wild.capt$Origin == "Wild",]
  wild <- wild[,-ncol(wild)]
  mean_wild <- sapply(wild,FUN=mean)
  mean_wild <- as.data.frame(mean_wild)
  sd_wild <- sapply(wild,FUN=sd)
  sd_wild <- as.data.frame(sd_wild)
  summary_wild <- merge(mean_wild,sd_wild, by = "row.names")
  captive <- wild.capt[wild.capt$Origin == "Captivity",]
  captive <- captive[,-ncol(captive)]
  mean_captive <- sapply(captive,FUN=mean)
  mean_captive <- as.data.frame(mean_captive)
  sd_captive <- sapply(captive,FUN=sd)
  sd_captive <- as.data.frame(sd_captive)
  summary_captive <- merge(mean_captive,sd_captive, by = "row.names")
  summary.rel <- merge(summary_wild,summary_captive, by = "Row.names")
  assign(paste("summary_rel_",code,sep=""),get("summary.rel"))
  write.table(get(paste("summary_rel_",code,sep="")),paste("Results/Relative_means/summary_rel_",code,".tsv",sep=""))
}

#### Shared taxa ####
shared.taxa <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  filtered.pecent <- tss(final.table)*100
  filtered.pecent.t <- as.data.frame(t(as.matrix(filtered.pecent)))
  table.W <- tibble::rownames_to_column(filtered.pecent.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild.capt <- wild.capt[,-c(1)]
  #all.columns <- colnames(wild.capt)
  #taxa.all <- all.columns[-length(all.columns)]
  wild <- wild.capt[wild.capt$Origin == "Wild",]
  wild <- wild[,-ncol(wild)]
  wild <- wild[,colSums(wild) > 0]
  wild.taxa <- colnames(wild)
  captive <- wild.capt[wild.capt$Origin == "Captivity",]
  captive <- captive[,-ncol(captive)]
  captive <- captive[,colSums(captive) > 0]
  captive.taxa <- colnames(captive)
  both.taxa <- intersect(wild.taxa,captive.taxa)
  wildonly.taxa <- wild.taxa[!wild.taxa %in% captive.taxa]
  captiveonly.taxa <- captive.taxa[! captive.taxa %in% wild.taxa]
  length(wildonly.taxa)
  length(both.taxa)
  length(captiveonly.taxa)
  shared <- cbind(length(wildonly.taxa),length(both.taxa),length(captiveonly.taxa))
  shared.taxa <- rbind(shared.taxa,shared)
}
colnames(shared.taxa) <- c("wild_only", "both", "captive_only")
rownames(shared.taxa) <- code.list
write.table(shared.taxa,"Results/shared_genus.tsv")

#### Beta diversity - incidence-based ####
summary_beta_all <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  final.table.inc <- to.incidence(final.table,hierarchy=hierarchy[,c(1,5)])
  divpart_q0 <- div_part(final.table.inc,qvalue=0)
  betadis_q0 <- beta_dis(beta=divpart_q0$Beta,qvalue=0,N=2)
  divpart_q1 <- div_part(final.table.inc,qvalue=1)
  betadis_q1 <- beta_dis(beta=divpart_q1$Beta,qvalue=1,N=2)
  summary_beta_all <- rbind(summary_beta_all,c(divpart_q0$Beta,divpart_q1$Beta,betadis_q0$CqN,betadis_q0$UqN,betadis_q1$CqN,betadis_q1$VqN,betadis_q1$SqN))
 }
colnames(summary_beta_all) <- c("beta_q0", "beta_q1","CqN_q0","UqN_q0","CqN_q1","VqN_q1","SqN_q1")
rownames(summary_beta_all) <- code.list
write.table(summary_beta_all,"Results/summary_beta_all.tsv")

#### Diversity analysis based on abundace at species level ####
summary_R_all <- c()
summary_REH_all <- c()
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  divqR <- div_test(final.table,qvalue=0,hierarchy=hierarchy[,c(1,5)])
  divqR_table <- divqR$data
  n_all <- nrow(divqR_table)
  mean_all <- mean(divqR_table[,c("Value")])
  sd_all <- sd(divqR_table[,c("Value")])
  summary_R_all <- rbind(summary_R_all,c(n_all,mean_all,sd_all))
}
colnames(summary_R_all) <- c("n_R", "mean_R", "sd_R")
rownames(summary_R_all) <- code.list
for (code in code.list){
  final.table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  rownames(final.table) <- gsub(" ","_",rownames(final.table))
  rownames(final.table) <- gsub("[","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("]","",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub("(","-",rownames(final.table),fixed = TRUE)
  rownames(final.table) <- gsub(")","-",rownames(final.table),fixed = TRUE)
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(final.table)),]
  tree_filtered <- match_data(final.table,capwild.tree,output="tree")
  divREH <- div_test(final.table,qvalue=1,hierarchy=hierarchy[,c(1,5)], tree=tree_filtered)
  divREH_table <- divREH$data
  n_all <- nrow(divREH_table)
  mean_all <- mean(divREH_table[,c("Value")])
  sd_all <- sd(divREH_table[,c("Value")])
  summary_REH_all <- rbind(summary_REH_all,c(n_all,mean_all,sd_all))
}
colnames(summary_REH_all) <- c("n_REH", "mean_REH", "sd_REH")
rownames(summary_REH_all) <- code.list
summary_all <- merge(summary_R_all,summary_REH_all, by= "row.names")
summary_all <- tibble::column_to_rownames(summary_all, "Row.names")
write.table(summary_all, "Results/Mean_diversity_species.tsv")

#### Species level analysis ####
#Merging all tables in one
#Load files at genus level
files.g = list.files(path="Tables",pattern="countfiltered_*")
#names.g = gsub("countfiltered_","",files.g)
#names.g = gsub(".tsv","",names.g)
filelist.g = lapply(paste("Tables",files.g,sep="/"), function(x)read.table(x, header=T))
rownameslist = lapply(filelist.g, function(x) rownames(x))
for (l in c(1:length(filelist.g))){
  filelist.g[[l]]$taxa <- rownames(filelist.g[[l]])
}
#Note some erroneous files can prevent all files from loading (remove erroneous ones)
#Transform list to matrix
filelist.g %>% purrr::reduce(full_join, by = "taxa") -> count.table.all.g
rownames(count.table.all.g) <- count.table.all.g[,"taxa"]
count.table.all.g <- count.table.all.g[,colnames(count.table.all.g) != "taxa"]
#Name columns
#colnames(count.table.all.g) <- names.g
#Replace NAs by 0s
count.table.all.g[is.na(count.table.all.g)] <- 0
write.table(count.table.all.g, "Tables/count_Genus_all.tsv")

#### Diversity analysis at species level ####
count.table.all.g <- read.table("Tables/count_Genus_all.tsv")
rownames(count.table.all.g) <- gsub(" ","_",rownames(count.table.all.g))
rownames(count.table.all.g) <- gsub("[","",rownames(count.table.all.g),fixed = TRUE)
rownames(count.table.all.g) <- gsub("]","",rownames(count.table.all.g),fixed = TRUE)
rownames(count.table.all.g) <- gsub("(","-",rownames(count.table.all.g),fixed = TRUE)
rownames(count.table.all.g) <- gsub(")","-",rownames(count.table.all.g),fixed = TRUE)
hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(count.table.all.g)),]

##Diversity comparison
#R
divR.all <- div_test(count.table.all.g,qvalue=0,hierarchy=hierarchy[,c(1,3)])
#REH
capwild.tree <- read.tree("Data/genustree.tre")
tree_filtered <- match_data(count.table.all.g,capwild.tree,output="tree")
divREH.all <- div_test(count.table.all.g,qvalue=1,hierarchy=hierarchy[,c(1,3)],tree = tree_filtered)
saveRDS(divREH.all,"Results/RDS/divREH.all.tree.RData")

##Pairwise (dis)similarity computation based on beta diversities
#R
pairdis.R.all <- pair_dis(count.table.all.g,qvalue=0,hierarchy=hierarchy[,c(1,3)])
saveRDS(pairdis.R.all,"Results/RDS/pairdis_R_all.RData")
#Plot it (remove outlier D12719)
pairdis.R.all <- readRDS("Results/RDS/pairdis_R_all.RData")
u0n <- pairdis.R.all$L1_UqN
u0n <- u0n[rownames(u0n) != "D12719",colnames(u0n) != "D12719"]
dis_nmds(u0n,hierarchy=hierarchy[hierarchy$Samples != "D12719",c(1,19)], centroids=TRUE, runs=200)
#REH
pairdis.REH.all <- pair_dis(count.table.all.g,qvalue=1,hierarchy=hierarchy[,c(1,3)], tree = tree_filtered)
saveRDS(pairdis.REH.all,"Results/RDS/pairdis_REH_all.RData")
#Plot it (remove outlier D12719)
pairdis.REH.all <- readRDS("Results/RDS/pairdis_REH_all.RData")
u1n.tree <- pairdis.REH.all$L1_UqN
u1n.tree <- u1n.tree[rownames(u1n.tree) != "D12719",colnames(u1n.tree) != "D12719"]
dis_nmds(u1n.tree,hierarchy=hierarchy[hierarchy$Samples != "D12719",c(1,19)], centroids=TRUE, runs=200)

#Permanova
##Filter the metadata
samples.kept <- colnames(count.table.all.g)
metadata.filtered <- metadata[rownames(metadata) %in% samples.kept,]
#R
pairdis.R.all <- readRDS("Results/RDS/pairdis_R_all.RData")
u0n <- pairdis.R.all$L1_UqN
u0n.dist <- as.dist(u0n)
ps.disper.u0n.species <- betadisper(u0n.dist, metadata.filtered$Species)
permutest(ps.disper.u0n.species, pairwise = TRUE)
adonis(u0n.dist ~ Species, data =metadata.filtered, permutations = 999)
#REH
pairdis.REH.all <- readRDS("Results/RDS/pairdis.REH.all")
u1n.tree <- pairdis.REH.all$L1_UqN
u1n.dist.tree <- as.dist(u1n.tree)
ps.disper.u1n.tree <- betadisper(u1n.dist.tree, metadata.filtered$Species)
permutest(ps.disper.u1n.tree, pairwise = TRUE)
adonis(u1n.dist.tree ~ Species, data = metadata.filtered, permutations = 999)

#Correlation: diversity vs sequencing depth
diversity1.tree <- divq1.all.tree$data
seq.depth <- colSums(count.table.all.g)
seq.depth1 <- as.data.frame(seq.depth)
seq.depth1 <- tibble::rownames_to_column(seq.depth1, "Sample")
table <- merge(diversity1.tree,seq.depth1, by="Sample")
div <- table$Value
seq <- table$seq.depth
smd <- as.numeric(as.character(merged$smd_value))
cor.test(div,seq)
plot(div,seq)
abline(lm(seq ~ div, data = table), col = "red")





#### Have we included this????
#### Sequencing depth by species ####
Wilcox_result <- c()
for (code in code.list){
  print(code)
  table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(table)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  #filtered.pecent <- tss(final.table)*100
  table.t <- as.data.frame(t(as.matrix(table)))
  table.W <- tibble::rownames_to_column(table.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild.capt <- tibble::column_to_rownames(wild.capt, "Samples")
  wild.capt <- wild.capt[,-ncol(wild.capt)]
  sequencing <- rowSums(wild.capt)
  sequencing <- as.data.frame(sequencing)
  sequencing <- tibble::rownames_to_column(sequencing, "Samples")
  sequencing <- merge(sequencing,hierarchy_species,by="Samples")
  sequencing <- sequencing[,-c(1)]
  colnames(sequencing) <- c("Sequencing_depth", "Origin")
  res.wilcox <- wilcox.test(Sequencing_depth ~ Origin, data = sequencing,
                            exact = FALSE, alternative = "less")
  Wilcox_result <- rbind(Wilcox_result,c(res.wilcox[[1]],Pvalue=res.wilcox[[3]]))
}
rownames(Wilcox_result) <- code.list


#### Sequencing depth means ####
summary.seq <- c()
for (code in code.list){
  print(code)
  table <- read.table(paste("Tables/countfiltered_",code,".tsv",sep=""))
  hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(table)),]
  hierarchy_species <- hierarchy[,c(1,5)]
  #filtered.pecent <- tss(final.table)*100
  table.t <- as.data.frame(t(as.matrix(table)))
  table.W <- tibble::rownames_to_column(table.t, "Samples")
  wild.capt <- merge(table.W,hierarchy_species,by="Samples")
  wild <- wild.capt[wild.capt$Origin == "Wild",]
  wild <- tibble::rownames_to_column(wild, "N")
  wild <- tibble::column_to_rownames(wild, "Samples")
  wild <- wild[,-c(1)]
  wild <- wild[,-ncol(wild)]
  wild.seq <- rowSums(wild)
  wild.seq <- as.data.frame(wild.seq)
  mean_wild <- sapply(wild.seq,FUN=mean)
  sd_wild <- sapply(wild.seq,FUN=sd)
  Captive <- wild.capt[wild.capt$Origin == "Captivity",]
  Captive <- tibble::rownames_to_column(Captive, "N")
  Captive <- tibble::column_to_rownames(Captive, "Samples")
  Captive <- Captive[,-c(1)]
  Captive <- Captive[,-ncol(Captive)]
  Captive.seq <- rowSums(Captive)
  Captive.seq <- as.data.frame(Captive.seq)
  mean_captive <- sapply(Captive.seq,FUN=mean)
  sd_captive <- sapply(Captive.seq,FUN=sd)
  summary.seq <- rbind(summary.seq, c(mean_wild,sd_wild,mean_captive,sd_captive))
}
colnames(summary.seq) <- c("Wild.Mean", "Wild.SD", "Captive.Mean", "Captive.SD")
rownames(summary.seq) <- code.list

##### Did captivity enrich some taxa?####
incid_result_table <- c()
incid_ratio_table <- c()
count.table.all <- read.table(paste("Tables/count_all.tsv",sep = ""))
hierarchy <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(count.table.all)),]
for (code in code.list){
  filtered.table.w <- count.table.all[,hierarchy[(hierarchy$Species == code) & (hierarchy$Origin == "Wild"),"Samples"]]
  wild <- to.incidence(filtered.table.w)
  filtered.table.c <- count.table.all[,hierarchy[(hierarchy$Species == code) & (hierarchy$Origin == "Captivity"),"Samples"]]
  captive <- to.incidence(filtered.table.c)
  incid_result <- cbind(wild,captive)
  incid_result1 <- as.data.frame(cbind(wild=incid_result[,1]*100/ncol(filtered.table.w),captive=incid_result[,2]*100/ncol(filtered.table.c)))
  #Ratio
  incid_ratio_vector <- incid_result1$captive / incid_result1$wild
  incid_ratio_table <- cbind(incid_ratio_table,incid_ratio_vector)
  #Data
  colnames(incid_result1) <- c(paste(code,"wild",sep="_"),paste(code,"captive",sep="_"))
  if (is.null(nrow(incid_result_table))){
    incid_result_table <- incid_result1
    }else{
    incid_result_table <- cbind(incid_result_table,incid_result1)
  }
}
rownames(incid_ratio_table) <- rownames(incid_result_table)
colnames(incid_ratio_table) <- code.list
incid_ratio_table2 <- incid_ratio_table
incid_ratio_table2[incid_ratio_table2 > 1] <- 1
incid_ratio_table2[incid_ratio_table2 < 1] <- 0
captive.enriched <- rowSums(incid_ratio_table2,na.rm=TRUE) / rowSums(!is.na(incid_ratio_table2)) * 100
captive.enriched[captive.enriched == 100]
